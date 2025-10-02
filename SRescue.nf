// Default parameter values
params.samplelist = "/home/m.wyczalkowski/Projects/Adhoc/2025/SRescue/SRescue/dat/samplelist.dat"
params.survivor = "/opt/SURVIVOR/Debug/SURVIVOR"
params.cutesv = "/opt/conda/envs/env/bin/cuteSV"
params.bgzip = "bgzip"
params.tabix = "tabix"
params.perl = "/usr/bin/perl"
params.ref = "ref.fa"
params.samtools = "samtools"
params.config = null

// Note: must not include the "/rdcw/fs2/home1/Active" prefix to home 
params.script_polish_survivor = "/home/m.wyczalkowski/Projects/Adhoc/2025/SRescue/SRescue/script/01_polish_survivor.pl"
params.script_get_only_svs = "/home/m.wyczalkowski/Projects/Adhoc/2025/SRescue/SRescue/script/02_get.only.pl"
params.script_vcf4fc = "/home/m.wyczalkowski/Projects/Adhoc/2025/SRescue/SRescue/script/03_makeVCF4forceCalling.pl"
params.script_finalfilter = "/home/m.wyczalkowski/Projects/Adhoc/2025/SRescue/SRescue/script/04_get.finalfilter_cuteSV_ins.pl"

params.reference = "/storage1/fs1/dinglab/Active/Projects/PECGS/ref_genome/GRCh38.d1.vd1.fa"

def lrvcf_fn="/storage1/fs1/dinglab/Active/Projects/yuweiz/projects/HTAN/ccRCC/Longread/t2t/01_compare2sr/00_start_allSample/SRescue/test_data/LR.chr2.vcf"
def srvcf_fn="/storage1/fs1/dinglab/Active/Projects/yuweiz/projects/HTAN/ccRCC/Longread/t2t/01_compare2sr/00_start_allSample/SRescue/test_data/SR.chr2.vcf"
def norm_fn="/storage1/fs1/dinglab/Active/Projects/yuweiz/projects/HTAN/ccRCC/Longread/t2t/01_compare2sr/00_start_allSample/SRescue/test_data/LR.normal.chr2.bam"
def norm_bai="/storage1/fs1/dinglab/Active/Projects/yuweiz/projects/HTAN/ccRCC/Longread/t2t/01_compare2sr/00_start_allSample/SRescue/test_data/LR.normal.chr2.bam.bai"
def tum_fn="/storage1/fs1/dinglab/Active/Projects/yuweiz/projects/HTAN/ccRCC/Longread/t2t/01_compare2sr/00_start_allSample/SRescue/test_data/LR.tumor.chr2.bam"
def tum_bai="/storage1/fs1/dinglab/Active/Projects/yuweiz/projects/HTAN/ccRCC/Longread/t2t/01_compare2sr/00_start_allSample/SRescue/test_data/LR.tumor.chr2.bam.bai"

// Process 1: Merge LR and SR VCFs using SURVIVOR
// requires SURVIVOR
process merge_lr_sr {
    label 'survivor'
    input:
    path(lrvcf)
    path(srvcf)

    output:
    path("sr2lr.vcf")

    script:
    """
    echo "${srvcf}" > lr_sr.vcflist
    echo "${lrvcf}" >> lr_sr.vcflist
    ${params.survivor} merge lr_sr.vcflist 1000 1 0 0 0 50 sr2lr.vcf
    """
}

// Process 2: Polish SURVIVOR output
// Requires PERL
process polish_survivor {
    label 'perl'
    input:
    path(merged_vcf)
    path(perl_script)

    output:
    path("sr2lr.polished.bedpe")
    path("sr2lr.polished.vcf")

    script:
    """
    ${params.perl} ${perl_script} ${merged_vcf} sr2lr.polished.bedpe sr2lr.polished.vcf
    """
}

// Process 3: Get only SVs
process get_only_svs {
    label 'perl'
    input:
    path(polished_bedpe)
    path(polished_vcf)
    path(perl_script)
    
    output:
//    path("sronly.tsv")
//    path("lronly.tsv")
//    path("shared.tsv")
    path("sronly.vcf")
    path("lronly.vcf")
    path("shared.vcf")
    
    script:
    """
    ${params.perl} ${perl_script} ${polished_bedpe} ${polished_vcf} ./
    """
}

// Process 4: Make VCF for force calling
process make_vcf_for_force_calling {
    label 'perl'
    input:
    path(sronly_vcf)
    path(perl_script)
    
    output:
    path("sronly4fc.vcf")
    
    script:
    """
    ${params.perl} ${perl_script} ${sronly_vcf}
    """
}


// Process 5: Run cuteSV 
process run_cutesv_normal {
    label 'cutesv'
    input:
    path bam
    path bai
    path reference
    path sronly4fc_vcf

    output:
    path("cutesv.call.vcf")
    
    script:
    """
    ${params.cutesv} ${bam} ${reference} cutesv.call.vcf ./ --min_mapq 10 -Ivcf ${sronly4fc_vcf} -t 2 -L -1 --report_readid
    """
}

// Process 6: run cuteSV on tumor
process run_cutesv_tumor {
    label 'cutesv'
    input:
    path bam
    path bai
    path reference
    path sronly4fc_vcf

    output:
    path("cutesv.call.vcf")
    
    script:
    """
    ${params.cutesv} ${bam} ${reference} cutesv.call.vcf ./ --min_mapq 10 -Ivcf ${sronly4fc_vcf} -t 2 -L -1 --report_readid
    """
}


// Process 7: Get final filter cuteSV insertions
process get_final_filter {
    label 'perl'
    input:
    path "cutesv_tumor.call.vcf"
    path "cutesv_normal.call.vcf"
    path(perl_script)
    
    output:
    path("finalfilter_cutesv.vcf")
    path("supread.loci.bed")
    path("supread.tsv")
    
    script:
    """
    ${params.perl} ${perl_script} ./
    """
}


// Main workflow
workflow {
    // Process 1
    sr2lr = merge_lr_sr(
        channel.fromPath(lrvcf_fn),
        channel.fromPath(srvcf_fn)
    )

    // TODO: should not be passing the scripts as a channel from here, should be in the process itself
    // Process 2
    (s2_bedpe, s2_vcf) = polish_survivor(sr2lr, params.script_polish_survivor)

    // Process 3
    (sro, lro, sha) = get_only_svs(s2_bedpe, s2_vcf, params.script_get_only_svs)

    // Process 4
    sro4fc = make_vcf_for_force_calling(sro, params.script_vcf4fc)

    // Process 5 - cuteSV on normal
    normal_cutesv = run_cutesv_normal(norm_fn, norm_bai, params.reference, sro4fc)

    // Process 6 - cuteSV on tumor
    tumor_cutesv = run_cutesv_tumor(tum_fn, tum_bai, params.reference, sro4fc)

    // Process 7
    (vcf, bed, tsv) = get_final_filter(normal_cutesv, tumor_cutesv, params.script_finalfilter)
}


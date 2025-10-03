// Default parameter values
params.samplelist = "/home/m.wyczalkowski/Projects/Adhoc/2025/SRescue/SRescue/dat/samplelist.dat"
params.survivor = "/opt/SURVIVOR/Debug/SURVIVOR"
params.cutesv = "/opt/conda/envs/env/bin/cuteSV"
params.bgzip = "/usr/bin/bgzip"
params.tabix = "/usr/bin/tabix"
//params.perl = "/usr/bin/perl"   
params.samtools = "/usr/local/bin/samtools"
params.config = null

// Note: must not include the "/rdcw/fs2/home1/Active" prefix to home 


/*
// Scripts 
// should not be passing paths to scripts on the command line.  A preferred approach seems to be binary modules
// documentation: https://www.nextflow.io/docs/latest/module.html#module-binaries
// example: https://github.com/davidmasp/hello-modbins/tree/master

// may need to add these to the bin directory
// Also: https://sateeshperi.github.io/nextflow_varcal/nextflow/nextflow_modules
// also here: https://github.com/nextflow-io/nextflow/discussions/4651
*/

// all perl scripts are added to the bin directory
// https://www.nextflow.io/docs/latest/sharing.html#the-bin-directory

params.reference = "/storage1/fs1/dinglab/Active/Projects/PECGS/ref_genome/GRCh38.d1.vd1.fa"

# TODO: how to nicely specify secondary / index file?
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

    output:
    path("sr2lr.polished.bedpe")
    path("sr2lr.polished.vcf")

    script:
    """
    01_polish_survivor.pl ${merged_vcf} sr2lr.polished.bedpe sr2lr.polished.vcf
    """
}

// Process 3: Get only SVs
process get_only_svs {
    label 'perl'
    input:
    path(polished_bedpe)
    path(polished_vcf)
    
    output:
//    path("sronly.tsv")
//    path("lronly.tsv")
//    path("shared.tsv")
    path("sronly.vcf")
    path("lronly.vcf")
    path("shared.vcf")
    
    script:
    """
    02_get.only.pl ${polished_bedpe} ${polished_vcf} ./
    """
}

// Process 4: Make VCF for force calling
process make_vcf_for_force_calling {
    label 'perl'
    input:
    path(sronly_vcf)
    
    output:
    path("sronly4fc.vcf")
    
    script:
    """
    03_makeVCF4forceCalling.pl ${sronly_vcf}
    """
}


// Process 5: Run cuteSV 
process run_cutesv_normal {
    label 'cutesv'
    input:
    path normal_bam
    path normal_bai
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
    path tumor_bam
    path tumor_bai
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
    
    output:
    path("finalfilter_cutesv.vcf")
    path("supread.loci.bed")
    path("supread.tsv")
    
    script:
    """
    04_get.finalfilter_cuteSV_ins.pl ./
    """
}

// Process 8: Process supporting reads
process process_supporting_reads {
    label 'samtools'
    
    input:
    path(supread_loci_bed)
    path(tumorbam)
    path("finalfilter_cutesv.vcf")
    path(supread_tsv)
    
    output:
    path("finalfilter_cutesv.vcf")  // same filename but different file than the input
    
    script:
    """
    ${params.samtools} depth -b ${supread_loci_bed} ${tumorbam} > supread.depth.tsv
    ${params.samtools} view -N ${supread_tsv} ${tumorbam} > supread.bam.tsv
    04_2_correct.rb.pl ./
    rm finalfilter_cutesv2.vcf && mv finalfilter_cutesv3.vcf finalfilter_cutesv.vcf
    """

//    04_2_correct.rb.pl 
//        reads - open(supread,"$workdir/supread.bam.tsv");
//        reads - open(file,"$workdir/finalfilter_cutesv.vcf");
//        reads - open(file,"$workdir/supread.depth.tsv");
//        writes - open(res,">$workdir/finalfilter_cutesv2.vcf");
//        writes - open(res,">$workdir/finalfilter_cutesv3.vcf");
}

// Process 9: Polish all VCFs
process polish_all_vcfs {
    label 'perl'
    
    input:
    path(sronly_vcf)
    path(lronly_vcf)
    path(shared_vcf)
    
    output:
    path("sronly.polished.bedpe")
    path("sronly.polished.vcf")
    path("lronly.polished.bedpe")
    path("lronly.polished.vcf")
    path("shared.polished.bedpe")
    path("shared.polished.vcf")
    
    script:
    """
    01_polish_survivor.pl ${sronly_vcf} sronly.polished.bedpe sronly.polished.vcf
    01_polish_survivor.pl ${lronly_vcf} lronly.polished.bedpe lronly.polished.vcf
    01_polish_survivor.pl ${shared_vcf} shared.polished.bedpe shared.polished.vcf
    """
}

// Process 10: Merge final VCFs
process merge_final_vcfs {
    label 'survivor'
    
    input:
    path(lrvcf)
    path(final_filter_vcf)

    output:
    path("final.sr2lr.sv.vcf")
    
    script:
    """
    echo "${lrvcf}" > final.vcflist
    echo "${final_filter_vcf}" >> final.vcflist
    ${params.survivor} merge final.vcflist 1000 1 0 0 0 50 final.sr2lr.sv.vcf
    """
}

// Process 11: Final polish and compress
process final_polish_and_compress {
    label 'perl'
    
    input:
    path(final_sr2lr_sv_vcf)
    path(shared_polished_bedpe)
    
    output:
    path("final.sr2lr.polished.vcf.gz")
    path("final.sr2lr.polished.vcf.gz.tbi")
    
    script:
    """
    05_final_polish.pl ${final_sr2lr_sv_vcf} ${shared_polished_bedpe} final.sr2lr.sv.bedpe final.sr2lr.polished.vcf
    ${params.bgzip} -c final.sr2lr.polished.vcf > final.sr2lr.polished.vcf.gz
    ${params.tabix} final.sr2lr.polished.vcf.gz
    """
}

// Main workflow

// TODO: stage outputs nicely
// final.sr2lr.polished.vcf.gz 
// finalfilter_cutesv.vcf
// shared.polished.vcf
// sronly.polished.vcf
// lronly.polished.vcf

workflow {
    // Process 1
    sr2lr = merge_lr_sr(lrvcf_fn, srvcf_fn)

    // TODO: should not be passing the scripts as a channel from here, should be in the process itself
    // Process 2
    (s2_bedpe, s2_vcf) = polish_survivor(sr2lr)

    // Process 3
    (sro, lro, sha) = get_only_svs(s2_bedpe, s2_vcf)

    // Process 4
    sro4fc = make_vcf_for_force_calling(sro)

    // Process 5 - cuteSV on normal
    normal_cutesv = run_cutesv_normal(norm_fn, norm_bai, params.reference, sro4fc)

    // Process 6 - cuteSV on tumor
    tumor_cutesv = run_cutesv_tumor(tum_fn, tum_bai, params.reference, sro4fc)

    // Process 7
    (vcf, bed, tsv) = get_final_filter(normal_cutesv, tumor_cutesv)

    // Process 8
    ffc_vcf = process_supporting_reads(bed, tum_fn, vcf, tsv)

    // Process 9
    (srp_bpe, srp_vcf, lrp_bpe, lrp_vcf, shp_bpe, shp_vcf) = polish_all_vcfs(sro, lro, sha)

    // Process 10
    sr2lr_sv_vcf = merge_final_vcfs(lrvcf_fn, ffc_vcf)

    // Process 11
    final_polish_and_compress(sr2lr_sv_vcf, shp_bpe)
}


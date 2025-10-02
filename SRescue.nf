// Default parameter values
params.samplelist = "/rdcw/fs2/home1/Active/home/m.wyczalkowski/Projects/Adhoc/2025/SRescue/SRescue/dat/samplelist.dat"
params.outdir = "results"
params.survivor = "/opt/SURVIVOR/Debug/SURVIVOR"
params.cutesv = "cuteSV"
params.bgzip = "bgzip"
params.tabix = "tabix"
params.perl = "/usr/bin/perl"
//params.script_path = "./script"
params.script_path = "/rdcw/fs2/home1/Active/home/m.wyczalkowski/Projects/Adhoc/2025/SRescue/SRescue/script"
params.ref = "ref.fa"
params.samtools = "samtools"
params.config = null

def lrvcf_fn="/storage1/fs1/dinglab/Active/Projects/yuweiz/projects/HTAN/ccRCC/Longread/t2t/01_compare2sr/00_start_allSample/SRescue/test_data/LR.chr2.vcf"
def srvcf_fn="/storage1/fs1/dinglab/Active/Projects/yuweiz/projects/HTAN/ccRCC/Longread/t2t/01_compare2sr/00_start_allSample/SRescue/test_data/SR.chr2.vcf"
def norm_fn="/storage1/fs1/dinglab/Active/Projects/yuweiz/projects/HTAN/ccRCC/Longread/t2t/01_compare2sr/00_start_allSample/SRescue/test_data/LR.normal.chr2.bam"
def tum_fn="/storage1/fs1/dinglab/Active/Projects/yuweiz/projects/HTAN/ccRCC/Longread/t2t/01_compare2sr/00_start_allSample/SRescue/test_data/LR.tumor.chr2.bam"

// Process 1: Merge LR and SR VCFs using SURVIVOR
// requires SURVIVOR
process merge_lr_sr {
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
    input:
    path(merged_vcf)

    output:
    path("sr2lr.polished.bedpe")
    path("sr2lr.polished.vcf")

    script:
    """
    ${params.perl} ${params.script_path}/01_polish_survivor.pl ${merged_vcf} sr2lr.polished.bedpe sr2lr.polished.vcf
    """
}

// Process 3: Get only SVs
process get_only_svs {
    input:
    path(polished_bedpe)
    path(polished_vcf)
    
    output:
    path("sronly.tsv")
    path("lronly.tsv")
    path("shared.tsv")
    path("sronly.vcf")
    path("lronly.vcf")
    path("shared.vcf")
    
    script:
    """
    ${params.perl} ${params.script_path}/02_get.only.pl ${polished_bedpe} ${polished_vcf} ${params.outdir}
    """
}


// Main workflow
workflow {
    // Process 1
    sr2lr = merge_lr_sr(
        channel.fromPath(lrvcf_fn),
        channel.fromPath(srvcf_fn)
    )

    // Process 2
    (s2_bedpe, s2_vcf) = polish_survivor(sr2lr)

    // Process 3
    r3 = get_only_svs(s2_bedpe, s2_vcf)
}


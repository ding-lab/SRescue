#!/usr/bin/env nextflow

/*
* SV Validation Pipeline
* This pipeline validates SVs detected by short-read WGS using Tumor-Normal long-read WGS.
* Author: Original by ywzhang0713@gmail.com, converted to Nextflow
*/

// Print help message if no parameters are provided
if (params.help) {
    helpMessage()
    exit 0
}

// Default parameter values
params.samplelist = null
params.outdir = "results"
params.survivor = "SURVIVOR"
params.cutesv = "cuteSV"
params.bgzip = "bgzip"
params.tabix = "tabix"
params.script_path = "${projectDir}/bin"
params.ref = "ref.fa"
params.samtools = "samtools"
params.config = null

// Function to print help message
def helpMessage() {
    log.info"""
    =========================================
    SV Validation Pipeline
    =========================================
    
    This pipeline was designed for validating the SVs detected by only short-read WGS using Tumor-Normal long-read WGS.
    Original Author: ywzhang0713@gmail.com Yuwei ZHANG
    
    Usage:
    nextflow run main.nf --config <config_file>
    
    Parameters:
    --config            Config file path (if provided, other parameters from this file will be used)
    --samplelist        Path to sample list file
    --outdir            Output directory (default: results)
    --survivor          Path to SURVIVOR executable (default: SURVIVOR)
    --cutesv            Path to cuteSV executable (default: cuteSV)
    --bgzip             Path to bgzip executable (default: bgzip)
    --tabix             Path to tabix executable (default: tabix)
    --script_path       Path to helper scripts (default: bin directory in workflow)
    --ref               Reference genome path (default: ref.fa)
    --samtools          Path to samtools executable (default: samtools)
    --help              Print this help message
    """.stripIndent()
}

// Load configuration from file if provided
if (params.config) {
    def configFile = file(params.config)
    if (!configFile.exists()) {
        exit 1, "Config file not found: ${params.config}"
    }
    
    configFile.eachLine { line ->
        if (!line.startsWith('#') && line.trim()) {
            def fields = line.split('\t')
            if (fields.size() >= 2) {
                switch(fields[0]) {
                    case 'SAMPLE':
                        params.samplelist = fields[1]
                        break
                    case 'OUTPATH':
                        params.outdir = fields[1]
                        break
                    case 'SURVIVOR':
                        params.survivor = fields[1]
                        break
                    case 'CUTESV':
                        params.cutesv = fields[1]
                        break
                    case 'BGZIP':
                        params.bgzip = fields[1]
                        break
                    case 'TABIX':
                        params.tabix = fields[1]
                        break
                    case 'PIPE_SCRIPT_PATH':
                        params.script_path = fields[1]
                        break
                    case 'REF':
                        params.ref = fields[1]
                        break
                    case 'SAMTOOLS':
                        params.samtools = fields[1]
                        break
                }
            }
        }
    }
}

// Check for required parameters
if (!params.samplelist) {
    exit 1, "Sample list file not specified. Provide it with --samplelist or in the config file."
}

// Create sample channels from the sample list file
Channel
    .fromPath(params.samplelist)
    .splitCsv(header: true, sep: '\t')
    .map { row -> 
        def sid = row[0]
        def lrvcf = file(row[1])
        def srvcf = file(row[2])
        def normalbam = file(row[3])
        def tumorbam = file(row[4])
        
        if (lrvcf.exists() && srvcf.exists() && normalbam.exists() && tumorbam.exists()) {
            return tuple(sid, lrvcf, srvcf, normalbam, tumorbam)
        } else {
            log.error "Error: check the availability of the files for sample ${sid}"
            return null
        }
    }
    .filter { it != null }
    .set { samples_ch }

// Process 1: Merge LR and SR VCFs using SURVIVOR
process merge_lr_sr {
    tag "${sid}"
    publishDir "${params.outdir}/${sid}", mode: 'copy'
    
    input:
    tuple val(sid), path(lrvcf), path(srvcf), path(normalbam), path(tumorbam)
    
    output:
    tuple val(sid), path("sr2lr.vcf"), path(lrvcf), path(srvcf), path(normalbam), path(tumorbam), emit: merged_output
    
    script:
    """
    echo "${srvcf}" > lr_sr.vcflist
    echo "${lrvcf}" >> lr_sr.vcflist
    ${params.survivor} merge lr_sr.vcflist 1000 1 0 0 0 50 sr2lr.vcf
    """
}

// Process 2: Polish SURVIVOR output
process polish_survivor {
    tag "${sid}"
    publishDir "${params.outdir}/${sid}", mode: 'copy'
    
    input:
    tuple val(sid), path(merged_vcf), path(lrvcf), path(srvcf), path(normalbam), path(tumorbam)
    
    output:
    tuple val(sid), path("sr2lr.polished.bedpe"), path("sr2lr.polished.vcf"), path(lrvcf), path(srvcf), path(normalbam), path(tumorbam), emit: polished_output
    
    script:
    """
    perl ${params.script_path}/01_polish_survivor.pl ${merged_vcf} sr2lr.polished.bedpe sr2lr.polished.vcf
    """
}

// Process 3: Get only SVs
process get_only_svs {
    tag "${sid}"
    publishDir "${params.outdir}/${sid}", mode: 'copy'
    
    input:
    tuple val(sid), path(polished_bedpe), path(polished_vcf), path(lrvcf), path(srvcf), path(normalbam), path(tumorbam)
    
    output:
    tuple val(sid), path("sronly.vcf"), path("lronly.vcf"), path("shared.vcf"), path(lrvcf), path(normalbam), path(tumorbam), emit: only_svs_output
    
    script:
    """
    perl ${params.script_path}/02_get.only.pl ${polished_bedpe} ./
    """
}

// Process 4: Make VCF for force calling
process make_vcf_for_force_calling {
    tag "${sid}"
    publishDir "${params.outdir}/${sid}", mode: 'copy'
    
    input:
    tuple val(sid), path(sronly_vcf), path(lronly_vcf), path(shared_vcf), path(lrvcf), path(normalbam), path(tumorbam)
    
    output:
    tuple val(sid), path("sronly4fc.vcf"), path(sronly_vcf), path(lronly_vcf), path(shared_vcf), path(lrvcf), path(normalbam), path(tumorbam), emit: vcf_for_fc_output
    
    script:
    """
    perl ${params.script_path}/03_makeVCF4forceCalling.pl ${sronly_vcf}
    """
}

// Process 5: Run cuteSV on normal sample
process run_cutesv_normal {
    tag "${sid}"
    publishDir "${params.outdir}/${sid}", mode: 'copy'
    
    input:
    tuple val(sid), path(sronly4fc_vcf), path(sronly_vcf), path(lronly_vcf), path(shared_vcf), path(lrvcf), path(normalbam), path(tumorbam)
    
    output:
    tuple val(sid), path("cutesv_normal.call.vcf"), path(sronly4fc_vcf), path(sronly_vcf), path(lronly_vcf), path(shared_vcf), path(lrvcf), path(normalbam), path(tumorbam), emit: cutesv_normal_output
    
    script:
    """
    if [ -e signatures ]; then
        rm -r signatures
    fi
    ${params.cutesv} ${normalbam} ${params.ref} cutesv_normal.call.vcf ./ --min_mapq 10 -Ivcf ${sronly4fc_vcf} -t 2 -L -1 --report_readid
    """
}

// Process 6: Run cuteSV on tumor sample
process run_cutesv_tumor {
    tag "${sid}"
    publishDir "${params.outdir}/${sid}", mode: 'copy'
    
    input:
    tuple val(sid), path(cutesv_normal_vcf), path(sronly4fc_vcf), path(sronly_vcf), path(lronly_vcf), path(shared_vcf), path(lrvcf), path(normalbam), path(tumorbam) 
    
    output:
    tuple val(sid), path("cutesv_tumor.call.vcf"), path(cutesv_normal_vcf), path(sronly_vcf), path(lronly_vcf), path(shared_vcf), path(lrvcf), path(normalbam), path(tumorbam), emit: cutesv_tumor_output
    
    script:
    """
    if [ -e signatures ]; then
        rm -r signatures
    fi
    ${params.cutesv} ${tumorbam} ${params.ref} cutesv_tumor.call.vcf ./ --min_mapq 10 -Ivcf ${sronly4fc_vcf} -t 2 -L -1 --report_readid
    """
}

// Process 7: Get final filter cuteSV insertions
process get_final_filter {
    tag "${sid}"
    publishDir "${params.outdir}/${sid}", mode: 'copy'
    
    input:
    tuple val(sid), path(cutesv_tumor_vcf), path(cutesv_normal_vcf), path(sronly_vcf), path(lronly_vcf), path(shared_vcf), path(lrvcf), path(normalbam), path(tumorbam)
    
    output:
    tuple val(sid), path("finalfilter_cutesv.vcf"), path("supread.loci.bed"), path("supread.tsv"), path(sronly_vcf), path(lronly_vcf), path(shared_vcf), path(lrvcf), path(normalbam), path(tumorbam), emit: final_filter_output
    
    script:
    """
    perl ${params.script_path}/04_get.finalfilter_cuteSV_ins.pl ./
    """
}

// Process 8: Process supporting reads
process process_supporting_reads {
    tag "${sid}"
    publishDir "${params.outdir}/${sid}", mode: 'copy'
    
    input:
    tuple val(sid), path(final_filter_vcf), path(supread_loci_bed), path(supread_tsv), path(sronly_vcf), path(lronly_vcf), path(shared_vcf), path(lrvcf), path(normalbam), path(tumorbam)
    
    output:
    tuple val(sid), path("finalfilter_cutesv.vcf"), path(sronly_vcf), path(lronly_vcf), path(shared_vcf), path(lrvcf), emit: processed_reads_output
    
    script:
    """
    ${params.samtools} depth -b ${supread_loci_bed} ${tumorbam} > supread.depth.tsv
    ${params.samtools} view -N ${supread_tsv} ${tumorbam} > supread.bam.tsv
    perl ${params.script_path}/04_2_correct.rb.pl ./
    rm finalfilter_cutesv2.vcf && mv finalfilter_cutesv3.vcf finalfilter_cutesv.vcf
    """
}

// Process 9: Polish all VCFs
process polish_all_vcfs {
    tag "${sid}"
    publishDir "${params.outdir}/${sid}", mode: 'copy'
    
    input:
    tuple val(sid), path(final_filter_vcf), path(sronly_vcf), path(lronly_vcf), path(shared_vcf), path(lrvcf)
    
    output:
    tuple val(sid), path("sronly.polished.bedpe"), path("sronly.polished.vcf"), 
          path("lronly.polished.bedpe"), path("lronly.polished.vcf"),
          path("shared.polished.bedpe"), path("shared.polished.vcf"),
          path(final_filter_vcf), path(lrvcf), emit: polished_all_output
    
    script:
    """
    perl ${params.script_path}/01_polish_survivor.pl ${sronly_vcf} sronly.polished.bedpe sronly.polished.vcf
    perl ${params.script_path}/01_polish_survivor.pl ${lronly_vcf} lronly.polished.bedpe lronly.polished.vcf
    perl ${params.script_path}/01_polish_survivor.pl ${shared_vcf} shared.polished.bedpe shared.polished.vcf
    """
}

// Process 10: Merge final VCFs
process merge_final_vcfs {
    tag "${sid}"
    publishDir "${params.outdir}/${sid}", mode: 'copy'
    
    input:
    tuple val(sid), path(sronly_polished_bedpe), path(sronly_polished_vcf), 
          path(lronly_polished_bedpe), path(lronly_polished_vcf),
          path(shared_polished_bedpe), path(shared_polished_vcf),
          path(final_filter_vcf), path(lrvcf)
    
    output:
    tuple val(sid), path("final.sr2lr.sv.vcf"), emit: merged_final_output
    
    script:
    """
    echo "${lrvcf}" > final.vcflist
    echo "${final_filter_vcf}" >> final.vcflist
    ${params.survivor} merge final.vcflist 1000 1 0 0 0 50 final.sr2lr.sv.vcf
    """
}

// Process 11: Final polish and compress
process final_polish_and_compress {
    tag "${sid}"
    publishDir "${params.outdir}/${sid}", mode: 'copy'
    
    input:
    tuple val(sid), path(final_sr2lr_sv_vcf)
    
    output:
    tuple val(sid), path("final.sr2lr.polished.vcf.gz"), path("final.sr2lr.polished.vcf.gz.tbi"), emit: final_output
    
    script:
    """
    perl ${params.script_path}/05_final_polish.pl ${final_sr2lr_sv_vcf} final.sr2lr.sv.bedpe final.sr2lr.polished.vcf
    ${params.bgzip} -c final.sr2lr.polished.vcf > final.sr2lr.polished.vcf.gz
    ${params.tabix} final.sr2lr.polished.vcf.gz
    """
}

// Main workflow
workflow {
    // Run the pipeline
    merge_lr_sr(samples_ch)
    polish_survivor(merge_lr_sr.out.merged_output)
    get_only_svs(polish_survivor.out.polished_output)
    make_vcf_for_force_calling(get_only_svs.out.only_svs_output)
    run_cutesv_normal(make_vcf_for_force_calling.out.vcf_for_fc_output)
    run_cutesv_tumor(run_cutesv_normal.out.cutesv_normal_output)
    get_final_filter(run_cutesv_tumor.out.cutesv_tumor_output)
    process_supporting_reads(get_final_filter.out.final_filter_output)
    polish_all_vcfs(process_supporting_reads.out.processed_reads_output)
    merge_final_vcfs(polish_all_vcfs.out.polished_all_output)
    final_polish_and_compress(merge_final_vcfs.out.merged_final_output)
}

// Display execution summary when the run completes
workflow.onComplete {
    log.info """
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    """
}
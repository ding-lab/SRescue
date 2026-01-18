# SRescue: Short- and Long-Read SV Integration Pipeline

A Nextflow pipeline for **integrating short-read (SR) and long-read (LR) structural variant (SV) calls**, rescuing SR-only SVs using LR data, and generating a polished, merged SV callset.

This pipeline is designed for **tumor–normal WGS** analyses and focuses on improving SV sensitivity by combining complementary sequencing technologies.

---

## Workflow Overview

1. **Merge LR and SR SVs**
   - Merge LR and SR VCFs using **SURVIVOR**

2. **Polish merged SVs**
   - Normalize and convert merged VCFs to BEDPE/VCF formats

3. **Classify SVs**
   - Identify:
     - SR-only SVs  
     - LR-only SVs  
     - Shared SVs

4. **Force calling SR-only SVs**
   - Generate VCF for force calling
   - Re-call SR-only SVs in normal and tumor BAMs using **cuteSV**

5. **Somatic filtering**
   - Compare tumor vs normal calls
   - Filter and retain high-confidence somatic SVs
   - Collect supporting read evidence

6. **Final integration**
   - Merge rescued SR SVs with LR SVs
   - Final polishing, sorting, compression, and indexing

---

## Inputs

- Long-read SV VCF  
- Short-read SV VCF  
- Tumor BAM (+ BAI)  
- Normal BAM (+ BAI)  
- Reference genome (FASTA)

---

## Outputs

- final.sr2lr.polished.vcf.gz
- final.sr2lr.polished.vcf.gz.tbi

---

## Dependencies

- **Nextflow**
- **SURVIVOR**
- **cuteSV**
- **samtools**
- **bgzip / tabix**
- **Perl** (custom polishing and filtering scripts)

---

## Usage

```bash
nextflow run main.nf \
--lrvcf_fn longread.vcf \
--srvcf_fn shortread.vcf \
--norm_fn normal.bam \
--tum_fn tumor.bam \
--reference reference.fa


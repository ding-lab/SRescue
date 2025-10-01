STORAGE_ALLOCATION="m.wyczalkowski"

NF="mnp_filter.nf"
CONFIG="mnp_filter.config"
OUTD="nextflow_out"
mkdir -p $OUTD


# /storage1/fs1/m.wyczalkowski/Active/Primary/CPTAC3.share/GDAN-GDC/GDC_import/data/c091d0df-6c94-4dd2-a363-6c9a8ff826c2/CTSP-B2ID-TTP1-A-1-0-D-A889-36.WholeGenome.RP-1329.bam

# /storage1/fs1/m.wyczalkowski/Active/Primary/CPTAC3.share/GDAN-GDC/GDC_import/data/c091d0df-6c94-4dd2-a363-6c9a8ff826c2/CTSP-B2ID-TTP1-A-1-0-D-A889-36.WholeGenome.RP-1329.bam


#LSF_DOCKER_VOLUMES="/scratch1/fs1/ris:/scratch1/fs1/ris /storage1/fs1/${STORAGE_ALLOCATION}/Active:/storage1/fs1/{$STORAGE_ALLOCATION}/Active /storage1/fs1/dinglab/Active:/storage1/fs1/dinglab/Active " \

LSF_DOCKER_VOLUMES="/scratch1/fs1/ris:/scratch1/fs1/ris /storage1/fs1/m.wyczalkowski/Active:/storage1/fs1/m.wyczalkowski/Active /storage1/fs1/dinglab/Active:/storage1/fs1/dinglab/Active" \
thpc-terminal bash -c "nextflow run $NF -c $CONFIG --outdir $OUTD $@"


#thpc-terminal bash -c "cd $OUTD; nextflow run $NF -profile test -c /scratch1/fs1/ris/application/nextflow/conf/rnaseq.config --outdir $OUTD"


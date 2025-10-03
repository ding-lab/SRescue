STORAGE_ALLOCATION="m.wyczalkowski"

NF="main.nf"
CONFIG="nextflow.config"
# OUTD="../nextflow_out" this is defined in the script at this time
mkdir -p $OUTD

#ARGS="--preview -with-dag Rescue.png -resume"
ARGS="-resume"

LSF_DOCKER_VOLUMES="/scratch1/fs1/ris:/scratch1/fs1/ris /storage1/fs1/m.wyczalkowski/Active:/storage1/fs1/m.wyczalkowski/Active /storage1/fs1/dinglab/Active:/storage1/fs1/dinglab/Active" \
thpc-terminal bash -c "nextflow run $NF -c $CONFIG $ARGS $@"
#thpc-terminal bash -c "nextflow run $NF -c $CONFIG --outdir $OUTD $@"




STORAGE_ALLOCATION="m.wyczalkowski"

NF="main.nf"
#CONFIG="nextflow-demo.config"
CONFIG="nextflow-HT970.config"
# OUTD="../nextflow_out" this is defined in the script at this time
#mkdir -p $OUTD

>&2 echo NOTE: not resuming deliberately
#ARGS="-resume"

# reports
RDIR="report"
REPORT_OUT="$RDIR/report.html"
TIMELINE_OUT="$RDIR/timeline.html"
FLOWCHART_OUT="$RDIR/flowchart.pdf"
# -with-trace writes out trace.txt
# -with-timeline [file name]
REPORT="-with-report $REPORT_OUT -with-trace -with-timeline $TIMELINE_OUT -with-dag $FLOWCHART_OUT"

LSF_DOCKER_VOLUMES="/scratch1/fs1/ris:/scratch1/fs1/ris /storage1/fs1/m.wyczalkowski/Active:/storage1/fs1/m.wyczalkowski/Active /storage1/fs1/dinglab/Active:/storage1/fs1/dinglab/Active" \
thpc-terminal bash -c "nextflow run $NF -c $CONFIG $ARGS $REPORT $@"
#thpc-terminal bash -c "nextflow run $NF -c $CONFIG --outdir $OUTD $@"




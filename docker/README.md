Image for running SRescue.  Includes the following dependencies

https://github.com/tjiangHIT/cuteSV
* requires the following:
    1. python3
    2. scipy
    2. pysam
    3. Biopython
    4. cigar
    5. numpy
    6. pyvcf3
    7. scikit-learn

-> it will be easier to make an image for just CuteSV 
here is an example: https://github.com/ncc-gap/module_box_aokad/blob/master/20220906-cuteSV/Dockerfile.txt

https://github.com/fritzsedlazeck/SURVIVOR
* C++ compiler needed

Also need
* bgzip
* samtools
* perl
* tabix

## perl 
based on TinDaisy/submodules/VEP_annotate/docker/Dockerfile



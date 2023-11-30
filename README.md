# enHunter
A pipeline to identify enhancer RNAs (eRNAs) from total RNA-seq data

## running enHunter on Terra 
see the /WDL for the .wdl scripts 

1. create Docker image using the scripts in /WDL/Docker.
docker build --tag name /WDL/Docker/script/Dockerfile

2. create Broad Methods Repository and import WDL script

## WDL scripts 
TNE_caller.wdl  -- detects transcribed noncoding elements 
bam2bigwig.19.wdl  -- converts alignment file (.bam) to normalized genomic coverage file (.bigwig)
combine_bigwig.24.wdl -- combine bigwig for multiple samples from the same group
star_splice_sites.wdl -- determine splice sites present in > 5 star junction files 

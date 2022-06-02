The bam files in the test folder were extracted from the SNCA gene region (chr4:90,616,701-90,787,997 in hg19) based on the BRAINcode data (Dong et al. Nature Neuroscience, 2018). Below is the code for genearting the files:
for i in HC*_SNDA_*; do echo $i; samtools view -b -o ~/neurogen/temp/scratch/$i.SNCA.bam $i/accepted_hits.bam chr4:90616701-90787997; done

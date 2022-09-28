# create a conda environment to install bedtools and samtools 

# to determine if bam2bigwig script is correctly splitting + and - strands  
# using the FOXP1 gene region as a test (- strand)

gsutil -u amp-pd-dawg-analyses cp gs://amp-pd-transcriptomics/samples/rnaseq/star/align-reads/BF-1001-SVM0_5T1/BF-1001-SVM0_5T1.star.bam ./
gsutil -u amp-pd-dawg-analyses cp gs://amp-pd-transcriptomics/samples/rnaseq/star/align-reads/BF-1001-SVM0_5T1/BF-1001-SVM0_5T1.star.bam.bai ./
samtools view -b -o FOXP1-BF-1001-SVM0_5T1.star.bam BF-1001-SVM0_5T1.star.bam chr3:70958139-71583756

# 94025343 is the number of PF_ALIGNED_READS for this sample 
total_mapped_reads=$((94025343 * 2))
RPMscale=$(bc <<< "scale=6;1000000/$total_mapped_reads") 

## previous script testing
bedtools genomecov -ibam FOXP1-BF-1001-SVM0_5T1.star.bam -bg -scale $RPMscale -split -strand + > FOXP1-BF-1001-SVM0_5T1.star.plus.bedgraph 
bedtools genomecov -ibam FOXP1-BF-1001-SVM0_5T1.star.bam -bg -scale $RPMscale -split -strand - > FOXP1-BF-1001-SVM0_5T1.star.minus.bedgraph 

gsutil -u amp-pd-dawg-analyses cp FOXP1-BF-1001-SVM0_5T1.star.plus.bedgraph gs://fc-secure-4300ea7f-8e77-4f20-a20b-38417286eaaf/
gsutil -u amp-pd-dawg-analyses cp FOXP1-BF-1001-SVM0_5T1.star.minus.bedgraph gs://fc-secure-4300ea7f-8e77-4f20-a20b-38417286eaaf/

## testing new script using deep tools 
#conda install -c bioconda deeptools

# changing total_mapped_reads num because only half the reads map to sense (other half map to antisense)
total_mapped_reads=94025343
RPMscale=$(bc <<< "scale=6;1000000/$total_mapped_reads") 

# deepTools requires an index 
samtools index FOXP1-BF-1001-SVM0_5T1.star.bam > FOXP1-BF-1001-SVM0_5T1.star.bam.bai

# note: split reads are natively handled by modules of deepTools 
# i.e. intron regions are given a signal of 0 
bamCoverage -b FOXP1-BF-1001-SVM0_5T1.star.bam --scaleFactor $RPMscale --filterRNAstrand forward --outFileFormat bedgraph -o dt-FOXP1-BF-1001-SVM0_5T1.star.plus.bedGraph
bamCoverage -b FOXP1-BF-1001-SVM0_5T1.star.bam --scaleFactor $RPMscale --filterRNAstrand reverse --outFileFormat bedgraph -o dt-FOXP1-BF-1001-SVM0_5T1.star.minus.bedGraph

# do I have to change binsize?


### testing plus strand gene: RXRA
#chr9:134,326,455-134,440,585 

samtools view -b -o RXRA-BF-1001-SVM0_5T1.star.bam BF-1001-SVM0_5T1.star.bam chr9:134326455-134440585 

total_mapped_reads=94025343
RPMscale=$(bc <<< "scale=6;1000000/$total_mapped_reads") 

# deepTools requires an index 
samtools index RXRA-BF-1001-SVM0_5T1.star.bam > RXRA-BF-1001-SVM0_5T1.star.bam.bai

# note: split reads are natively handled by modules of deepTools 
# i.e. intron regions are given a signal of 0 
bamCoverage -b RXRA-BF-1001-SVM0_5T1.star.bam --scaleFactor $RPMscale --filterRNAstrand forward --outFileFormat bedgraph -bs 10 -o dt-RXRA-BF-1001-SVM0_5T1.star.plus.bedGraph
bamCoverage -b RXRA-BF-1001-SVM0_5T1.star.bam --scaleFactor $RPMscale --filterRNAstrand reverse --outFileFormat bedgraph -bs 10 -o dt-RXRA-BF-1001-SVM0_5T1.star.minus.bedGraph

gsutil -u amp-pd-dawg-analyses cp dt-RXRA-BF-1001-SVM0_5T1.star.plus.bedGraph gs://fc-secure-4300ea7f-8e77-4f20-a20b-38417286eaaf/
gsutil -u amp-pd-dawg-analyses cp dt-RXRA-BF-1001-SVM0_5T1.star.minus.bedGraph gs://fc-secure-4300ea7f-8e77-4f20-a20b-38417286eaaf/


####### testing entire genome strand splitting 
total_mapped_reads=94025343
RPMscale=$(bc <<< "scale=6;1000000/$total_mapped_reads") 

# but first.. bin size testing
bamCoverage -b BF-1001-SVM0_5T1.star.bam --scaleFactor $RPMscale --filterRNAstrand forward --outFileFormat bigwig -bs 20 -o BF-1001-SVM0_5T1.star.bs20.plus.bw
bamCoverage -b BF-1001-SVM0_5T1.star.bam --scaleFactor $RPMscale --filterRNAstrand reverse --outFileFormat bigwig -bs 20 -o BF-1001-SVM0_5T1.star.bs20.minus.bw


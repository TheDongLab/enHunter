#!/bin/bash

## eRNA RPM files 
input=/data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs
#minus/eRNA.meanRPM.xls
#plus/eRNA.meanRPM.xls

## making the headers file 
# headers from both eRNA.meanRPM files should be the same 
head -n 1 $input/minus/eRNA.meanRPM.xls | sed 's/\t/\n/g' | sed '1d' > $input/all.RNAseq.samples.txt
#8258
awk 'OFS="\t" {split($1, a, "-"); print a[1]"-"a[2]}' all.RNAseq.samples.txt | sort | uniq > $input/all.RNAseq.participants.txt
#3207

# headers from Ruifeng's PPMI and PDBF files 
# participants in PPMI and PDBF files 
head -n 1 $input/eQTL/PDBF_SNP_Matrix.txt | sed 's/\t/\n/g' | sed '1,6d' | sort | uniq > all.PDBF.txt

# select for the participants from the wgs matrix for which you have RNAseq data for 
scripts=/data/bioinformatics/projects/donglab/AMPPD_eRNA/bin
# PDBF
bsub -q normal -n 1 -M 5000 "awk -f $scripts/select-matrix-cols.awk $input/all.RNAseq.participants.txt $input/eQTL/PDBF_SNP_Matrix.txt > eRNA.PDBF_SNP_Matrix.txt"

# PPMI
bsub -q normal -n 1 -M 5000 "awk -f $scripts/select-matrix-cols.awk $input/all.RNAseq.participants.txt $input/eQTL/PPMI_SNP_Matrix.txt > eRNA.PPMI_SNP_Matrix.txt"

#### select only one sample per participant ####
# file : participants with RNAseq samples and WGS genotyping = all.genotyping.with.RNAseq.txt
# file : participants with RNAseq samples = all.RNAseq.participants.txt # 3207 lines 

# look for baseline or BioFind samples
grep "BLM0T1\|SVM0_5T1" all.RNAseq.samples.txt > tmp.1 # 3110

# get the participants represented in the first step 
awk 'OFS="\t" {split($1, a, "-"); print a[1]"-"a[2]}' tmp.1 | sort | uniq > tmp.ppl

# get missing participants 
diff tmp.ppl all.RNAseq.participants.txt | sed 's/> //' | grep P > missing.ppl # 97

# check the RNAseq samples present for the missing ppl 
grep -f missing.ppl -F all.RNAseq.samples.txt > tmp.samples

# select samples with the next earliest time point i.e. SVM6T1 # 76
grep "SVM6T1" tmp.samples > tmp.2

# find the participants for these samples
awk 'OFS="\t" {split($1, a, "-"); print a[1]"-"a[2]}' tmp.2 | sort | uniq > tmp.2.ppl

# find the participants still missing after this filter 
diff tmp.2.ppl missing.ppl | sed 's/> //' | grep P > missing.ppl.2

# check the RNAseq samples present for the missing ppl pt 2 
grep -f missing.ppl.2 -F all.RNAseq.samples.txt | sort > tmp.samples.2

# manually pick out the rest of the samples 
cp tmp.samples.2 tmp.samples.3 # 21

# tmp.1 + tmp.2 + tmp.samples.3
cat tmp.1 tmp.2 tmp.samples.3 > RNAseq.samples.to.select.txt

# combine PDBF and PPMI datasets 
# check to see if the order of the snps are the same 
# bsub -q normal -n 1 -M 1000 "diff <(cut -f 1 eRNA.PDBF_SNP_Matrix.txt) <(cut -f 1 eRNA.PPMI_SNP_Matrix.txt) > PDBF.PPMI.diff.txt"
bsub -q normal -n 1 -M 1000 "cut -f 1 eRNA.PDBF_SNP_Matrix.txt > eRNA.PDBF_SNPs"
bsub -q normal -n 1 -M 1000 "cut -f 1 eRNA.PPMI_SNP_Matrix.txt > eRNA.PPMI_SNPs"

diff eRNA.PPMI_SNPs eRNA.PDBF_SNPs # output is empty! there is no difference 

bsub -q normal -n 1 -M 1000 "cut --complement -f 1 eRNA.PDBF_SNP_Matrix.txt > tmp" 
bsub -q normal -n 1 -M 1000 "paste eRNA.PPMI_SNP_Matrix.txt tmp > eRNA.SNP_Matrix.txt"

#### create one file with both mRNA and TNE expression values ####

# select for the RNAseq samples in eRNA and mRNA matrix 
bsub -q normal -n 1 -M 1000 "awk -f $scripts/select-eRNA-matrix-cols.awk $input/RNAseq.samples.to.select.txt $input/eRNA.merged.readCounts.v2.xls > eRNA.Matrix.txt"
bsub -q normal -n 1 -M 1000 "awk -f $scripts/select-eRNA-matrix-cols.awk $input/RNAseq.samples.to.select.txt $input/gene_exp/gene_expr_matrix_tpm_row_genes.txt > gene.Matrix.txt"

# combine mRNA and TNE matrixes 
cat gene.Matrix.txt <(sed '1d' eRNA.Matrix.txt) > eRNA.gene.Matrix.txt

# select for participants with genotyping data 
head -n 1 eRNA.gene.Matrix.txt | sed 's/\t/\n/g' | awk 'OFS="\t" {split($1, a, "-"); print a[1]"-"a[2]}' | tr '\n' '\t' > participant.header 
sed '1d' eRNA.gene.Matrix.txt > tmp 
cat participant.header tmp > eRNA.gene.participant.Matrix.txt

bsub -q normal -n 1 -M 1000 "awk -f $scripts/select-eRNA-matrix-cols.awk $input/all.genotyping.with.RNAseq.txt $input/eQTL/eRNA.gene.participant.Matrix.txt > eRNA.gene.genotyping.participant.Matrix.txt"

#### create covariates file for all tested participants ####
## moved all files to the amppd_data directory 
# reference: https://unix.stackexchange.com/questions/169995/rows-to-column-conversion-of-file 
cut -f 1,2,4,12 participant_info_full_results_04212021_unique_id.tsv | sed 's/Male/0/g'| sed 's/Female/1/g' | sed 's/Control/0/g' | sed 's/Case/1/g' | sed 's/Other/2/g' | awk 'BEGIN {OFS="\t"; FS="\t"}{ for (i=1; i<=NF; i++) RtoC[i]= (i in RtoC?RtoC[i] OFS :"") $i; } END{ for (i=1; i<=NF; i++) print RtoC[i] }' > covariate.matrix 

##### getting the positions #####
# snpid chr pos 
bsub -q normal -n 1 -M 1000 "cut -f 1,2,4 PDBF_SNP_Matrix.txt > all.snps"
awk 'OFS="\t" {print $2, "chr"$1, $3}' all.snps > snps.coordinates

# geneid chr left right 
# based on AMP PD gencode v29
join <(sort -k1,1 $input/gene_exp/genes_annotation_rh.tsv) <(cut -f 1 gene.Matrix.txt | sort) | awk 'OFS="\t" {print $1, $4, $5, $6}' > gene.coordinates
cut -f 1 eRNA.Matrix.txt |awk 'OFS="\t" {split($1,a,"_"); print $1, a[1], a[2], a[3]}' | sed '1d' > eRNA.coordinates 

cat gene.coordinates eRNA.coordinates > gene.eRNA.coordinates 

#grep -f all.genotyping.with.RNAseq.txt -F participant_info_full_results_04212021_unique_id.tsv | sort -k 1,1 > genotyping.with.RNAseq.participants.txt 
#2984
#diff all.genotyping.with.RNAseq.txt <(cut -f 1 genotyping.with.RNAseq.participants.txt )

bsub -q normal -n 1 -M 10000 Rscript matrixeQTL.R

## files:
## snps = eRNA.SNP_Matrix.txt
## snp coordinates = snps.coordinates
## gene expression = eRNA.gene.genotyping.participant.Matrix.txt 
## gene coordinates = gene.eRNA.coordinates
## covariates = covariate.matrix 

###### REORDER gene expression and covariate samples to be the same order as snps files ######
awk -f $scripts/select-eRNA-matrix-cols.awk participant.order.txt covariate.matrix > covariate.reordered.txt 
bsub -q normal -n 1 -M 1000 "awk -f $scripts/select-eRNA-matrix-cols.awk participant.order.txt genes/eRNA.gene.Matrix.txt > genes/eRNA.gene.Matrix.reordered.txt"

###### remove extra tab at the end of each line ######
#https://superuser.com/questions/217463/removing-the-last-tab-of-a-line-in-a-file 
bsub -q normal -n 1 -M 1000 "sed '1d' eRNA.SNP_Matrix.txt | sed 's/\t$//' | cat header - > eRNA.SNP_Matrix.revised.txt"







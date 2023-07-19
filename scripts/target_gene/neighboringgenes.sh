#!/bin/bash
#### usage #####
# ./neighboring_genes.sh LITAF_loci.bed 1000000 /data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs/gene_exp/gene_expr_matrix_tpm_row_genes.txt LITAF_locus_expr

module load ucsc/default
module load bedtools/2.23.0
module load R/3.1.0

##### inputs ######
## 1. genomic coordinates 
## 2. number of bp upstream and downstream 
###################

##### outputs ######
## list of genes located x amount of bp up and downstream of the provided coordinate 
####################
inputbed=$1
dist=$2
matrix=$3
output=$4

index=hg38
GENOME=/data/bioinformatics/referenceGenome/Homo_sapiens/UCSC/$index

## step 1. find the genes which are within the search area (up and downstream of the loci)
awk -v dist=$dist 'OFS="\t" {print $1, $2 - dist, $3 + dist, $4 }' $inputbed | intersectBed -a - -b <(sortBed -i $GENOME/Annotation/Genes/gencode.v37.annotation.genes.bed | grep chr ) -wb | cut -f 8 | awk -F'___' '{split($2,a,"."); print a[1]}'> neighboring_genes.tmp
### 1,000,000 bp up and down stream 
### chr16 11608526  11643272 eRNA

## step 2. for each gene found in step 1, find the sample + gene expression values 
cat <(head -n 1 $matrix) <(grep -f neighboring_genes.tmp $matrix) > $output.txt


# while read gene; do
#   row=<(echo "$gene"| cut -f1 -d".")
#   cat <(cut -f 1 $matrix) <(grep $row $matrix) 
# done < neighboring_genes.txt
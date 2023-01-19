#!/bin/bash

# getting the eRNA and eRNA_pval
# on erisone
cd /data/bioinformatics/projects/donglab/AMPPD_eRNA/DE_eRNA/PPMI/All

cut -f1,7 DEresult.padj05.xls

sed '1d' FINAL.xls | awk 'BEGIN{OFS="\t"; FS="\t"} {$2 == "+"?strand="plus":strand="minus"; print $1"_"strand, $8, $10}'



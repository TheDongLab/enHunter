#!/bin/bash
# bedtools clustering distance tests 

module load bedtools/2.26.0

# input data for cluster must be presorted by chromosome and start position 
inputfile=/data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs/eRNA_stranded_sorted.bed

bedtools cluster -i $inputfile -d 300 > eRNA_cluster_300bp.bed

bedtools cluster -i $inputfile -d 1000 > eRNA_cluster_1kbp.bed

bedtools cluster -i $inputfile -d 3000 > eRNA_cluster_3kbp.bed

bedtools cluster -i $inputfile -d 5000 > eRNA_cluster_5kbp.bed

bedtools cluster -i $inputfile -d 6000 > eRNA_cluster_6kbp.bed

bedtools cluster -i $inputfile -d 7000 > eRNA_cluster_7kbp.bed

bedtools cluster -i $inputfile -d 8000 > eRNA_cluster_8kbp.bed

bedtools cluster -i $inputfile -d 200000 > eRNA_cluster_200kbp.bed

module load bedtools/2.20.1

bedtools groupby -i eRNA_cluster_5kbp.bed -g 7 -c 7 -o count > eRNA_cluster_5kb.grouped
bedtools groupby -i eRNA_cluster_8kbp.bed -g 7 -c 7 -o count > eRNA_cluster_8kb.grouped

join <(bedtools groupby -i eRNA_cluster_12.5kbp.bed -g 7,1 -c 3 -o max) <(bedtools groupby -i eRNA_cluster_12.5kbp.bed -g 7 -c 2 -o min) | awk 'OFS="\t" {print $1, $2, $4, $3, $3 - $4}' > eRNA_cluster_12.5kbp.distance

join <(bedtools groupby -i eRNA_cluster_50kbp.bed -g 7,1 -c 3 -o max) <(bedtools groupby -i eRNA_cluster_50kbp.bed -g 7 -c 2 -o min) | awk 'OFS="\t" {print $1, $2, $4, $3, $3 - $4}' > eRNA_cluster_50kbp.distance

join <(bedtools groupby -i eRNA_cluster_200kbp.bed -g 7,1 -c 3 -o max) <(bedtools groupby -i eRNA_cluster_200kbp.bed -g 7 -c 2 -o min) | awk 'OFS="\t" {print $1, $2, $4, $3, $3 - $4}' > eRNA_cluster_200kbp.distance

join <(bedtools groupby -i eRNA_cluster_100kbp.bed -g 7,1 -c 3 -o max) <(bedtools groupby -i eRNA_cluster_100kbp.bed -g 7 -c 2 -o min) | awk 'OFS="\t" {print $1, $2, $4, $3, $3 - $4}' > eRNA_cluster_100kbp.distance

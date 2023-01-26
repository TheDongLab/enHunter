#!/bin/bash
### this code is *copied* from bidir_pair_matrix.sh
### script for getting merged raw read count number for TNE pairs 
### found in bidirectional_pairs_all (no distance filter)

# bsub -q normal -J bidirmatrix -o %J.out -e %J.err -R 'rusage [mem=8000]' ./bidirectional_pairs_matrix.sh

cd /data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs

head -n 1 eRNA.merged.readCounts.v2.xls > eRNA.bidirectional_pairs.readCounts.xls

cat bidirectional_pairs/bidirectional_pairs_all.txt| sed 's/+/plus/g; s/-/minus/g' | 
while read TNE1 TNE1_strand TNE2 TNE2_strand dist 
do 
  TNE1+="_"$TNE1_strand # plus strand
  TNE2+="_"$TNE2_strand # minus strand
  combined=$TNE1"_"$TNE2
  
  bsub -q vshort -n 1 -J $combined bash merge_bidir_TNEs.sh $TNE1 $TNE2 $combined
done 


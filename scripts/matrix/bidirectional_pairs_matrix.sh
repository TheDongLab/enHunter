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

### running less number of pairs 
head -n 1 eRNA.bidirectional_pairs_all.readCounts.xls > eRNA.bidirectional_pairs.readCounts.xls

cat bidirectional_pairs_700.txt| sed 's/+/plus/g; s/-/minus/g' | 
while IFS=$'\t' read TNE1 TNE1_strand TNE2 TNE2_strand dist 
do 
  TNE1+="_"$TNE1_strand # plus strand
  TNE2+="_"$TNE2_strand # minus strand
  combined=$TNE1"_"$TNE2
  
  #bsub -q vshort -n 1 -J $combined bash merge_bidir_TNEs.sh $TNE1 $TNE2 $combined
  TNE2_vals=$(grep -F $TNE2 /data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs/eRNA.merged.readCounts.v2.xls)
  IFS=$'\t' read -r -a array_minus <<< "$TNE2_vals"
  TNE1_vals=$(grep -F $TNE1 /data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs/eRNA.merged.readCounts.v2.xls)
  IFS=$'\t' read -r -a array_plus <<< "$TNE1_vals"
  # gets the length of array 
  # TNE1_vals and TNE2_vals should be of equal length! 
  # aka array_minus and array_plus
  len=${#array_plus[@]}
  echo "${#array_plus[@]}"
  for ((i=1; i <$len; i++));do c+=(`expr ${array_plus[$i]} + ${array_minus[$i]}`);done
  (
    IFS=$'\t'
    echo -e "$combined\t${c[*]}"
  ) >> eRNA.bidirectional_pairs.readCounts.xls

done 


### testing why this doesn't work perfectly ###
cat bidirectional_pairs_700.txt| sed 's/+/plus/g; s/-/minus/g' | 
while IFS=$'\t' read TNE1 TNE1_strand TNE2 TNE2_strand dist 
do 
  TNE1+="_"$TNE1_strand # plus strand
  TNE2+="_"$TNE2_strand # minus strand
  combined=$TNE1"_"$TNE2
  
  (
  echo -e "$combined"
  )  >> test.txt
done 

### the above code works as expected ... it probably has to do with submitting to the cluster multiple times... 
### going to try without submitting to cluster 

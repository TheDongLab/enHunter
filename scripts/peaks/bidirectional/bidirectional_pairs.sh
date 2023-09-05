#!/bin/bash
# usage: ./bidirectional_pairs.sh input_files/closest/plus.minus.bed input_files/closest/minus.plus.bed 900

##### 1. both TNEs are in minus to plus orientation (negative distance)
##### 2. summit distance is < x bp 

plusminus=$1
minusplus=$2
dist=$3

mkdir scripts/peaks/bidirectional/pairs/$dist
cd scripts/peaks/bidirectional/pairs/$dist

##### 1. both TNEs are in minus to plus orientation (negative distance)

# plus.minus.bed 
#input_files/closest/plus.minus.bed
awk '$21<0 {print $0}' $plusminus > plus.minus.step1.bed
wc -l plus.minus.step1.bed
# there are 37613 pairs in minus to plus orientation for plus.minus.bed 

# input_files/closest/minus.plus.bed
awk '$21<0 {print $0}' $minusplus > minus.plus.step1.bed
wc -l minus.plus.step1.bed
# there are 24078 pairs in minus to plus orientation for minus.plus.bed 

##### 2. summit distance is < dist bp 

# plus.minus.bed 
awk -v dist="${dist}" '$21>-dist {print $0}' plus.minus.step1.bed > plus.minus.step2.bed
wc -l plus.minus.step2.bed

# minus.plus.bed 
awk -v dist="${dist}" '$21>-dist {print $0}' minus.plus.step1.bed > minus.plus.step2.bed
wc -l minus.plus.step2.bed

## merge minus.plus.step2.bed and plus.minus.step2.bed, checking for duplicate pairs 
## reorganize minus.plus.step2.bed into plus minus TNE listing
#system(paste0("cd scripts/peaks/bidirectional/pairs/",dist,"bp") )
paste <(cut -f11-20 minus.plus.step2.bed) <(cut -f1-10 minus.plus.step2.bed) <(cut -f21 minus.plus.step2.bed) > plus.minus.step2.tmp
cat plus.minus.step2.tmp plus.minus.step2.bed | sort | uniq > merged.bed
wc -l merged.bed

cut -f4,6,14,16,21 merged.bed > bidirectional_pairs.txt

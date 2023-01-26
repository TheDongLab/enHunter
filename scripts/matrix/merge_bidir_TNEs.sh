#!/bin/bash

# inputs 
TNE1=$1 # plus strand
TNE2=$2 # minus strand
combined=$3
  
TNE2_vals=$(grep -F $TNE2 eRNA.merged.readCounts.v2.xls)
IFS=$'\t' read -r -a array_minus <<< "$TNE2_vals"
  
TNE1_vals=$(grep -F $TNE1 eRNA.merged.readCounts.v2.xls)
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
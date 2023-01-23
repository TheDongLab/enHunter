## ====================================================== ##
##### raw read count matrix for bi directional pairs 
##### minus strand read counts + plus strand read counts
## ====================================================== ##

# TNE pairs 
awk 'OFS="\t" {
$2 == "+"?strand_1="plus":strand_1="minus";
$4 == "+"?strand_2="plus":strand_2="minus";
print $1"_"strand_1"_"$3"_"strand_2
}' scripts/peaks/bidirectional/pairs/bidirectional_pairs.txt > bidir_pairs.txt
## there are 372 pairs total 

# 
sed 's/+/plus/g; s/-/minus/g' scripts/peaks/bidirectional/pairs/bidirectional_pairs_all.txt | 
while read TNE1 TNE1_strand TNE2 TNE2_strand dist 
do 
  TNE1+="_"$TNE1_strand # plus strand
  TNE2+="_"$TNE2_strand # minus strand
  combined=$TNE1"_"$TNE2
  
  grep $TNE1 eRNA.minus.readCounts.xls 
done 

## for testing purposes 
head  -n 1 bidirectional_pairs/bidirectional_pairs_all.txt| sed 's/+/plus/g; s/-/minus/g' | 
while read TNE1 TNE1_strand TNE2 TNE2_strand dist 
do 
  TNE1+="_"$TNE1_strand # plus strand
  TNE2+="_"$TNE2_strand # minus strand
  combined=$TNE1"_"$TNE2
  
  TNE2_vals=$(grep -F $TNE2 eRNA.merged.readCounts.v2.xls)
  IFS=$'\t' read -r -a array_minus <<< "$TNE2_vals"
  
  TNE1_vals=$(grep -F $TNE1 eRNA.merged.readCounts.v2.xls)
  IFS=$'\t' read -r -a array_plus <<< "$TNE1_vals"
  
  # gets the length of array 
  # TNE1_val and TNE2_val should be of equal length! 
  len=${#array_plus[@]}
  echo "${#array_plus[@]}"
  
  for ((i=1; i <$len; i++));do c+=(`expr ${array_plus[$i]} + ${array_minus[$i]}`);done
  
  (
    IFS=$'\t'
    echo "${c[*]}"
  ) > test.txt
done 



### double checking to make sure the column names are the same for plus and minus eRNA files 
### run on erisone 
### /data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs

plus=$(awk 'NR==1' minus/eRNA.minus.readCounts.xls) 
minus=$(awk 'NR==1' plus/eRNA.plus.readCounts.xls) 

if [[ $plus == $minus ]]; then 
echo "the same"
fi

#!/bin/bash

##### dopamine neurons RPM matrix 

curl https://www.humanbraincode.org/download/data/BRAINcode.TNE_SNDA.meanRPM.allSamples.xls.gz -o BRAINcode.TNE_SNDA.meanRPM.allSamples.xls

cut -f1 BRAINcode.TNE_SNDA.meanRPM.allSamples.xls | sed '1d' | awk 'BEGIN{OFS="\t"; FS="_"} {print $1, $2, $3, $1"_"$2"_"$3}' > BRAINcode.dopamine.TNE.hg19.bed

# convert from hg19 to hg38
module load ucsc/default
liftOver BRAINcode.dopamine.TNE.hg19.bed /data/bioinformatics/external_data/externalData/hg19ToHg38.over.chain.gz BRAINcode.dopamine.TNE.hg38.bed unMapped 

# overlapping dopamine and blood TNEs
module load bedtools/2.26.0 
bedtools intersect -a BRAINcode.dopamine.TNE.hg38.bed -b /data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs/eRNA_stranded_sorted.bed -wa -wb > TNE.dopamine.and.blood.txt

# getting the unique number of dopamine TNEs which were also identified in blood 
cut -f 4 TNE.dopamine.and.blood.txt | sort | uniq | wc -l # 6533

### finding out class information 
join <(cut -f 1,28 braincode.eRNA.characterize.xls | sort -k1,1) <(cut -f 4,8 TNE.dopamine.and.blood.txt) > tmp.txt

join <(cut -f 1,24 /data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs/eRNA.characterize.feature.color.xls | sort -k1,1) <(sort -k3,3 tmp.txt) -1 1 -2 3 > TNE.dopamine.and.blood.class.txt
# note: after this join, blood TNE name and class are moved to the front 

# adding a header
#sed -i '1s/^/dopamine\tblood\n/'

## getting class count information 
module load bedtools/2.20.1 
bedtools groupby -i <(awk 'OFS="\t" {print $1,$2,$3,$4}' TNE.dopamine.and.blood.class.txt | sort -k2,2 -k4,4) -g 2,4 -c 1 -o count > TNE.dopamine.and.blood.class.counts 








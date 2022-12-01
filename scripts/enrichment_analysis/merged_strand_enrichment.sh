#!/bin/bash 
### script for running enrichment analysis by combining plus and minus strand TNEs and seperate by class 
cd /data/bioinformatics/projects/donglab/AMPPD_eRNA/output/

# select for unique SNP and tissue pairs 
cat ./minus/eRNA.minus.f16.GWASDisease.intersect ./plus/eRNA.plus.f16.GWASDisease.intersect > ./merged/eRNA.f16.GWASDisease.intersect
cat ./minus/eRNA.minus.f18.eSNP.Disease.intersect ./plus/eRNA.plus.f18.eSNP.Disease.intersect > ./merged/eRNA.f18.eSNP.Disease.intersect
cat ./minus/eRNA.minus.f18.eSNP.gtexCaviarDisease.intersect ./plus/eRNA.plus.f18.eSNP.gtexCaviarDisease.intersect > ./merged/eRNA.f18.eSNP.gtexCaviarDisease.intersect
cat ./minus/eRNA.minus.f18.eSNP.gtexDapgDisease.intersect ./plus/eRNA.plus.f18.eSNP.gtexCaviarDisease.intersect > ./merged/eRNA.f18.eSNP.gtexCaviarDisease.intersect

## GWAS disease enrichment 
# If one SNP is overlapped by two TNEs (this is the max num bc TNEs in the same strand can not overlap so one must be plus and the other is minus), that SNP should still be counted only once 
grep -F -f ../inputs/class/class1_TNEs.txt ./merged/eRNA.f16.GWASDisease.intersect > ./merged/class1.eRNA.f16.GWASDisease.intersect
grep -F -f ../inputs/class/class2_TNEs.txt ./merged/eRNA.f16.GWASDisease.intersect > ./merged/class2.eRNA.f16.GWASDisease.intersect
grep -F -f ../inputs/class/class3_TNEs.txt ./merged/eRNA.f16.GWASDisease.intersect > ./merged/class3.eRNA.f16.GWASDisease.intersect

mkdir ./merged/eRNA.enrichment
module load bedtools/2.23.0
cut -f8,9  ./merged/class1.eRNA.f16.GWASDisease.intersect |sort -k 2,2 | bedtools groupby -g 9 -c 9 -o count | sort -k 1 | bedtools groupby -g 1 -c 2 -o sum > ./merged/eRNA.enrichment/class1.eRNA.f16.GWASDisease.counts
sort -k 9  ./merged/class2.eRNA.f16.GWASDisease.intersect | bedtools groupby -g 9 -c 9 -o count | sort -k 1 | bedtools groupby -g 1 -c 2 -o sum > ./merged/eRNA.enrichment/class2.eRNA.f16.GWASDisease.counts
sort -k 9  ./merged/class3.eRNA.f16.GWASDisease.intersect | bedtools groupby -g 9 -c 9 -o count | sort -k 1 | bedtools groupby -g 1 -c 2 -o sum > ./merged/eRNA.enrichment/class3.eRNA.f16.GWASDisease.counts

scp rw552@erisone.partners.org:/data/bioinformatics/projects/donglab/AMPPD_eRNA/output/merged/eRNA.enrichment/class'*'.eRNA.f16.GWASDisease.counts /Users/rw552/Documents/enHunter/input_files/characterization/feature.enrichment/counts/merged

## eSNP disease enrichment 
grep -F -f ../inputs/class/class1_TNEs.txt ./merged/eRNA.f18.eSNP.Disease.intersect > ./merged/class1.eRNA.f18.eSNP.Disease.intersect
grep -F -f ../inputs/class/class2_TNEs.txt ./merged/eRNA.f18.eSNP.Disease.intersect > ./merged/class2.eRNA.f18.eSNP.Disease.intersect
grep -F -f ../inputs/class/class3_TNEs.txt ./merged/eRNA.f18.eSNP.Disease.intersect > ./merged/class3.eRNA.f18.eSNP.Disease.intersect

mkdir ./merged/eRNA.enrichment
module load bedtools/2.23.0
sort -k 9  ./merged/class1.eRNA.f18.eSNP.Disease.intersect | bedtools groupby -g 9 -c 9 -o count > ./merged/eRNA.enrichment/class1.eRNA.f18.eSNP.Disease.counts
sort -k 9  ./merged/class2.eRNA.f18.eSNP.Disease.intersect| bedtools groupby -g 9 -c 9 -o count > ./merged/eRNA.enrichment/class2.eRNA.f18.eSNP.Disease.counts
sort -k 9  ./merged/class3.eRNA.f18.eSNP.Disease.intersect | bedtools groupby -g 9 -c 9 -o count > ./merged/eRNA.enrichment/class3.eRNA.f18.eSNP.Disease.counts

scp rw552@erisone.partners.org:/data/bioinformatics/projects/donglab/AMPPD_eRNA/output/merged/eRNA.enrichment/class'*'.eRNA.f16.GWASDisease.counts /Users/rw552/Documents/enHunter/input_files/characterization/feature.enrichment/counts/merged

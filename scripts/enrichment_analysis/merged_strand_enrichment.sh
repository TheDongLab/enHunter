#!/bin/bash 
### script for running enrichment analysis by combining plus and minus strand TNEs and seperate by class 
cd /data/bioinformatics/projects/donglab/AMPPD_eRNA/output/

# select for unique SNP and tissue pairs 
cat ./minus/eRNA.minus.f16.GWASDisease.intersect ./plus/eRNA.plus.f16.GWASDisease.intersect > ./merged/eRNA.f16.GWASDisease.intersect
cat ./minus/eRNA.minus.f18.eSNP.Disease.intersect ./plus/eRNA.plus.f18.eSNP.Disease.intersect > ./merged/eRNA.f18.eSNP.Disease.intersect
cat ./minus/eRNA.minus.f18.eSNP.gtexCaviarDisease.intersect ./plus/eRNA.plus.f18.eSNP.gtexCaviarDisease.intersect > ./merged/eRNA.f18.eSNP.gtexCaviarDisease.intersect
cat ./minus/eRNA.minus.f18.eSNP.gtexDapgDisease.intersect ./plus/eRNA.plus.f18.eSNP.gtexDapgDisease.intersect  > ./merged/eRNA.f18.eSNP.gtexDapgDisease.intersect

## GWAS disease enrichment 
grep -F -f ../inputs/class/class1_TNEs.txt ./merged/eRNA.f16.GWASDisease.intersect > ./merged/class1.eRNA.f16.GWASDisease.intersect
grep -F -f ../inputs/class/class2_TNEs.txt ./merged/eRNA.f16.GWASDisease.intersect > ./merged/class2.eRNA.f16.GWASDisease.intersect
grep -F -f ../inputs/class/class3_TNEs.txt ./merged/eRNA.f16.GWASDisease.intersect > ./merged/class3.eRNA.f16.GWASDisease.intersect

mkdir ./merged/eRNA.enrichment
module load bedtools/2.23.0
# If one SNP is overlapped by two TNEs (this is the max num bc TNEs in the same strand can not overlap so one must be plus and the other is minus), that SNP should still be counted only once 
# search for unique SNP-disease pairs 

cut -f8,9  ./merged/class1.eRNA.f16.GWASDisease.intersect | sort -k 2,2 | uniq | bedtools groupby -g 2 -c 2 -o count | sort -k 1 | bedtools groupby -g 1 -c 2 -o sum > ./merged/eRNA.enrichment/class1.eRNA.f16.GWASDisease.counts
cut -f8,9  ./merged/class2.eRNA.f16.GWASDisease.intersect | sort -k 2,2 | uniq | bedtools groupby -g 2 -c 2 -o count | sort -k 1 | bedtools groupby -g 1 -c 2 -o sum > ./merged/eRNA.enrichment/class2.eRNA.f16.GWASDisease.counts
cut -f8,9  ./merged/class3.eRNA.f16.GWASDisease.intersect | sort -k 2,2 | uniq | bedtools groupby -g 2 -c 2 -o count | sort -k 1 | bedtools groupby -g 1 -c 2 -o sum > ./merged/eRNA.enrichment/class3.eRNA.f16.GWASDisease.counts

scp rw552@erisone.partners.org:/data/bioinformatics/projects/donglab/AMPPD_eRNA/output/merged/eRNA.enrichment/class'*'.eRNA.f16.GWASDisease.counts .

#wd=/Users/rosanwang/Documents/college/dong_lab/code/enHunter/input_files/characterization/feature.enrichment/counts
wd=/Users/rw552/Documents/amp-pd/enHunter/input_files/characterization/feature.enrichment/counts
Rscript fisher_test.R $wd/merged/class1.eRNA.f16.GWASDisease.counts $wd/GWAS_20220810.v1.02.counts.v2 GWAS_class1.pdf "GWAS class 1"
Rscript fisher_test.R $wd/merged/class2.eRNA.f16.GWASDisease.counts $wd/GWAS_20220810.v1.02.counts.v2 GWAS_class2.pdf "GWAS class 2"
Rscript fisher_test.R $wd/merged/class3.eRNA.f16.GWASDisease.counts $wd/GWAS_20220810.v1.02.counts.v2 GWAS_class3.pdf "GWAS class 3"

#/Users/rw552/Documents/enHunter/input_files/characterization/feature.enrichment/counts/merged

## eSNP disease enrichment 
# caviar
grep -F -f ../inputs/class/class1_TNEs.txt ./merged/eRNA.f18.eSNP.gtexCaviarDisease.intersect > ./merged/class1.eRNA.f18.eSNP.gtexCaviarDisease.intersect
grep -F -f ../inputs/class/class2_TNEs.txt ./merged/eRNA.f18.eSNP.gtexCaviarDisease.intersect > ./merged/class2.eRNA.f18.eSNP.gtexCaviarDisease.intersect
grep -F -f ../inputs/class/class3_TNEs.txt ./merged/eRNA.f18.eSNP.gtexCaviarDisease.intersect > ./merged/class3.eRNA.f18.eSNP.gtexCaviarDisease.intersect

module load bedtools/2.23.0
cut -f8,9 ./merged/class1.eRNA.f18.eSNP.gtexCaviarDisease.intersect | sort -k 2,2 | uniq | bedtools groupby -g 2 -c 2 -o count > ./merged/eRNA.enrichment/class1.eRNA.f18.eSNP.gtexCaviarDisease.counts
cut -f8,9 ./merged/class2.eRNA.f18.eSNP.gtexCaviarDisease.intersect | sort -k 2,2 | uniq | bedtools groupby -g 2 -c 2 -o count > ./merged/eRNA.enrichment/class2.eRNA.f18.eSNP.gtexCaviarDisease.counts
cut -f8,9 ./merged/class3.eRNA.f18.eSNP.gtexCaviarDisease.intersect | sort -k 2,2 | uniq | bedtools groupby -g 2 -c 2 -o count > ./merged/eRNA.enrichment/class3.eRNA.f18.eSNP.gtexCaviarDisease.counts

scp rw552@erisone.partners.org:/data/bioinformatics/projects/donglab/AMPPD_eRNA/output/merged/eRNA.enrichment/class'*'.eRNA.f18.eSNP.gtexCaviarDisease.counts .

Rscript fisher_test.R $wd/merged/class1.eRNA.f18.eSNP.gtexCaviarDisease.counts $wd/snp_gtexCaviar.counts gtexCaviar_class1.pdf "gtexCaviar class 1"
Rscript fisher_test.R $wd/merged/class2.eRNA.f18.eSNP.gtexCaviarDisease.counts $wd/snp_gtexCaviar.counts gtexCaviar_class2.pdf "gtexCaviar class 2"
Rscript fisher_test.R $wd/merged/class3.eRNA.f18.eSNP.gtexCaviarDisease.counts $wd/snp_gtexCaviar.counts gtexCaviar_class3.pdf "gtexCaviar class 3"

# dapg
grep -F -f ../inputs/class/class1_TNEs.txt ./merged/eRNA.f18.eSNP.gtexDapgDisease.intersect > ./merged/class1.eRNA.f18.eSNP.gtexDapgDisease.intersect
grep -F -f ../inputs/class/class2_TNEs.txt ./merged/eRNA.f18.eSNP.gtexDapgDisease.intersect > ./merged/class2.eRNA.f18.eSNP.gtexDapgDisease.intersect
grep -F -f ../inputs/class/class3_TNEs.txt ./merged/eRNA.f18.eSNP.gtexDapgDisease.intersect > ./merged/class3.eRNA.f18.eSNP.gtexDapgDisease.intersect

module load bedtools/2.23.0
cut -f8,9 ./merged/class1.eRNA.f18.eSNP.gtexDapgDisease.intersect | sort -k 2,2 | uniq | bedtools groupby -g 2 -c 2 -o count > ./merged/eRNA.enrichment/class1.eRNA.f18.eSNP.gtexDapgDisease.counts
cut -f8,9 ./merged/class2.eRNA.f18.eSNP.gtexDapgDisease.intersect | sort -k 2,2 | uniq | bedtools groupby -g 2 -c 2 -o count > ./merged/eRNA.enrichment/class2.eRNA.f18.eSNP.gtexDapgDisease.counts
cut -f8,9 ./merged/class3.eRNA.f18.eSNP.gtexDapgDisease.intersect | sort -k 2,2 | uniq | bedtools groupby -g 2 -c 2 -o count > ./merged/eRNA.enrichment/class3.eRNA.f18.eSNP.gtexDapgDisease.counts

scp rw552@erisone.partners.org:/data/bioinformatics/projects/donglab/AMPPD_eRNA/output/merged/eRNA.enrichment/class'*'.eRNA.f18.eSNP.gtexDapgDisease.counts .

Rscript fisher_test.R $wd/merged/class1.eRNA.f18.eSNP.gtexDapgDisease.counts $wd/snp_gtexDapg.counts gtexDapg_class1.pdf "gtexDapg class 1"
Rscript fisher_test.R $wd/merged/class2.eRNA.f18.eSNP.gtexDapgDisease.counts $wd/snp_gtexDapg.counts gtexDapg_class2.pdf "gtexDapg class 2"
Rscript fisher_test.R $wd/merged/class3.eRNA.f18.eSNP.gtexDapgDisease.counts $wd/snp_gtexDapg.counts gtexDapg_class3.pdf "gtexDapg class 3"

#### bidirectional pairs (144 pairs)
grep -F -f ../inputs/class/all_bidirectional_class_TNEs.txt ./merged/eRNA.f18.eSNP.gtexDapgDisease.intersect > ./merged/bidir.eRNA.f18.eSNP.gtexDapgDisease.intersect
grep -F -f ../inputs/class/all_bidirectional_class_TNEs.txt ./merged/eRNA.f18.eSNP.gtexCaviarDisease.intersect > ./merged/bidir.eRNA.f18.eSNP.gtexCaviarDisease.intersect
grep -F -f ../inputs/class/all_bidirectional_class_TNEs.txt ./merged/eRNA.f16.GWASDisease.intersect > ./merged/bidir.eRNA.f16.GWASDisease.intersect

module load bedtools/2.23.0
cut -f8,9 ./merged/bidir.eRNA.f18.eSNP.gtexDapgDisease.intersect | sort -k 2,2 | uniq | bedtools groupby -g 2 -c 2 -o count > ./merged/eRNA.enrichment/bidir.eRNA.f18.eSNP.gtexDapgDisease.counts
cut -f8,9 ./merged/bidir.eRNA.f18.eSNP.gtexCaviarDisease.intersect | sort -k 2,2 | uniq | bedtools groupby -g 2 -c 2 -o count > ./merged/eRNA.enrichment/bidir.eRNA.f18.eSNP.gtexCaviarDisease.counts
cut -f8,9  ./merged/bidir.eRNA.f16.GWASDisease.intersect | sort -k 2,2 | uniq | bedtools groupby -g 2 -c 2 -o count | sort -k 1 | bedtools groupby -g 1 -c 2 -o sum > ./merged/eRNA.enrichment/bidir.eRNA.f16.GWASDisease.counts

scp rw552@erisone.partners.org:/data/bioinformatics/projects/donglab/AMPPD_eRNA/output/merged/eRNA.enrichment/bidir.eRNA.'*'.counts .

wd=/Users/rosanwang/Documents/college/dong_lab/code/enHunter/input_files/characterization/feature.enrichment/counts
Rscript fisher_test.R $wd/merged/bidir.eRNA.f18.eSNP.gtexDapgDisease.counts $wd/snp_gtexDapg.counts gtexDapg_bidir.pdf "gtexDapg bidir"
Rscript fisher_test.R $wd/merged/bidir.eRNA.f18.eSNP.gtexCaviarDisease.counts $wd/snp_gtexCaviar.counts gtexCaviar_bidir.pdf "gtexCaviar bidir"
Rscript fisher_test.R $wd/merged/bidir.eRNA.f16.GWASDisease.counts $wd/GWAS_20220810.v1.02.counts.v2 GWAS_bidir.pdf "GWAS bidir"


##### on ALL TNEs 
module load bedtools/2.23.0
cut -f8,9 ./merged/eRNA.f18.eSNP.gtexDapgDisease.intersect | sort -k 2,2 | uniq | bedtools groupby -g 2 -c 2 -o count > ./merged/eRNA.enrichment/all.eRNA.f18.eSNP.gtexDapgDisease.counts
cut -f8,9 ./merged/eRNA.f18.eSNP.gtexCaviarDisease.intersect | sort -k 2,2 | uniq | bedtools groupby -g 2 -c 2 -o count > ./merged/eRNA.enrichment/all.eRNA.f18.eSNP.gtexCaviarDisease.counts
cut -f8,9  ./merged/eRNA.f16.GWASDisease.intersect | sort -k 2,2 | uniq | bedtools groupby -g 2 -c 2 -o count | sort -k 1 | bedtools groupby -g 1 -c 2 -o sum > ./merged/eRNA.enrichment/all.eRNA.f16.GWASDisease.counts

scp rw552@erisone.partners.org:/data/bioinformatics/projects/donglab/AMPPD_eRNA/output/merged/eRNA.enrichment/all.eRNA.'*'.counts .

wd=/Users/rosanwang/Documents/college/dong_lab/code/enHunter/input_files/characterization/feature.enrichment/counts
Rscript fisher_test.R $wd/merged/all.eRNA.f18.eSNP.gtexDapgDisease.counts $wd/snp_gtexDapg.counts gtexDapg_all.pdf "gtexDapg all"
Rscript fisher_test.R $wd/merged/all.eRNA.f18.eSNP.gtexCaviarDisease.counts $wd/snp_gtexCaviar.counts gtexCaviar_all.pdf "gtexCaviar all"
Rscript fisher_test.R $wd/merged/all.eRNA.f16.GWASDisease.counts $wd/GWAS_20220810.v1.02.counts.v2 GWAS_all.pdf "GWAS all"


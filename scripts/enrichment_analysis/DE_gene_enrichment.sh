

pipeline_path=/data/bioinformatics/projects/donglab/AMPPD_eRNA
source $pipeline_path/config.txt

STRAND=$1

EXTERNAL_FEATURE=/data/bioinformatics/external_data/externalData
inputbed=$pipeline_path/inputs/$STRAND/eRNA.bed

# need a new bedtools version to run 
module load bedtools/2.20.1 

# reports all features from inputbed even if the feature doesn't overlap with any genes
intersectBed -a $inputbed -b $GENOME/Annotation/Genes/gencode.v37.annotation.genes.bed -wao

cut -f1 DEresult.padj05.xls | sed '1d' > DEeRNAs.txt

cd /data/bioinformatics/projects/donglab/AMPPD_eRNA/output 
mkdir merged/DE_analysis 
cat <(awk 'OFS="\t" {print $1"_plus", $2}' plus/eRNA.plus.f19.Hostgene.txt) <(awk 'OFS="\t" {print $1"_minus", $2}' minus/eRNA.minus.f19.Hostgene.txt) > merged/DE_analysis/eRNA.f19.Hostgene.txt

grep -F -f /data/bioinformatics/projects/donglab/AMPPD_eRNA/DE_eRNA/PPMI/All/DEeRNAs.txt merged/DE_analysis/eRNA.f19.Hostgene.txt > DEeRNAs.f19.Hostgene.txt 
# note: 
# KI270718.1_2810_3380_minus and KI270718.1_290_830_minus from DEeRNAs.txt are not included in DEeRNAs.f19.Hostgene.txt

grep -v -F -f /data/bioinformatics/projects/donglab/AMPPD_eRNA/DE_eRNA/PPMI/All/DEeRNAs.txt merged/DE_analysis/eRNA.f19.Hostgene.txt > nonDEeRNAs.f19.Hostgene.txt 

module load bedtools/2.20.1 
sort -k2,2 DEeRNAs.f19.Hostgene.txt | bedtools groupby -g 2 -c 1 -o count > DEeRNAs.f19.Hostgene.counts
sort -k2,2 nonDEeRNAs.f19.Hostgene.txt | bedtools groupby -g 2 -c 1 -o count > nonDEeRNAs.f19.Hostgene.counts

# REDOING IT :D
########### 
pipeline_path=/data/bioinformatics/projects/donglab/AMPPD_eRNA
source $pipeline_path/config.txt

cut -f1 DEresult.padj05.xls | sed '1d' > DEgenes.txt 
zcat DEresult.all.xls.gz | cut -f1 | sed '1d' |grep -v -F -f DEgenes.txt > nonDEgenes.txt

grep -F -f DEgenes.txt $GENOME/Annotation/Genes/gencode.v37.annotation.genes.bed > DEgenes.bed
grep -F -f nonDEgenes.txt $GENOME/Annotation/Genes/gencode.v37.annotation.genes.bed > nonDEgenes.bed
# grep commands line above was giving me the wrong number of lines 
grep -F -f nonDEgenes.txt $GENOME/Annotation/Genes/gencode.v37.annotation.genes.bed > nonDEgenes.bed

# testing...  
grep -F -v -f <(cut -f4 DEgenes.bed | awk 'FS="___" {print $2}' | sort | sed '1d') DEgenes.txt | sort
cut -f4 gencode.v37.annotation.genes.bed | awk 'FS="___" {print $2}' | sort 

# seems like it is because the version of the DE genes is older than the gene in gencode.v37.annotation.genes.bed
cut -f1 -d . DEgenes.txt > DEgenes_no_versions.txt

grep -F -f DEgenes_no_versions.txt $GENOME/Annotation/Genes/gencode.v37.annotation.genes.bed > DEgenes_noversion.bed

## still 3 DE genes aren't there .. . 
grep -F -v -f <(cut -f4 DEgenes_noversion.bed | awk 'FS="___" {print $2}' | sort | sed '1d' | cut -f1 -d .) DEgenes_no_versions.txt | sort

#ENSG00000230839 lncRNA
#ENSG00000270544 lincRNA
#ENSG00000274068 

# ENSG00000169084 is repeated twice 

grep -F -f nonDEgenes_no_versions.txt $GENOME/Annotation/Genes/gencode.v37.annotation.genes.bed > nonDEgenes_noversion.bed
# 34846

# should be 35026 -> many missing, i think they are lncRNA

# making bed file of DE eRNAs 
awk 'OFS="\t" {split($0, a, "_"); print a[1], a[2], a[3], $1}' DEeRNAs.txt  > DEeRNAs.bed

zcat DEresult.all.xls.gz | cut -f1 | grep -v -F -f DEeRNAs.txt > nonDEeRNAs.txt
awk 'OFS="\t" {split($0, a, "_"); print a[1], a[2], a[3], $1}' nonDEeRNAs.txt  > nonDEeRNAs.bed

### note I don't believe Ruifeng ran DE analysis on all genes? 


####### anyways..... going to run bedtools intersect 
cd /data/bioinformatics/projects/donglab/AMPPD_eRNA/output/merged/DE_analysis

# AB: DE genes with DE eRNAs inside & nAB: DE genes without DE eRNAs inside 
bedtools intersect -a /data/bioinformatics/projects/donglab/AMPPD_eRNA/DE_genes/PPMI/DEgenes_noversion.bed  -b sorted_DEeRNAs.bed -c > DE_gene_DE_eRNA.bed

# AnB: not DE gene with DE eRNAs inside & nAnB: not DE gene without DE eRNAs inside 
bedtools intersect -a /data/bioinformatics/projects/donglab/AMPPD_eRNA/DE_genes/PPMI/nonDEgenes_noversion.bed  -b sorted_DEeRNAs.bed -c > nonDE_gene_DE_eRNA.bed

# column 7 is the count of TNEs 

> enrichment.txt 

# AB: 
awk '$7 > 0' DE_gene_DE_eRNA.bed | wc -l >> enrichment.txt 

# AnB 
awk '$7 == 0' DE_gene_DE_eRNA.bed | wc -l >> enrichment.txt 

# nAB: 
awk '$7 > 0' nonDE_gene_DE_eRNA.bed | wc -l >> enrichment.txt 

# nAnB 
awk '$7 == 0' nonDE_gene_DE_eRNA.bed | wc -l >> enrichment.txt 


### checking LITAF 

bedtools intersect -a /data/bioinformatics/projects/donglab/AMPPD_eRNA/DE_genes/PPMI/DEgenes_noversion.bed  -b sorted_DEeRNAs.bed -wb | grep HGNC:16841___ENSG00000189067.14___LITAF 

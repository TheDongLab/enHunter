# creating merge table for ALL eQTL and GWAS snp information 

# usage: Rscript TNE-eQTL-SNP-TFBS-table.R caviar.table dapg.table
# eRNA.plus.f18.eSNP.gtexDapgDisease.intersect eRNA.minus.f18.eSNP.gtexDapgDisease.intersect eRNA.plus.f18.eSNP.gtexCaviarDisease.intersect eRNA.minus.f18.eSNP.gtexCaviarDisease.intersect eRNA.plus.f16.GWASDisease.intersect eRNA.minus.f16.GWASDisease.intersect 

args<-commandArgs(trailingOnly=TRUE)

#plus strand Dapg SNPs
system(paste0("awk 'OFS=\"", "\t\"", " {print $4, \"+\", $8, $5\"_\"$6\"_\"$7, $9, \"Yes\"}' eRNA.plus.f18.eSNP.gtexDapgDisease.intersect > eRNA.f18.plus.eSNP.gtexDapg.txt"))

# minus strand Dapg SNPs
system(paste0("awk 'OFS=\"\t\" {print $4, \"-\", $8, $5\"_\"$6\"_\"$7, $9, \"Yes\"}' eRNA.minus.f18.eSNP.gtexDapgDisease.intersect > eRNA.f18.minus.eSNP.gtexDapg.txt"))

system(paste0("cat eRNA.f18.plus.eSNP.gtexDapg.txt eRNA.f18.minus.eSNP.gtexDapg.txt > eRNA.f18.eSNP.gtexDapg.txt"))

#plus strand Caviar SNPs
system(paste0("awk 'OFS=\"\t\" {print $4, \"+\", $8, $5\"_\"$6\"_\"$7, $9, \"Yes\"}' eRNA.plus.f18.eSNP.gtexCaviarDisease.intersect > eRNA.f18.plus.eSNP.gtexCaviar.txt") )

# minus strand Caviar SNPs
system(paste0("awk 'OFS=\"\t\" {print $4, \"-\", $8, $5\"_\"$6\"_\"$7, $9, \"Yes\"}' eRNA.minus.f18.eSNP.gtexCaviarDisease.intersect > eRNA.f18.minus.eSNP.gtexCaviar.txt") )

system(paste0("cat eRNA.f18.plus.eSNP.gtexCaviar.txt eRNA.f18.minus.eSNP.gtexCaviar.txt > eRNA.f18.eSNP.gtexCaviar.txt") )

# read in the Caviar and Dapg table 
# V19: tissue 
# V1: chromosome 
# V7: eQTL start 
# V8: eQTL end 
# V14: rsid 
# V16: gene name 
caviar <- read.table(args[1])
caviar <- caviar[, c(1,7,8,14,16,19)]
colnames(caviar) <- c("chr", "start", "end", "rsid", "gene", "tissue")
caviar$eQTL_pos <- paste(caviar$chr, caviar$start, caviar$end, sep="_")
#caviar[,-c("chr", "start", "end")]
caviar[,-c(1:3)]

dapg <- read.table(args[2])
dapg <- dapg[, c(1,7,8,14,16,19)]
colnames(dapg) <- c("chr", "start", "end", "rsid", "gene", "tissue")
dapg$eQTL_pos <- paste(dapg$chr, dapg$start, dapg$end, sep="_")
#dapg <- dapg[, -c("chr", "start", "end")]
dapg <- dapg[, -c(1:3)]

TNE_Dapg <- read.table("eRNA.f18.eSNP.gtexDapg.txt")
colnames(TNE_Dapg) <- c("TNE", "strand", "rsid", "eQTL_pos", "tissue", "dapg")

TNE_Dapg <- merge(TNE_Dapg, dapg, by=c("eQTL_pos", "tissue", "rsid"), all.x = TRUE)

TNE_Caviar <- read.table("eRNA.f18.eSNP.gtexCaviar.txt")
colnames(TNE_Caviar) <- c("TNE", "strand", "rsid", "eQTL_pos", "tissue", "caviar")

TNE_Dapg <- merge(TNE_Caviar , caviar, by=c("eQTL_pos", "tissue", "rsid"), all.x = TRUE)

all_eQTL <- merge(TNE_Dapg, TNE_Caviar, by=c("TNE", "strand", "rsid" ,"eQTL_pos", "tissue"), all=TRUE)

write.table(all_eQTL, "eRNA.eQTL.snps", sep="\t", quote = F, col.names = T, row.names = F)

## combining the GWAS information 
system(paste0("awk 'OFS=\"\t\"; FS=\"\t\" {print $4, \"+\", $8, $5\"_\"$6\"_\"$7, $9}' eRNA.plus.f16.GWASDisease.intersect > GWAS_plus.txt"))
system(paste0("awk 'OFS=\"\t\"; FS=\"\t\" {print $4, \"-\", $8, $5\"_\"$6\"_\"$7, $9}'  eRNA.minus.f16.GWASDisease.intersect > GWAS_minus.txt"))

system("cat GWAS_minus.txt GWAS_plus.txt > GWAS.txt")

GWAS <- read.table("GWAS.txt")
colnames(GWAS) <- c("TNE", "strand", "rsid", "eQTL_pos", "disease")

eQTL_GWAS <- merge(all_eQTL, GWAS, sep="\t", all=T)

write.table(eQTL_GWAS, "eRNA.eQTL.snps.GWAS.xls", sep="\t", quote = F, col.names = T, row.names = F)
### TODO update this script with script on erisone 

##### adding TFBS table information 

# processing the TFBS intersect file 
# 1,2,3 : binding site of TF in the TNE 
# 4: TNE name
# 8: TF name
#awk 'OFS="\t" {print $1, $2, $3, $4, $8}' eRNA.minus.f06.TFBS.intersect > eRNA.minus.TFBS.bed 
#awk 'OFS="\t" {print $1, $2, $3, $4, $8}' eRNA.plus.f06.TFBS.intersect > eRNA.plus.TFBS.bed 

# getting the SNP positions from eQTL_GWAS in a bed file 
# field 4 : eQTL chrom position 
# cut -f 4 eRNA.eQTL.snps.GWAS.xls |sed '1d' | awk 'BEGIN{FS="_"; OFS="\t"} {print $1, $2, $3, $0}' | sort -k1,1 -k2,2n | uniq > eQTL_GWAS_snps.bed

# module load bedtools/2.26.0
#bedtools intersect -a /data/bioinformatics/projects/donglab/AMPPD_eRNA/output/table/eQTL_GWAS_snps.bed -b /data/bioinformatics/projects/donglab/AMPPD_eRNA/output/minus/eRNA.TFBS/eRNA.minus.TFBS.bed  -wa -wb > eRNA.minus.TF_SNP.intersect 
#bedtools intersect -a /data/bioinformatics/projects/donglab/AMPPD_eRNA/output/table/eQTL_GWAS_snps.bed -b /data/bioinformatics/projects/donglab/AMPPD_eRNA/output/plus/eRNA.TFBS/eRNA.plus.TFBS.bed  -wa -wb > eRNA.plus.TF_SNP.intersect 

##### output format 
# 1,2,3: SNP coordinates 
# 4: SNP eQTL pos 
# 5,6,7: TNE TFBS coordinates 
# 8 : TNE name 
# 9: TF name 

#module load bedtools/2.20.1 
# cat eRNA.minus.TF_SNP.intersect eRNA.plus.TF_SNP.intersect | cut -f4,9 | sort -k4,4 | groupBy -g 1 -c 2 -o collapse > grouped_eRNA.TF_SNPs.txt 

library(data.table)
TFBS_TNEs <- fread("grouped_eRNA.TF_SNPs.txt", header = FALSE)
colnames(TFBS_TNEs) <- c("eQTL_pos", "TF")

eQTL_GWAS <- fread("eRNA.eQTL.snps.GWAS.xls")
eQTL_GWAS_TFBS <- merge(TFBS_TNEs, eQTL_GWAS, by="eQTL_pos", all.y=T)

fwrite(eQTL_GWAS_TFBS, "eQTL_GWAS_TFBS.xls", sep="\t", quote = F, col.names = T, row.names = F)

### add TNE class information to table 
library(data.table)
library(dplyr)
eQTL_GWAS_TFBS <- fread("eQTL_GWAS_TFBS.xls")
eRNA_class <- fread("./input_files/characterization/eRNA.characterize.feature.color.xls")  %>% select(V1, class, strand)
colnames(eRNA_class) <- c("TNE", "class", "strand")
eRNA_class_eQTL_GWAS_TFBS <- merge(eQTL_GWAS_TFBS, eRNA_class, by=c("TNE", "strand"), all.x=TRUE)

write.table(eRNA_class_eQTL_GWAS_TFBS, "eRNA.class.eQTL.snps.GWAS.TFBS.xls", sep="\t", quote = F, col.names = T, row.names = F)

### add DE gene information 
library(data.table)
gencode_genes <- fread("/data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs/gencode.genes.no_version.txt")
colnames(gencode_genes) <- c("ENSG", "Gene")

eRNA_class_eQTL_GWAS_TFBS <- fread("eRNA.class.eQTL.snps.GWAS.TFBS.xls")

merge(eRNA_class_eQTL_GWAS_TFBS, gencode_genes)




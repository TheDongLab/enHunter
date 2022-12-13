# creating merge table for ALL eQTL and GWAS snp information 

# usage: Rscript TNE-eQTL-SNP-TFBS-table.R caviar.table dapg.table

# eRNA.plus.f18.eSNP.gtexDapgDisease.intersect eRNA.minus.f18.eSNP.gtexDapgDisease.intersect eRNA.plus.f18.eSNP.gtexCaviarDisease.intersect eRNA.minus.f18.eSNP.gtexCaviarDisease.intersect eRNA.plus.f16.GWASDisease.intersect eRNA.minus.f16.GWASDisease.intersect 
# module load R/4.0.2 

#args<-commandArgs(trailingOnly=TRUE)

#plus strand Dapg SNPs
#system(paste0("awk 'OFS=\"", "\t\"", " {print $4, \"+\", $8, $5\"_\"$6\"_\"$7, $9, \"Yes\"}' eRNA.plus.f18.eSNP.gtexDapgDisease.intersect > eRNA.f18.plus.eSNP.gtexDapg.txt"))

# minus strand Dapg SNPs
#system(paste0("awk 'OFS=\"\t\" {print $4, \"-\", $8, $5\"_\"$6\"_\"$7, $9, \"Yes\"}' eRNA.minus.f18.eSNP.gtexDapgDisease.intersect > eRNA.f18.minus.eSNP.gtexDapg.txt"))

#system(paste0("cat eRNA.f18.plus.eSNP.gtexDapg.txt eRNA.f18.minus.eSNP.gtexDapg.txt > eRNA.f18.eSNP.gtexDapg.txt"))

#plus strand Caviar SNPs
#system(paste0("awk 'OFS=\"\t\" {print $4, \"+\", $8, $5\"_\"$6\"_\"$7, $9, \"Yes\"}' eRNA.plus.f18.eSNP.gtexCaviarDisease.intersect > eRNA.f18.plus.eSNP.gtexCaviar.txt") )

# minus strand Caviar SNPs
#system(paste0("awk 'OFS=\"\t\" {print $4, \"-\", $8, $5\"_\"$6\"_\"$7, $9, \"Yes\"}' eRNA.minus.f18.eSNP.gtexCaviarDisease.intersect > eRNA.f18.minus.eSNP.gtexCaviar.txt") )

#system(paste0("cat eRNA.f18.plus.eSNP.gtexCaviar.txt eRNA.f18.minus.eSNP.gtexCaviar.txt > eRNA.f18.eSNP.gtexCaviar.txt") )

library(data.table)
library(dplyr)
# read in the Caviar and Dapg table 
# V19: tissue 
# V1: chromosome 
# V7: eQTL start 
# V8: eQTL end 
# V14: rsid 
# V16: gene name 
print("Reading in caviar data table")
caviar <- fread("/data/bioinformatics/external_data/externalData/GTEx_p_value/GTEx_hg38_UCSC_track/gtexCaviar.bed")
#caviar <- fread(args[1])
caviar <- caviar[, c(1,7,8,14,16,17,19)] ### contains duplicates
colnames(caviar) <- c("chr", "start", "end", "rsid", "caviar_gene", "caviar_ensg", "tissue")
caviar$eQTL_pos <- paste(caviar$chr, caviar$start, caviar$end, sep="_")
#caviar[,-c("chr", "start", "end")]
caviar <- caviar[,-c(1:3)]

print("Reading in dapg data table")
dapg <- fread("/data/bioinformatics/external_data/externalData/GTEx_p_value/GTEx_hg38_UCSC_track/gtexDapg.bed")
dapg <- dapg[, c(1,7,8,14,16,17,19)]
colnames(dapg) <- c("chr", "start", "end", "rsid", "dapg_gene", "dapg_ensg","tissue")
dapg$eQTL_pos <- paste(dapg$chr, dapg$start, dapg$end, sep="_")
#dapg <- dapg[, -c("chr", "start", "end")]
dapg <- dapg[,-c(1:3)]

print("Reading in TNEs")
TNE_Dapg <- fread("eRNA.f18.eSNP.gtexDapg.txt")
colnames(TNE_Dapg) <- c("TNE", "strand", "rsid", "eQTL_pos", "tissue", "dapg")

print("merging Dapg")
# this is the line adding more rows .. because one eQTL_pos/rsid/tissue can have multiple genes
# also, dapg and caviar dataframes contain some duplicate rows ... 
# note: all.x = TRUE and all.x = FALSE does not change the number of rows 
TNE_Dapg <- merge(TNE_Dapg, dapg, by=c("eQTL_pos", "tissue", "rsid"), all.x = TRUE) %>% 
  group_by(eQTL_pos, tissue, rsid, TNE, strand) %>% 
  mutate(dapg_genes = paste0(dapg_gene, collapse=",")) %>% 
  mutate(dapg_ensgs = paste0(dapg_ensg, collapse=",")) %>% 
  distinct(across(contains("dapg_genes")), .keep_all = TRUE)

TNE_Dapg <- TNE_Dapg[, -c(7,8)]

TNE_Caviar <- read.table("eRNA.f18.eSNP.gtexCaviar.txt")
colnames(TNE_Caviar) <- c("TNE", "strand", "rsid", "eQTL_pos", "tissue", "caviar")

print("merging Caviar...")

TNE_Caviar <- merge(TNE_Caviar , caviar, by=c("eQTL_pos", "tissue", "rsid"), all.x = TRUE) %>% 
  group_by(eQTL_pos, tissue, rsid, TNE, strand) %>% 
  mutate(caviar_genes = paste0(caviar_gene, collapse=",")) %>% 
  mutate(caviar_ensgs = paste0(caviar_ensg, collapse=",")) %>% 
  distinct(across(contains("caviar_genes")), .keep_all = TRUE)

# distinct removes the duplicate rows with same caviar_genes values 
# nrow of TNE_Caviar now matches the eRNA.f18.eSNP.gtexCaviar table 

TNE_Caviar <- TNE_Caviar[, -c(7,8)]

# groupby the eQTL_pos, tissue, rsid -> multiple genes into a list 
all_eQTL <- merge(TNE_Dapg, TNE_Caviar, by=c("TNE", "strand", "rsid" ,"eQTL_pos", "tissue"), all=TRUE)
# i believe the only NAs should only be from the the caviar and dapg column 
all_eQTL[is.na(all_eQTL)] <- "NA"

print("writing table")
fwrite(all_eQTL, "eRNA.eQTL.snps.xls", sep="\t", quote = F, col.names = T, row.names = F)

## combining the GWAS information 
#system(paste0("awk 'OFS=\"\t\"; FS=\"\t\" {print $4, \"+\", $8, $5\"_\"$6\"_\"$7, $9}' eRNA.plus.f16.GWASDisease.intersect > GWAS_plus.txt"))
#system(paste0("awk 'OFS=\"\t\"; FS=\"\t\" {print $4, \"-\", $8, $5\"_\"$6\"_\"$7, $9}'  eRNA.minus.f16.GWASDisease.intersect > GWAS_minus.txt"))


#system("cat GWAS_minus.txt GWAS_plus.txt > GWAS.txt")

GWAS <- fread("GWAS.txt", header = FALSE)
colnames(GWAS) <- c("TNE", "strand", "rsid", "eQTL_pos", "disease")

eQTL_GWAS <- merge(all_eQTL, GWAS, sep="\t", all=T)

write.table(eQTL_GWAS, "eRNA.eQTL.snps.GWAS.xls", sep="\t", quote = F, col.names = T, row.names = F)

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

# module load bedtools/2.20.1 
# cat eRNA.minus.TF_SNP.intersect eRNA.plus.TF_SNP.intersect | cut -f4,9 | sort -k4,4 | groupBy -g 1 -c 2 -o collapse > grouped_eRNA.TF_SNPs.txt 

library(data.table)
TFBS_TNEs <- fread("grouped_eRNA.TF_SNPs.txt", header = FALSE)
colnames(TFBS_TNEs) <- c("eQTL_pos", "TF")

eQTL_GWAS <- fread("eRNA.eQTL.snps.GWAS.xls")
eQTL_GWAS_TFBS <- merge(TFBS_TNEs, eQTL_GWAS, by="eQTL_pos", all.y=T)

fwrite(eQTL_GWAS_TFBS, "eQTL_GWAS_TFBS.xls", sep="\t", quote = F, col.names = T, row.names = F)

### add DE TNE information 
library(stringr)
DE_TNE <- fread("/data/bioinformatics/projects/donglab/AMPPD_eRNA/DE_eRNA/PPMI/All/DEeRNAs.txt", header = F)
DE_TNE <- cbind(DE_TNE, str_split_fixed(DE_TNE$V1, "_", 4))
colnames(DE_TNE) <- c("TNE", "chr", "start", "end", "strand")
DE_TNE$strand[DE_TNE$strand == "minus"] <- "-"
DE_TNE$strand[DE_TNE$strand == "plus"] <- "+"
DE_TNE <- DE_TNE %>% mutate(TNE = paste(chr, start, end, sep = "_"))

DE_TNE$DE_TNE <- TRUE
DE_TNE <- DE_TNE[,-c(2:4)]

new_test <- merge(eQTL_GWAS_TFBS, DE_TNE, by = c("TNE", "strand"), all.x = TRUE)

write.table(new_test, "eRNA_class_eQTL_GWAS_TFBS_DETNE.xls", sep="\t", quote = F, col.names = T, row.names = F)

### add TNE class and host_gene information to table 
library(data.table)
library(dplyr)
eQTL_GWAS_TFBS <- fread("eRNA_class_eQTL_GWAS_TFBS_DETNE.xls")
eRNA_class <- fread("./input_files/characterization/eRNA.characterize.feature.color.xls")  %>% select(V1, class, strand, f19.Hostgene)
colnames(eRNA_class) <- c("TNE", "class", "strand", "hostgene")
eRNA_class_eQTL_GWAS_TFBS <- merge(eQTL_GWAS_TFBS, eRNA_class, by=c("TNE", "strand"), all.x=TRUE)

write.table(eRNA_class_eQTL_GWAS_TFBS, "eRNA.DETNE.class.hostgene.eQTL.snps.GWAS.TFBS.xls", sep="\t", quote = F, col.names = T, row.names = F)

### add PCHi-C data 

# eRNA.minus.f22.TNE.PCHiCPromoters format 
# 1,2,3 : TNE coordinates 
# 4 : TNE name 
# 5,6,7 : OE coordinates 
# 8 : OE ID 
# OE : promoter targets 

# cat /data/bioinformatics/projects/donglab/AMPPD_eRNA/output/minus/eRNA.PCHiC/eRNA.minus.f22.TNE.PCHiCPromoters /data/bioinformatics/projects/donglab/AMPPD_eRNA/output/plus/eRNA.PCHiC/eRNA.plus.f22.TNE.PCHiCPromoters | cut -f4-9 > eRNA.TNE.PCHiCPromoters
# change this input file so that it is grouped by TNE (i.e. combine multiple lines of promoter together)
# module load bedtools/2.20.1  
# sort -k1,1 eRNA.TNE.PCHiCPromoters | groupBy -g 1 -c 6 -o distinct > eRNA.TNE.PCHiCPromoters.distinct.txt
library(data.table)
library(dplyr)
HiC <- fread("eRNA.TNE.PCHiCPromoters.distinct.txt", header = FALSE)
colnames(HiC) <- c("TNE", "Promoters")
eRNA_class_eQTL_GWAS_TFBS <- fread("eRNA.DETNE.class.hostgene.eQTL.snps.GWAS.TFBS.xls")
test <- merge(eRNA_class_eQTL_GWAS_TFBS, HiC, by=c("TNE"), all.x = TRUE)

write.table(test, "eRNA_class_eQTL_GWAS_TFBS_PCHiC.xls", sep="\t", quote = F, col.names = T, row.names = F)

## doing some more clean up with the promoter Hi-C table 
PCHiC_table <- fread("eRNA_class_eQTL_GWAS_TFBS_PCHiC.xls")
PCHiC_table$Prom <- lapply(PCHiC_table$Promoters, function(x) {
  strsplit(x, ";|,") %>% unlist() %>% unique() %>% paste0(collapse=",") 
  })

PCHiC_table <- PCHiC_table %>% select(-Promoters)

fwrite(PCHiC_table, "eRNA_class_eQTL_GWAS_TFBS_PCHiC_v2.xls", sep="\t", quote = F, col.names = T, row.names = F)

### add DE gene information 
# this will be tricky 
library(data.table)
gencode_genes <- fread("gencode.genes.no_version.txt", header= F)
colnames(gencode_genes) <- c("ENSG", "Gene")

v1_table <- fread("eRNA_class_eQTL_GWAS_TFBS_PCHiC_v2.xls")
DE_genes <- fread("DEgenes_no_versions.txt", header = FALSE)

## this is so slow ... 
Prom_ensgs <- lapply(unique(v1_table$Prom), function(x){
  genes <- unlist(strsplit(x, ",", fixed=T))
  ensg <- gencode_genes %>% filter(gencode_genes$Gene %in% genes) %>% select(ENSG) %>% unlist()
  list(Prom=x, Prom_ensgs=paste0(ensg, collapse=",") )
})

prom_ensgs_table <- rbindlist(Prom_ensgs)

test <- merge(v1_table, prom_ensgs_table, by=c("Prom"), all.x = TRUE) %>% 
  select(TNE, strand, DE_TNE, class, eQTL_pos, rsid, tissue, dapg_genes, 
         dapg_ensgs, caviar_genes, caviar_ensgs, TF, disease, hostgene, Prom, Prom_ensgs)

# make temp column of all the genes (dapg, caviar and HiC) in ENSG format
new_df <- test  
new_df$DE_ensgs <- apply(new_df, 1, function(x) {
  all <- c(unlist(sub("[.][0-9]+$", "", unlist(strsplit(x[["dapg_ensgs"]], ",", fixed=TRUE )))), 
           unlist(sub("[.][0-9]+$", "", unlist(strsplit(x[["caviar_ensgs"]], ",", fixed=TRUE )))), 
           unlist(sub("[.][0-9]+$", "", unlist(strsplit(x[["Prom_ensgs"]], ",", fixed=TRUE )))))
  all <- all[all != ""]
  DE_all <- paste0(unlist(DE_genes[DE_genes$V1 %in% all]), collapse=",")
  DE_all
})

new_df$DE_genes <- lapply(new_df$DE_ensgs,function(x)  {
  DE_all <- unlist(strsplit(x, ",", fixed = TRUE))
  gencode_genes[gencode_genes$ENSG %in% DE_all]$Gene 
})

fwrite(new_df, "FINAL.xls", sep="\t", quote = F, col.names = T, row.names = F)



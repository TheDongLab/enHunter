library(data.table)
library(dplyr)

setwd("~/Documents/college/dong_lab/code/playground/GWAS-eQTL")

gwas <- fread("eRNA.minus.f16.GWASDisease.txt")
gwas <- gwas[,4:9]
colnames(gwas) <- c("TNE", "chr", "start", "end", "rsid", "disease")

eqtl <- fread("eRNA.minus.f18.eSNP.gtexCaviarDisease.txt")
eqtl <- eqtl[,4:9]
colnames(eqtl) <- c("TNE", "chr", "start", "end", "rsid", "tissue")

# TODO add class information of the TNE 
all_merged_snps <- merge(gwas, eqtl, by=c("chr", "start", "end", "rsid", "TNE"))
blood_snps <- all_merged_snps[tissue == "Whole_Blood"]

# read in the Caviar table 
# V19: tissue 
# V1: chromosome 
# V7: eQTL start 
# V8: eQTL end 
# V14: rsid 
# V16: gene name 
caviar <- fread("gtexCaviar.bed")
caviar <- caviar[, c(1,7,8,14,16,19)]
colnames(caviar) <- c("chr", "start", "end", "rsid", "gene", "tissue")

blood_snps_genes <- merge(blood_snps, caviar, by=c("chr", "start", "end", "rsid", "tissue"))

## checking for eqtl snps which have multiple genes 
dups <- duplicated(blood_snps_genes[,1:7])
multi_genes <- blood_snps_genes[dups]

write.csv(blood_snps_genes %>% select(rsid, TNE, disease, gene), 
            file="blood_gene_snps.csv", row.names=FALSE)

# RERE
# STRIP1
# JAKMIP3
# C1QTNF4
# RAB30-AS1
# RP11-85K15.3
# SERPINA1
# SLC12A1
# DUT 
# CHSY1
# CMIP
# RP11-391L3.5
# NATD1
# TRAF4 
# LRRC37A
# RP11-798G7.6
# ARL17A
# HOXB4 





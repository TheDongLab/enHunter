library(data.table)
library(dplyr)

gwas <- fread("./input_files/characterization/feature.enrichment/eRNA.minus.f16.GWASDisease.txt")
gwas <- gwas[,4:9]
colnames(gwas) <- c("TNE", "chr", "start", "end", "rsid", "disease")

eqtl <- fread("./input_files/characterization/feature.enrichment/eRNA.minus.f18.eSNP.gtexCaviarDisease.txt")
eqtl <- eqtl[,4:9]
colnames(eqtl) <- c("TNE", "chr", "start", "end", "rsid", "tissue")

all_merged_snps <- merge(gwas, eqtl, by=c("chr", "start", "end", "rsid", "TNE"))
blood_snps <- all_merged_snps[tissue == "Whole_Blood"]

# read in the Caviar table 
# V19: tissue 
# V1: chromosome 
# V7: eQTL start 
# V8: eQTL end 
# V14: rsid 
# V16: gene name 
caviar <- fread("./input_files/characterization/externalData/gtexCaviar.bed")
caviar <- caviar[, c(1,7,8,14,16,19)]
colnames(caviar) <- c("chr", "start", "end", "rsid", "gene", "tissue")

blood_snps_genes <- merge(blood_snps, caviar, by=c("chr", "start", "end", "rsid", "tissue"))

## checking for eqtl snps which have multiple genes 
dups <- duplicated(blood_snps_genes[,1:7])
multi_genes <- blood_snps_genes[dups]

# add class information of the TNE 
eRNA_class <- fread("./input_files/characterization/eRNA.minus.characterize.xls")  %>% select(V1, class)
blood_snps_genes <- merge(blood_snps_genes, eRNA_class, by.x="TNE", by.y ="V1")

write.csv(blood_snps_genes %>% select(rsid, TNE, disease, gene), 
            file="blood_gene_snps.csv", row.names=FALSE)

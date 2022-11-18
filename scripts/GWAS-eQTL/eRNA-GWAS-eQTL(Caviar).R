library(data.table)
library(dplyr)

gwas <- fread("./input_files/characterization/feature.enrichment/eRNA.minus.f16.GWASDisease.intersect")
gwas <- gwas[,4:9]
colnames(gwas) <- c("TNE", "chr", "start", "end", "rsid", "disease")

eqtl <- fread("./input_files/characterization/feature.enrichment/eRNA.minus.f18.eSNP.gtexCaviarDisease.intersect")
eqtl <- eqtl[,4:9]
colnames(eqtl) <- c("TNE", "chr", "start", "end", "rsid", "tissue")
# the start and end coordinates are of the SNP 

all_merged_snps <- merge(gwas, eqtl, by=c("chr", "start", "end", "rsid", "TNE"))
#blood_snps <- all_merged_snps[tissue == "Whole_Blood"]

# read in the Caviar table 
# V19: tissue 
# V1: chromosome 
# V7: eQTL start 
# V8: eQTL end 
# V14: rsid 
# V16: gene name 
caviar <- fread("./input_files/characterization/externalData/gtexCaviar.bed")
caviar <- caviar[, c(1,7,8,13,14,16,19)]
colnames(caviar) <- c("chr", "start", "end", "eqtlPos", "rsid", "gene", "tissue")

# need to get the gene information 
snps_genes <- merge(all_merged_snps, caviar, by=c("chr", "start", "end", "rsid", "tissue"))

# some SNPs have multiple genes, causing the TNE-SNP pair to appear twice 
## checking for eqtl snps which have multiple genes 
#dups <- duplicated(snps_genes[,1:7])
#multi_genes <- snps_genes[dups]

# add class information of the TNE 
eRNA_class <- fread("./input_files/characterization/eRNA.minus.characterize.xls")  %>% select(V1, class)
snps_genes <- merge(snps_genes, eRNA_class, by.x="TNE", by.y ="V1")

write.csv(snps_genes %>% select(-c("chr", "start", "end")), file="blood_gene_snps.csv", row.names=FALSE)

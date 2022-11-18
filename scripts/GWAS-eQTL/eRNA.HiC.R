library(data.table)
library(dplyr)

HiC <- fread("./input_files/characterization/feature.enrichment/eRNA.minus.f22.TNE.PCHiCPromoters")[, 4:9]
colnames(HiC) <- c("TNE", "oeChr", "oeStart", "oeEnd", "oeID", "Promoter")

GWAS_eQTLs <- fread("./scripts/GWAS-eQTL/blood_gene_snps.csv")

test <- merge(GWAS_eQTLs, HiC, by="TNE")

write.csv(test, file="HiC_eQTL_GWAS_snps.csv", row.names=FALSE)

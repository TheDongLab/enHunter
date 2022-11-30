library(data.table)
library(dplyr)

HiC <- fread("./input_files/characterization/feature.enrichment/eRNA.minus.f22.TNE.PCHiCPromoters")[, 4:9]
colnames(HiC) <- c("TNE", "oeChr", "oeStart", "oeEnd", "oeID", "Promoter")

GWAS_eQTLs <- fread("./scripts/GWAS-eQTL/blood_gene_snps.csv")

test <- merge(GWAS_eQTLs, HiC, by="TNE")
test2 <- merge(GWAS_eQTLs, HiC, by="TNE", all.x = TRUE)

test2$PCHiC_eQTL_gene <- apply(test2, 1, function(x) {
  grepl(x["gene"], x["Promoter"])
})

write.csv(test, file="HiC_eQTL_GWAS_snps.csv", row.names=FALSE)
write.csv(test2, file="HiC_eQTL_GWAS_snps_v2.csv", row.names=FALSE)



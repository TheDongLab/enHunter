# usage: R script to create the following table
####### table output #######
# columns
# 1. TNE
# 2. target_gene 
# 3. source of target_gene evidence 
# 4. pval for evidence association 
# 5. host gene

library(data.table)
library(dplyr)

TNE_hostgene <- fread("/data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs/eRNA.characterize.feature.color.xls") %>% 
  select(V1, strand, f19.Hostgene)

colnames(TNE_hostgene) <- c("TNE", "strand", "Hostgene")

#### adding eQTL dapg and caviar information ####
# read in the Caviar and Dapg table 
# V1: chromosome 
# V7: eQTL start 
# V8: eQTL end 
# V14: rsid 
# V16: gene name 
# V19: tissue 
print("Reading in caviar data table")
caviar <- fread("/data/bioinformatics/external_data/externalData/GTEx_p_value/GTEx_hg38_UCSC_track/gtexCaviar.bed")
caviar <- caviar[, c(1,7,8,14,16,17,19,20)] ### contains duplicates
colnames(caviar) <- c("chr", "start", "end", "rsid", "caviar_gene", "caviar_ensg", "tissue", "cpp")
caviar$eQTL_pos <- paste(caviar$chr, caviar$start, caviar$end, sep="_")
#caviar[,-c("chr", "start", "end")]
caviar <- caviar[,-c(1:3)]

#"eRNA.f18.eSNP.gtexCaviar.txt": contains all the SNPs located within TNEs 
TNE_Caviar <- read.table("/data/bioinformatics/projects/donglab/AMPPD_eRNA/output/snp_table/eRNA.f18.eSNP.gtexCaviar.txt")
colnames(TNE_Caviar) <- c("TNE", "strand", "rsid", "eQTL_pos", "tissue", "caviar")

print("merging Caviar...")

# merge caviar table (with gene and cpp information) with TNE caviar SNP information
TNE_Caviar <- merge(TNE_Caviar , caviar, by=c("eQTL_pos", "tissue", "rsid"), all.x = TRUE) %>% 
  group_by(eQTL_pos, tissue, rsid, TNE, strand) %>% 
  mutate(caviar_genes = paste0(caviar_gene, collapse=",")) %>% 
  mutate(caviar_ensgs = paste0(caviar_ensg, collapse=",")) %>% 
  distinct(across(contains("caviar_genes")), .keep_all = TRUE)

# distinct removes the duplicate rows with same caviar_genes values 
# nrow of TNE_Caviar now matches the eRNA.f18.eSNP.gtexCaviar table 

TNE_Caviar <- TNE_Caviar[, -c(7,8)]
TNE_Caviar$source <- "caviar"

### merge with hostgene information table 
TNEs <- merge(TNE_Caviar, TNE_hostgene, by=c("TNE", "strand"), all.y = T)
TNEs <- TNEs %>% select(-eQTL_pos, -rsid, -caviar) %>% rename(gene = caviar_genes, ensgs = caviar_ensgs)


print("Reading in dapg data table")
dapg <- fread("/data/bioinformatics/external_data/externalData/GTEx_p_value/GTEx_hg38_UCSC_track/gtexDapg.bed")
dapg <- dapg[, c(1,7,8,14,16,17,19,20)]
colnames(dapg) <- c("chr", "start", "end", "rsid", "dapg_gene", "dapg_ensg","tissue", "pip")
dapg$eQTL_pos <- paste(dapg$chr, dapg$start, dapg$end, sep="_")
dapg <- dapg[,-c(1:3)]

print("Reading in TNEs")
TNE_Dapg <- fread("/data/bioinformatics/projects/donglab/AMPPD_eRNA/output/snp_table/eRNA.f18.eSNP.gtexDapg.txt")
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

TNE_Dapg <- TNE_Dapg[, -c(7,8)] %>% select(-eQTL_pos, -rsid, -dapg) %>% rename(gene = dapg_genes, ensgs = dapg_ensgs)
TNE_Dapg$source <- "dapg"

### merge with hostgene information table 
TNEs <- rbind(TNE_Dapg, TNEs, by=c("TNE", "strand"))

print("writing table")
fwrite(TNEs, "TNEs.test.xls", sep="\t", quote = F, col.names = T, row.names = F)


# groupby the eQTL_pos, tissue, rsid -> multiple genes into a list 
# i believe the only NAs should only be from the the caviar and dapg column 
# all_eQTL[is.na(all_eQTL)] <- "NA"
# 
# print("writing table")



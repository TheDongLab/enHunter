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
TNE_Caviar <- TNE_Caviar %>% ungroup() %>% select(-c(eQTL_pos, rsid, caviar)) %>% rename(gene = caviar_genes, ensgs = caviar_ensgs)

### merge with hostgene information table 
## TODO I think this line is adding extra empty lines 
# TNEs <- merge(TNE_Caviar, TNE_hostgene, by=c("TNE", "strand"), all.y = T)

#TNEs <- TNEs %>% select(-eQTL_pos, -rsid, -caviar) %>% rename(gene = caviar_genes, ensgs = caviar_ensgs)

###### TNEs ######
# 1. TNE 
# 2. strand 
# 3. tissue 
# 4. cpp 
# 5. gene 
# 6. ensgs 
# 7. source 
# 8. Hostgene 

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

TNE_Dapg <- TNE_Dapg %>% ungroup() %>% select(-c(dapg_gene, dapg_ensg, eQTL_pos, rsid, dapg) ) %>% 
  rename(gene = dapg_genes, ensgs = dapg_ensgs, cpp = pip)
TNE_Dapg$source <- "dapg"

### merge with hostgene information table 
TNE_eQTLs <- rbind(TNE_Dapg, TNE_Caviar)

## TODO make sure to only merge with host gene at the END of getting all the individual target gene TNEs 
testing <- merge(TNE_eQTLs, TNE_hostgene, by=c("TNE", "strand"), all.y = T)
# 83331 TNEs missing 
# 299266 + 83331 = 382597

print("writing table")
fwrite(TNEs, "TNEs.test.xls", sep="\t", quote = F, col.names = T, row.names = F)


# groupby the eQTL_pos, tissue, rsid -> multiple genes into a list 
# i believe the only NAs should only be from the the caviar and dapg column 
# all_eQTL[is.na(all_eQTL)] <- "NA"
# 
# print("writing table")
testing <- fread("TNEs.test.xls")

minus_PCHiC <- fread("/data/bioinformatics/projects/donglab/AMPPD_eRNA/output/minus/eRNA.PCHiC/eRNA.minus.f22.PCHiCPromoters.bait.score.txt")
colnames(minus_PCHiC) <- c("TNE", "strand", "tissue", "cpp", "gene")
minus_PCHiC$source <- "PCHiC"
minus_PCHiC$ensgs <- ""

plus_PCHiC <- fread("/data/bioinformatics/projects/donglab/AMPPD_eRNA/output/plus/eRNA.PCHiC/eRNA.plus.f22.PCHiCPromoters.bait.score.txt")
colnames(plus_PCHiC) <- c("TNE", "strand", "tissue", "cpp", "gene")
plus_PCHiC$source <- "PCHiC"
plus_PCHiC$ensgs <- ""

PCHiC <- rbind(minus_PCHiC, plus_PCHiC)

# keep all the PCHiC rows (bait and oe interactions)
PCHiC_hostgene <- merge(TNE_hostgene, PCHiC, by=c("TNE", "strand"), all.y = T)

fwrite(PCHiC_hostgene, "TNE.PCHiC.xls", sep="\t", quote = F, col.names = T, row.names = F)

final <- rbind(PCHiC_hostgene, testing)
fwrite(final, "TNE.target.gene.xls", sep="\t", quote = F, col.names = T, row.names = F)







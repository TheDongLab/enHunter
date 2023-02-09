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
caviar <- caviar[, c(1,7,8,14,16,17,19,20)] 
colnames(caviar) <- c("chr", "start", "end", "rsid", "caviar_gene", "caviar_ensg", "tissue", "cpp")
caviar$eQTL_pos <- paste(caviar$chr, caviar$start, caviar$end, sep="_")
caviar <- caviar[,-c(1:3)]

# does NOT contain duplicates -> maybe this is because I added cpp score?
any(duplicated(caviar)) # FALSE

# it IS because you added cpp score 
# gene/tissue/rsid pairs can appear duplicated because the gene version is actually different 

#"eRNA.f18.eSNP.gtexCaviar.txt": contains all the SNPs located within TNEs 
TNE_Caviar <- read.table("/data/bioinformatics/projects/donglab/AMPPD_eRNA/output/snp_table/eRNA.f18.eSNP.gtexCaviar.txt")
colnames(TNE_Caviar) <- c("TNE", "strand", "rsid", "eQTL_pos", "tissue", "caviar")

print("merging Caviar...")

# merge caviar table (with gene and cpp information) with TNE caviar SNP information
TNE_Caviar <- merge(TNE_Caviar , caviar, by=c("eQTL_pos", "tissue", "rsid"), all.x = TRUE) 

# can not group by rsid/tissue/TNE here because then it messes up the cpp/pip score 
# cpp/pip score is different for each rsid/gene pair  

#  group_by(eQTL_pos, tissue, rsid, TNE, strand) %>% 
#  mutate(caviar_genes = paste0(caviar_gene, collapse=",")) %>% 
#  mutate(caviar_ensgs = paste0(caviar_ensg, collapse=",")) %>% 
#  distinct(across(contains("caviar_genes")), .keep_all = TRUE)

# tester <- TNE_Caviar %>% group_by(TNE, strand, rsid, eQTL_pos, tissue, caviar_gene, caviar_ensg) %>% mutate(importance = max(cpp))

# the "duplicate rows" only occur if you groupby gene instead of ensgs 
# there are multiple ensgs for each gene name 
#filter(tester, cpp != importance) %>% nrow() # 0 

tester <- TNE_Caviar %>% group_by(TNE, strand, rsid, eQTL_pos, tissue, caviar_gene, caviar_ensg) %>% mutate(importance = max(cpp))

TNE_Caviar$source <- "caviar"
TNE_Caviar <- TNE_Caviar %>% select(-c(eQTL_pos, rsid, caviar)) %>% rename(gene = caviar_gene, ensgs = caviar_ensg)

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

# note: all.x = TRUE and all.x = FALSE does not change the number of rows 
TNE_Dapg <- merge(TNE_Dapg, dapg, by=c("eQTL_pos", "tissue", "rsid"), all.x = TRUE)

TNE_Dapg <- TNE_Dapg %>% select(-c(eQTL_pos, rsid, dapg)) %>% 
  rename(gene = dapg_gene, ensgs = dapg_ensg, cpp = pip)
TNE_Dapg$source <- "dapg"

### merge with hostgene information table 
TNE_eQTLs <- rbind(TNE_Dapg, TNE_Caviar)

testing <- merge(TNE_eQTLs, TNE_hostgene, by=c("TNE", "strand"), all.y = T)
# 83331 TNEs missing -> do not have target gene information from eQTL snps 
# 299266 + 83331 = 382597

# group by TNE/gene/tissue
# if multiple snps associated with the same gene in one TNE/tissue, then add the cpp scores together 
TNEs <- testing %>% group_by(TNE, strand, ensgs, tissue, source) %>% summarise(importance = sum(cpp)) %>% distinct()

print("writing table")
fwrite(TNEs, "TNEs.test.xls", sep="\t", quote = F, col.names = T, row.names = F)

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

# if gene = NA (baitName) then drop the row 
PCHiC[!is.na(PCHiC$gene), ]

# keep all the PCHiC rows (bait and oe interactions)
PCHiC_hostgene <- merge(TNE_hostgene, PCHiC, by=c("TNE", "strand"), all.y = T)

fwrite(PCHiC_hostgene, "TNE.PCHiC.xls", sep="\t", quote = F, col.names = T, row.names = F)

final <- rbind(PCHiC_hostgene, testing)
fwrite(final, "TNE.target.gene.xls", sep="\t", quote = F, col.names = T, row.names = F)







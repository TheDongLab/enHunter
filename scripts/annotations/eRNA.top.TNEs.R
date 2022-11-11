# finds the top ~10 TNE candidates 
library(data.table)
library(dplyr)

eRNAs <- fread("input_files/characterization/eRNA.characterize.xls", header=T)

# columns used to calculate the features score 
# select=c(f06.TFBS, f07.P300, f08.CAGEbloodenhancer, f09.chromHMM_blood, f12.DNaseROADMAP, f15.HCNE) 
# where TFBS > 5

# some additional features to take into account 
# "f13.phyloP", 
# "f16.GWAS",
# "f18.eSNP" 
# "f10.VISTA", (not the best indicator but maybe??)

# make values binary 
eRNA_scores <- eRNAs %>% filter(features >= max(eRNAs$features)-1 & f13.phyloP > 0 & f16.GWAS > 0 & f18.eSNP > 0) 
eRNA_scores <- eRNA_scores[order(f10.VISTA), ]


write.table(eRNA_scores, file="top_TNES.xls", quote=F, sep="\t", row.names = F)

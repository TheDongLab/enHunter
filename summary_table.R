library(data.table)
library(dplyr)

# create a summary table for the number of TNEs, PCHiC regions, eQTLs for each class 
characterization <- fread("./input_files/characterization/eRNA.characterize.feature.color.xls")

test <- characterization %>% group_by(class) %>% summarise(num_TNEs = n())
caviar <- characterization %>% group_by(class) %>% summarise(num_Caviar = sum(f18.eSNP.gtexCaviar))
dapg <- characterization %>% group_by(class) %>% summarise(num_Dapg = sum(f18.eSNP.gtexDapg))
PCHiC <- characterization %>% group_by(class) %>% summarise(num_PCHiC = sum(f22.PCHiC))



final <- Reduce(merge, list(test, caviar, dapg, PCHiC ))

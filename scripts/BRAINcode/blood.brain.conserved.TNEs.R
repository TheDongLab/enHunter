library(data.table)
library(dplyr)
library(ggplot2)

dopamine <- fread("./input_files/BRAINcode/TNE.dopamine.and.blood.txt") %>% select(V4, V8, V10)
colnames(dopamine) <- c("dopamine", "blood", "strand")

pyramidal <- fread("./input_files/BRAINcode/TNE.pyramdial.and.blood.txt") %>% select(V4, V8, V10)
colnames(pyramidal) <- c("pyramidal", "blood", "strand")

nonneuronal <- fread("./input_files/BRAINcode/TNE.nonneuronal.and.blood.txt") %>% select(V4, V8, V10)
colnames(nonneuronal) <- c("nonneuronal", "blood", "strand")

# total number of TNEs were taken from the files for each cell type found at https://www.humanbraincode.org/download/
# number of TNEs mapped over to hg38 (may differ slightly from number of TNEs found in hg19)
total.dop <- 70996
total.pyr <- 37001 
total.nn <- 19686

# percent of cell type TNEs overlapping blood TNE 
dopamine.percent.overlap <- (dopamine$dopamine %>% unique() %>% length()) / total.dop 
pyramidal.percent.overlap <- (pyramidal$pyramidal %>% unique() %>% length()) / total.pyr
nonneuronal.percent.overlap <- (nonneuronal$nonneuronal%>% unique() %>% length()) / total.nn

cells <- c("dopamine", "pyramidal", "nonneuronal")
percent <- c(dopamine.percent.overlap, pyramidal.percent.overlap, nonneuronal.percent.overlap)

df <- data.frame(cells, percent)

percent <- ggplot(df, aes(x = percent, y = reorder(cells, percent, decreasing = T), fill = cells)) + 
  geom_bar(stat = "identity") + xlim(0, 1) + 
  geom_text(aes(label = percent %>% round(digits = 2))) + theme_bw()

ggsave("braincode.celltype.overlap.pdf", percent, device = "pdf")

#==============================================================#
#### identifying interesting blood + braincode TNE examples ####
#==============================================================#

dopamine <- fread("./input_files/BRAINcode/TNE.dopamine.and.blood.txt")
colnames(dopamine) <- c("dopamine.chr", "dopamine.start", "dopamine.end", "dopamine.name" ,"blood.chr", 
                        "blood.start", "blood.end", "blood.name", "blood.score", "blood.strand")

pyramidal <- fread("./input_files/BRAINcode/TNE.pyramdial.and.blood.txt") 
colnames(pyramidal) <- c("pyramidal.chr", "pyramidal.start", "pyramidal.end", "pyramidal.name" ,"blood.chr", 
                        "blood.start", "blood.end", "blood.name", "blood.score", "blood.strand")

nonneuronal <- fread("./input_files/BRAINcode/TNE.nonneuronal.and.blood.txt")
colnames(nonneuronal) <- c("nonneuronal.chr", "nonneuronal.start", "nonneuronal.end", "nonneuronal.name" ,"blood.chr", 
                        "blood.start", "blood.end", "blood.name", "blood.score", "blood.strand")

eRNA.blood.characterize <- fread("./input_files/characterization/eRNA.characterize.feature.color.xls") %>% select("V1", "class", "features", 
                                                                                                                  "f08.CAGEbloodenhancer", "strand")

dopamine.features <- merge( dopamine, eRNA.blood.characterize, by.y = c("V1", "strand"), by.x = c("blood.name", "blood.strand"))

pyramidal.features <- merge(pyramidal, eRNA.blood.characterize, by.y = c("V1", "strand"), by.x = c("blood.name", "blood.strand"))

nonneuronal.features <- merge(nonneuronal, eRNA.blood.characterize, by.y = c("V1", "strand"), by.x = c("blood.name", "blood.strand"))



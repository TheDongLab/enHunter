#######################################
# Usage: ./bidirectional_classes.R eRNA.characterize.feature.color.xls bidirectional_pairs.txt
# add class information to bidirectional pairs
#######################################
library(data.table)
library(dplyr)

args<-commandArgs(trailingOnly=TRUE)

class <- args[1]
pairs <- args[2]

class <- fread(class) %>% select(V1, strand, class)
pairs <- fread(pairs)
colnames(pairs) <- c("TNE1", "strand1", "TNE2", "strand2", "distance")

merge1 <- merge( pairs, class, by.y = c("V1", "strand"), by.x = c("TNE1", "strand1"))
colnames(merge1) <- c("TNE1", "strand1", "TNE2", "strand2", "distance", "class1")

merge2 <- merge(merge1, class, by.x = c("TNE2", "strand2"), by.y = c("V1", "strand") ) %>% rename(class2 = class)

fwrite(merge2, "bidirectional_pairs_class.csv")

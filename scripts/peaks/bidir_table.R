# script to calculate the directionality score for a given pair of TNEs 
# outputs in table format 
library(data.table)
library(dplyr)

#TODO potentially change this 
setwd("~/Documents/college/dong_lab/code/playground/peaks/attempt-2")

## read in the data

#select only the columns you need 
#name, RPM, pair, RPM, dis
system("cut -f4,7,14,17,21 plus.minus.bed > cut.plus.minus.bed")
system("cut -f4,7,14,17,21 minus.plus.bed > cut.minus.plus.bed")

minus.plus <- fread("cut.minus.plus.bed", sep = "\t")
plus.minus <- fread("cut.plus.minus.bed", sep = "\t")

## remove duplicate V1 chromosome ids and select for the neg distance
remove_dups <- function(df) {
  dups <- df$V1[duplicated(df$V1)]
  
  selector <- sapply(1:nrow(df), function(x) {
    row <- df[x, ]
    if(row$V1 %in% dups) {
      row$V5 < 0
    } else {
      TRUE
    }
  })
  
  final <- df[selector]
  
  #testing
  #anyDuplicated(final$V1)
  #filter(final, V1 %in% dups)
  final
}

minus.plus <- remove_dups(minus.plus)
plus.minus <-remove_dups(plus.minus)

## calculate the directionality score 
# see https://anderssonlab.org/tag/erna/ 

directionality <- function(f, r) {
  d <- (f - r) / (f + r)
}

final.minus.plus <- minus.plus %>% mutate(d= directionality(V4, V2))
final.plus.minus <- plus.minus %>% mutate(d= directionality(V2, V4))

write.table(final.minus.plus %>% select(V1, V3, V5, d), file="final.minus.plus", quote=FALSE, 
            row.names=FALSE, col.names=FALSE, sep="\t")

write.table(final.plus.minus %>% select(V1, V3, V5, d), file="final.plus.minus", quote=FALSE, 
            row.names=FALSE, col.names=FALSE, sep="\t")

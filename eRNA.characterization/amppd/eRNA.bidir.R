# ================================
### calculate the bidriectional transcription distance and 
## Usage: Rscript eRNA.bidir.R cut.plus.minus.bed "plus.minus"
# ================================
library(data.table)
library(dplyr)

args<-commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)<2) {
  stop("Two arguments must be supplied")
} 

## read in the data
bed <- fread(args[1], sep = "\t")

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

bed <- remove_dups(bed)

## calculate the directionality score 
# see https://anderssonlab.org/tag/erna/ 

directionality <- function(f, r) {
  d <- (f - r) / (f + r)
}

if (args[2] == "plus.minus") {
    final.plus.minus <- plus.minus %>% mutate(d= directionality(V2, V4))
    write.table(final.plus.minus %>% select(V1, V3, V5, d), file="eRNA.f29.bidir.txt", quote=FALSE, 
            row.names=FALSE, col.names=FALSE, sep="\t")

} else if (args[2] == "minus.plus") {
    final.minus.plus <- minus.plus %>% mutate(d= directionality(V4, V2))
    write.table(final.minus.plus %>% select(V1, V3, V5, d), file="eRNA.f29.bidir.txt", quote=FALSE, 
            row.names=FALSE, col.names=FALSE, sep="\t")
}

library(data.table)
library(dplyr)

setwd("~/Documents/college/dong_lab/code/playground/peaks/attempt-3")

###### function to generate bed file for UCSC genome visualization 

ucsc_narrowPeak <- function(df) {
  bed.df <- data.table(chrom=character(), start=numeric(), end=numeric(), 
                       name=character(), score=numeric(), strand=character(), 
                       thickStart=numeric(), thickEnd=numeric())
  
  #to iterate through all rows do 1:nrow(df)
  for(i in 1:nrow(df)) {
    row <- df[i, ]
    first_vals <- unlist(strsplit(row$V4, "_"))
    
    first_bed <- data.table(chrom=row$V1, start=first_vals[2], end= first_vals[3], 
                            name=row$V4, score=row$V5, strand=row$V6, 
                            thickStart=row$V2, thickEnd=row$V3)
    
    bed.df <- rbindlist(list(bed.df, first_bed))
  }
  
  bed.df
}

ucsc_closest_bed <- function(df) {
  bed.filtered.df <- data.table(chrom=character(), 
                                start=numeric(), 
                                end=numeric(), 
                                name=character(), 
                                score=numeric(), 
                                strand=character(), 
                                thickStart=numeric(), 
                                thickEnd=numeric())
  
  #to iterate through all rows do 1:nrow(df)
  for(i in 1:nrow(df)) {
    row <- df[i, ]
    first_vals <- unlist(strsplit(row$V4, "_"))

    first_bed <- data.table(chrom=row$V1, start=first_vals[2], end= first_vals[3], 
                            name=row$V4, score=row$V5, strand=row$V6, 
                            thickStart=row$V2, thickEnd=row$V3)
    
    sec_vals <- unlist(strsplit(row$V14, "_"))
    sec_bed <- data.table(chrom=row$V11, start=sec_vals[2], end=sec_vals[3], 
                          name=row$V14, score=row$V15, strand=row$V16, 
                          thickStart=row$V12, thickEnd=row$V13)
    
    bed.filtered.df <- rbindlist(list(bed.filtered.df, first_bed))
    bed.filtered.df <- rbindlist(list(bed.filtered.df, sec_bed))
  }
  
  bed.filtered.df
}

###### importing narrowPeak files for plus and minus strands 
plus.narrow <- fread("./inputs/peaks.n200.plus.narrowPeak", sep = "\t")
minus.narrow <- fread("./inputs/peaks.n200.minus.narrowPeak", sep = "\t")

plus <- ucsc_narrowPeak(plus.narrow)
minus <- ucsc_narrowPeak(minus.narrow)

plus.and.minus <- rbind(plus, minus)
write.table(plus.and.minus, file="ucsc.plus.and.minus.bed", quote=FALSE, 
            row.names=FALSE, col.names=FALSE, sep="\t")

###### with the plus minus orientation and POSITIVE distance only 
plus.minus <- fread("./inputs/plus.minus.bed", sep = "\t")

pos.plus.minus <- plus.minus[plus.minus$V21 > 0]
print(paste0("percentage of TNEs where + is before - ", nrow(pos.plus.minus) / nrow(plus.minus) ))
# ~69%

ucsc.plus.minus <- ucsc_closest_bed(pos.plus.minus)
write.table(ucsc.plus.minus, file="ucsc.plus.minus.pos.dist.bed", quote=FALSE, 
            row.names=FALSE, col.names=FALSE, sep="\t")


###### with the minus plus orientation and POSITIVE distance only 
minus.plus <- fread("./inputs/minus.plus.bed", sep = "\t")

pos.minus.plus <- minus.plus[minus.plus$V21 > 0]
print(paste0("percentage of TNEs where + is before - ", nrow(pos.minus.plus) / nrow(minus.plus) ))
# ~69%

ucsc.minus.plus <- ucsc_closest_bed(pos.minus.plus)
write.table(ucsc.minus.plus, file="ucsc.minus.plus.pos.dist.bed", quote=FALSE, 
            row.names=FALSE, col.names=FALSE, sep="\t")



###### with the plus minus orientation 
plus.minus <- fread("./inputs/plus.minus.bed", sep = "\t")

ucsc.plus.minus <- ucsc_closest_bed(plus.minus)
write.table(ucsc.plus.minus, file="ucsc.plus.minus.bed", quote=FALSE, 
            row.names=FALSE, col.names=FALSE, sep="\t")


###### with the minus plus orientation
minus.plus <- fread("./inputs/minus.plus.bed", sep = "\t")

ucsc.minus.plus <- ucsc_closest_bed(minus.plus)
write.table(ucsc.minus.plus, file="ucsc.minus.plus.bed", quote=FALSE, 
            row.names=FALSE, col.names=FALSE, sep="\t")


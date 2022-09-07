###### function to generate bed file for UCSC genome visualization 
# from the narrowPeak file format given by input bidirectional bedtools closest analysis for TNEs
# aka a weird narrowPeak format that I wrote 

plus.narrow <- fread("./inputs/peaks.n200.plus.narrowPeak", sep = "\t")
minus.narrow <- fread("./inputs/peaks.n200.minus.narrowPeak", sep = "\t")

ucsc_narrowPeak <- function(df) {
  bed.df <- data.table(chrom=character(), start=numeric(), end=numeric(), 
                       name=character(), score=numeric(), strand=character(), 
                       thickStart=numeric(), thickEnd=numeric())
  
  #to iterate through all rows do 1:nrow(df)
  for(i in 1:100) {
    row <- df[i, ]
    first_vals <- unlist(strsplit(row$V4, "_"))
    
    first_bed <- data.table(chrom=row$V1, start=first_vals[2], end= first_vals[3], 
                            name=row$V4, score=row$V5, strand=row$V6, 
                            thickStart=row$V2, thickEnd=row$V3)
    
    bed.df <- rbindlist(list(bed.df, first_bed))
  }
  
  bed.df
}

plus <- ucsc_narrowPeak(plus.narrow)
minus <- ucsc_narrowPeak(minus.narrow)

plus.and.minus <- rbind(plus, minus)
write.table(plus.and.minus, file="ucsc.plus.and.minus.bed", quote=FALSE, 
            row.names=FALSE, col.names=FALSE, sep="\t")

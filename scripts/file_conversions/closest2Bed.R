# convert the output from bedtools closest to bed file 
# peaks are represented in thickStart and thickEnd 
# usage: Rscript closest2Bed.R plus.minus.bed

args <- commandArgs()

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

plus.minus <- fread(args[1], sep = "\t")

ucsc.plus.minus <- ucsc_closest_bed(plus.minus)
write.table(ucsc.plus.minus, file=paste0("ucsc.", args[1]), quote=FALSE, 
            row.names=FALSE, col.names=FALSE, sep="\t")



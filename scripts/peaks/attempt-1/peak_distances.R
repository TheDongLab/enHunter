#!/usr/bin/env Rscript

#NOTE: 
#output distances from this script are negative if - strand is before + strand
#(regardless of the order strand files are inputed in)

suppressPackageStartupMessages(library(argparser))
suppressPackageStartupMessages(library(data.table))

#creating a parser 
p <- arg_parser("Calculate the genomic distance between expression peaks on plus and minus strands")

#add command line arguments 
p <- add_argument(p, "o", help = "order of the strands in bed file; can be either \"plus-minus\" or \"minus-plus\"")
p <- add_argument(p, "bed", help = "output of the bedtools closest function")
p <- add_argument(p, "minus", help = "narrowPeak file for minus strand")
p <- add_argument(p, "plus", help = "narrowPeak file for the plus strand")

args <- parse_args(p)

if ( !(args$o %in% c("plus-minus", "minus-plus")) ) {
  stop("--orientation provided is invalid")
}


#expect distances to be negative 
plus.minus <- fread(args$bed, sep = "\t")
  #read.table("plus.minus.bed", sep = "\t")

#expect distances to be positive 
#minus.plus <- read.table("minus.plus.bed", sep = "\t")

minus.peaks <- fread(args$minus, sep = "\t")

  #read.table("smaller.peaks.n200.minus.narrowPeak")
plus.peaks <- fread(args$plus, sep = "\t")
  #read.table("smaller.peaks.n200.plus.narrowPeak")


final_df <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(final_df) <- c("", "disPeak", "closestTNE")

for(i in 1:nrow(plus.minus)) {

  cur_row <- plus.minus[i, ]
  
  closest_tne <- cur_row[["V8"]]
  
  close_peaks <- if(args$o == "plus-minus") {
                      minus.peaks
                  } else {
                      plus.peaks
                  }
  
  close_peak_row <- close_peaks[V4 == closest_tne, ]
  close_peak <- close_peak_row[["V2"]] + close_peak_row[["V10"]]
  

  cur_tne <- cur_row[["V4"]]
  
  cur_peaks <- if(args$o == "plus-minus") {
    plus.peaks
  } else {
    minus.peaks
  }
  
  cur_peak_row <- cur_peaks[V4 == cur_tne, ]
  cur_peak <- cur_peak_row[["V2"]] + cur_peak_row[["V10"]]
  
  if (args$o == "plus-minus") {
    dist <- close_peak - cur_peak
  } else {
    dist <- cur_peak - close_peak
  }

  temp_df <- data.frame(cur_tne, dist, closest_tne)
  colnames(temp_df) <- c("", "disPeak", "closestTNE")
  
  final_df <- rbind(final_df, temp_df)
  
  #print(paste(cur_tne, closest_tne, ":" , dist))
}

write.table( final_df, file = paste0(args$o, ".bed"), sep = "\t", row.names = FALSE, col.names= FALSE)


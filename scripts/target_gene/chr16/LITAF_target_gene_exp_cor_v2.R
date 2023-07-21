# =============================
## calculates the pearson correlation between input eRNAs and genes 
## produces a ggplot graph 
## v2: revised to make sure values for each sample are plotted together 
## usage: Rscript target_gene_exp_cor.R gene gene_name eRNA eRNA_name
# =============================

library(data.table)
library(dplyr)
library(ggplot2)
library(GGally)
library(stringr)

# reads in the data table and converts it to the right format 
make_table <- function(file) {
  table <- fread(file = file) %>% melt(id = 1, variable.name = "sample" )
  name <- table[1,1][[1]]
  table <- table[, -c(1)]
  colnames(table) <- c("sample", name)
  
  return(table)
}

# eRNA 1 minus
chr16_11611980_11612400 <- make_table("./input_files/eRNAs/tables/chr16_11611980_11612400_exp_table.txt")

# eRNA 1 plus
### gets a warning because 0 is an integer and not double
chr16_11612780_11613560 <- make_table("./input_files/eRNAs/tables/chr16_11612780_11613560_exp_table.txt")

# eRNA 2 minus
chr16_11613470_11613780 <- make_table("./input_files/eRNAs/tables/chr16_11613470_11613780_exp_table.txt") 

# eRNA 2 plus
chr16_11613950_11614560 <- make_table("./input_files/eRNAs/tables/chr16_11613950_11614560_exp_table.txt") 

# eRNA 3 minus 
chr16_11639850_11640300 <- make_table("./input_files/eRNAs/tables/chr16_11639850_11640300_exp_table.txt") 

# eRNA 3 plus 
chr16_11640470_11641120 <- make_table("./input_files/eRNAs/tables/chr16_11640470_11641120_exp_table.txt") 

# gene 
LITAF <-  make_table("./input_files/eRNAs/tables/LITAF_gene_exp_table.txt") 

df <- Reduce(merge,list(chr16_11611980_11612400,chr16_11612780_11613560,chr16_11613470_11613780,
                  chr16_11613950_11614560,chr16_11639850_11640300,chr16_11640470_11641120, 
                  LITAF))

df <- df[-c(1)]
df <- df[is.finite(rowSums(df)),]

eRNA_gene <- ggpairs(df)
eRNA_gene

### make pairs ###
make_pair <- function(df1, df2, colname) {
  combined <- merge(df1, df2)
  combined[colname] <- combined[2] + combined[3]
  
  return(combined[-c(2,3)])
}

eRNA1 <- make_pair(chr16_11611980_11612400, chr16_11612780_11613560, "eRNA1")
eRNA2 <- make_pair(chr16_11613470_11613780, chr16_11613950_11614560, "eRNA2")
eRNA3 <- make_pair(chr16_11639850_11640300, chr16_11640470_11641120, "eRNA3")

eRNA_paired <- Reduce(merge, list(eRNA1, eRNA2, eRNA3, LITAF, CLEC16A))

eRNA_paired <- eRNA_paired[-c(1)]
eRNA_paired <- eRNA_paired[is.finite(rowSums(eRNA_paired)),]

eRNA_pairs <- ggpairs(eRNA_paired)
eRNA_pairs

### transcripts #### 
transcripts <- read.table("/Users/rw552/Documents/amp-pd/enHunter/input_files/eRNAs/LITAF_transcripts.tsv")
split_by_transcript <- split(transcripts, f = transcripts$Name)  

for (df in split_by_transcript) {
  
  df[["TPM"]] <- df[["TPM"]] %>% log(base=10) 
  transcript <- df[["Name"]][[1]]
  TPM <- df[-c(1)]
  
  colnames(TPM) <- c("sample", transcript)
  #TPM[["sample"]] <- str_replace_all(TPM[["sample"]], "-", ".")
  
  full_df <- Reduce(merge, list(eRNA1, eRNA2, eRNA3, TPM))
  
  full_df <- full_df[-c(1)]
  full_df <- full_df[is.finite(rowSums(full_df)),]
  
  eRNA_transcript <- ggpairs(full_df)
  ggsave(paste0(transcript, ".pdf"), plot = eRNA_transcript, width = 10, height = 10)
  
}

#ENST00000571627.5



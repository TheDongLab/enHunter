# =============================
## calculates the pearson correlation between input eRNAs and genes 
## produces a ggplot graph 
## usage: Rscript target_gene_exp_cor.R gene gene_name eRNA eRNA_name
# =============================

library(data.table)
library(dplyr)
library(ggplot2)
library(GGally)

args<-commandArgs(trailingOnly=TRUE)

read_vector <- function(file, name) {
  vector <- scan(file = file, sep = "\t", nlines = 1, na.strings = name)[-1] %>% log10()
  vector[is.infinite(vector)] <- NA
  return(vector)
}

gene <- read_vector(args[1], args[2])
eRNAs <- list()

for (x in seq(from = 3, to = length(args), by = 2) ){
  append(eRNAs, read_vector(args[x], args[x+1]) )
}

df <- as.data.frame(do.call(cbind, eRNAs))
df <- cbind(df, gene)

eRNA_gene <- ggpairs(df)

ggsave(paste0(args[2], ".pdf"), plot = eRNA_gene, width = 13, height = 13)


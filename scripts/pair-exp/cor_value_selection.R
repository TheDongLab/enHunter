library(data.table)
library(dplyr)
library(ggplot2)
library(GGally)
library(stringr)

### copied over from "target_gene_exp_cor_v2.R"
# reads in the data table and converts it to the right format 
make_table <- function(file) {
  table <- fread(file = file) %>% melt(id = 1, variable.name = "sample" )
  name <- table[1,1][[1]]
  table <- table[, -c(1)]
  colnames(table) <- c("sample", name)
  #table[, eval(name) := log10(table[[name]])]
  #table[name] <- log10(table[[name]])
  return(table)
}

# eRNA 1 minus
chr16_11611980_11612400 <- make_table("input_files/eRNAs/tables/chr16_11611980_11612400_exp_table.txt")

# eRNA 1 plus
chr16_11612780_11613560 <- make_table("input_files/eRNAs/tables/chr16_11612780_11613560_exp_table.txt")

# eRNA 2 minus
chr16_11613470_11613780 <- make_table("input_files/eRNAs/tables/chr16_11613470_11613780_exp_table.txt") 

# eRNA 2 plus
chr16_11613950_11614560 <- make_table("input_files/eRNAs/tables/chr16_11613950_11614560_exp_table.txt") 

# eRNA 3 minus 
chr16_11639850_11640300 <- make_table("input_files/eRNAs/tables/chr16_11639850_11640300_exp_table.txt") 

# eRNA 3 plus 
chr16_11640470_11641120 <- make_table("input_files/eRNAs/tables/chr16_11640470_11641120_exp_table.txt") 


make_pair <- function(df1, df2, colname) {
  combined <- merge(df1, df2)
  combined[, eval(colname) := combined[[2]] + combined[[3]]]
  #combined[colname] <- combined[[2]] + combined[[3]]
  
  return(combined[, -c(2,3)])
}

eRNA1 <- make_pair(chr16_11611980_11612400, chr16_11612780_11613560, "eRNA1")

eRNA2 <- make_pair(chr16_11613470_11613780, chr16_11613950_11614560, "eRNA2")

eRNA3 <- make_pair(chr16_11639850_11640300, chr16_11640470_11641120, "eRNA3")

LITAF <- make_table("input_files/LITAF/LITAF_gene_exp_table.txt") 

# log transform the gene TPM values 
LITAF[, ENSG00000189067.12 := log10(ENSG00000189067.12)]

df <- Reduce(merge,list(eRNA1, eRNA2, eRNA3, LITAF))

df <- df[, -c(1)]
df <- df[is.finite(rowSums(df)),]

eRNA_gene <- ggpairs(df)
eRNA_gene

#### testing transcripts correlation with eRNA pairs ####
# eRNA1[, eRNA1 := log10(eRNA1)]
# eRNA2[, eRNA2 := log10(eRNA2)]
# eRNA3[, eRNA3 := log10(eRNA3)]
# 
# transcripts <- read.table("input_files/LITAF/LITAF_transcripts.tsv")
# split_by_transcript <- split(transcripts, f = transcripts$Name)  
# 
# for (df in split_by_transcript) {
#   
#   df[["TPM"]] <- df[["TPM"]] %>% log(base=10) 
#   transcript <- df[["Name"]][[1]]
#   TPM <- df[, -c(1)]
#   
#   colnames(TPM) <- c("sample", transcript)
#   #TPM[["sample"]] <- str_replace_all(TPM[["sample"]], "-", ".")
#   
#   full_df <- Reduce(merge, list(eRNA1, eRNA2, eRNA3, TPM))
#   
#   full_df <- full_df[, -c(1)]
#   full_df <- full_df[is.finite(rowSums(full_df)),]
#   
#   eRNA_transcript <- ggpairs(full_df)
#   ggsave(paste0(transcript, ".pdf"), plot = eRNA_transcript, width = 10, height = 10)
#   
# }


#### combining eRNA information for all eRNAs to get a single correlation value for a gene ####
make_eRNA <- function(df1, df2) {
  combined <- merge(df1, df2)
  combined[, eRNA := combined[[2]] + combined[[3]]]
  #combined[colname] <- combined[[2]] + combined[[3]]
  return(combined[, c("sample","eRNA")])
}

# makes the combined log RPM values for eRNA1, 2 and 3 
eRNA <- Reduce(make_eRNA ,list(eRNA1, eRNA2, eRNA3))
eRNA[, eRNA := log10(eRNA)]

# merges the eRNA df with the gene to be tested log(TPM values)
df <- Reduce(merge,list(eRNA, LITAF))

df <- df[, -c(1)]
df <- df[is.finite(rowSums(df)),]
ggpairs(df)

correlation <- cor.test(df[["eRNA"]], df[["ENSG00000189067.12"]],  method = "pearson")
correlation$estimate 
correlation$p.value

### reading in the TPM values of ALL the neighboring genes in one table ###
original <- fread("input_files/LITAF_locus_expr.txt", header = FALSE)

transposed <- t(original)
colnames(transposed) <- replace(transposed[1, ], 1, "sample")
new <- transposed[-1, ] %>% as.data.frame() 
final <- merge(eRNA, new)


neighboring_genes <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(neighboring_genes) <- c('gene', 'estimate', 'p.value')

for (i in 3:ncol(final)) {
  eRNA <- final[["eRNA"]]
  gene <- final[, ..i][[1]] %>% as.numeric() %>% log10()

  # add eRNA and current gene being studied into separate df 
  temp <- data.frame(eRNA, gene)
  # drop rows with inf 
  temp <- temp[is.finite(rowSums(temp)),]
  
  correlation <- cor.test(temp[[1]], temp[[2]], method = "pearson")
  
  neighboring_genes[nrow(neighboring_genes) + 1, ] <- c(colnames(final)[[i]], correlation$estimate, correlation$p.value)
  
}


## how to visualize this information? 
fwrite(neighboring_genes, file = "LITAF_neighboring_gene_signif.csv")
  
## generating ggpairs plots for top genes 

## first, log transform the eRNA pair values 
eRNA1[, eRNA1 := log10(eRNA1)]
eRNA2[, eRNA2 := log10(eRNA2)]
eRNA3[, eRNA3 := log10(eRNA3)]


make_ggpair <- function(gene, df) {
  gene_df <- df[, .(sample, get(gene))]
  colnames(gene_df) <- c("sample", gene)
  gene_df[, eval(gene) := log10(as.numeric(gene_df[[gene]])) ]
  
  gene_df_eRNA <- Reduce(merge,list(eRNA1, eRNA2, eRNA3, gene_df))
  
  gene_df_eRNA <- gene_df_eRNA[, -c(1)]
  gene_df_eRNA <- gene_df_eRNA[is.finite(rowSums(gene_df_eRNA)),]
  ggsave(paste0(gene, "_pairs.pdf"), plot = ggpairs(gene_df_eRNA), width = 8, height = 8)
}

#ENSG00000277369.1 ( a lncRNA)
make_ggpair("ENSG00000277369.1", final)

# lncRNA <- final[, c("sample", "ENSG00000277369.1")]
# lncRNA$ENSG00000277369.1 <- as.numeric(lncRNA$ENSG00000277369.1) %>% log10()
# 
# lncRNA_eRNA <- Reduce(merge,list(eRNA1, eRNA2, eRNA3, lncRNA))
# 
# lncRNA_eRNA <- lncRNA_eRNA[, -c(1)]
# lncRNA_eRNA <- lncRNA_eRNA[is.finite(rowSums(lncRNA_eRNA)),]
# ggsave("ENSG00000277369.1_pairs.pdf", plot = ggpairs(lncRNA_eRNA))

# ENSG00000189067.12 (LITAF)
make_ggpair("ENSG00000189067.12", final)

# LITAF <- final[, c("sample", "ENSG00000189067.12")]
# LITAF$ENSG00000189067.12 <- as.numeric(LITAF$ENSG00000189067.12) %>% log10()
# LITAF_eRNA <- Reduce(merge,list(eRNA1, eRNA2, eRNA3, LITAF))
# 
# LITAF_eRNA <- LITAF_eRNA[, -c(1)]
# LITAF_eRNA <- LITAF_eRNA[is.finite(rowSums(LITAF_eRNA)),]
# ggpairs(LITAF_eRNA)
# ggsave("ENSG00000189067.12_pairs.pdf", plot = ggpairs(LITAF_eRNA))

# ENSG00000171490.12 
make_ggpair("ENSG00000171490.12", final)

# ENSG00000184602.5
make_ggpair("ENSG00000184602.5", final)



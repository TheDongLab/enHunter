library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)

############ FUNCTIONS ############
make_table <- function(file) {
  table <- fread(file = file) %>% melt(id = 1, variable.name = "sample" )
  name <- table[1,1][[1]]
  table <- table[, -c(1)]
  colnames(table) <- c("sample", name)
  return(table)
}

make_pair <- function(df1, df2, colname) {
  combined <- merge(df1, df2)
  combined[, eval(colname) := combined[[2]] + combined[[3]]]

  return(combined[, -c(2,3)])
}
####################################

####################################
########## read in values ##########

# eRNA 1 minus
chr16_11611980_11612400 <- make_table("input_files/eRNAs/tables/chr16_11611980_11612400_exp_table.txt")

# eRNA 1 plus
### gets a warning because 0 is an integer and not double
chr16_11612780_11613560 <- make_table("input_files/eRNAs/tables/chr16_11612780_11613560_exp_table.txt")

# eRNA 2 minus
chr16_11613470_11613780 <- make_table("input_files/eRNAs/tables/chr16_11613470_11613780_exp_table.txt") 

# eRNA 2 plus
chr16_11613950_11614560 <- make_table("input_files/eRNAs/tables/chr16_11613950_11614560_exp_table.txt") 

# eRNA 3 minus 
chr16_11639850_11640300 <- make_table("input_files/eRNAs/tables/chr16_11639850_11640300_exp_table.txt") 

# eRNA 3 plus 
chr16_11640470_11641120 <- make_table("input_files/eRNAs/tables/chr16_11640470_11641120_exp_table.txt") 

# values are not log10 transformed
eRNA1 <- make_pair(chr16_11611980_11612400, chr16_11612780_11613560, "eRNA1")

eRNA2 <- make_pair(chr16_11613470_11613780, chr16_11613950_11614560, "eRNA2")

eRNA3 <- make_pair(chr16_11639850_11640300, chr16_11640470_11641120, "eRNA3")

eRNA1[, eRNA1 := log10(eRNA1)]
eRNA2[, eRNA2 := log10(eRNA2)]
eRNA3[, eRNA3 := log10(eRNA3)]

# read in gene expression values
LITAF <- make_table("./input_files/eRNAs/tables/LITAF_gene_exp_table.txt") 

# log transform the gene TPM values 
LITAF[, ENSG00000189067.12 := log10(ENSG00000189067.12)]

df <- Reduce(merge,list(eRNA1, eRNA2, eRNA3, LITAF))

# drops the sample names
df <- df[, -c(1)]
df <- df[is.finite(rowSums(df)),]

############### check to see that samples are normally distributed ##############
ggplot(eRNA1, aes(x = eRNA1)) + geom_histogram()
ggplot(eRNA2, aes(x = eRNA2)) + geom_histogram()
ggplot(eRNA3, aes(x = eRNA3)) + geom_histogram()

######## graphing ########
## eRNA1 
eRNA1.correlation <- cor.test(df[["eRNA1"]], df[["ENSG00000189067.12"]],  method = "spearman")
eRNA1.cor.label <- paste0("correlation: ", round(eRNA1.correlation$estimate, 2))
eRNA1.pval.label <- paste0("p-value: ", signif(eRNA1.correlation$p.val, 3))

# t value is the distribution

eRNA1_plot <- ggplot(df, aes(x=eRNA1, ENSG00000189067.12)) + geom_point(colour = "blue") + 
  geom_smooth(method="lm", colour = "red") + 
  annotate(x = -1.05, y = 3,  geom = "text", label = eRNA1.cor.label, size = 5) + 
  annotate(x = -1, y = 2.9,  geom = "text", label = eRNA1.pval.label, size = 5) + 
  theme_bw()

ggsave("./scripts/target_gene/chr16/eRNA1.cor.pdf", plot = eRNA1_plot, device = "pdf")

## eRNA2
eRNA2.correlation <- cor.test(df[["eRNA2"]], df[["ENSG00000189067.12"]],  method = "pearson")
eRNA2.cor.label <- paste0("correlation: ", round(eRNA2.correlation$estimate, 2))
eRNA2.pval.label <- paste0("p-value: ", signif(eRNA2.correlation$p.val, 5))

eRNA2_plot <- ggplot(df, aes(x=eRNA2, ENSG00000189067.12)) + geom_point(colour = "blue") + 
  geom_smooth(method="lm", colour = "red") + 
  annotate(x = -1.05, y = 3,  geom = "text", label = eRNA2.cor.label, size = 5) + 
  annotate(x = -1, y = 2.9,  geom = "text", label = eRNA2.pval.label, size = 5) + 
  theme_bw()

ggsave("./scripts/target_gene/chr16/eRNA2.cor.pdf", plot = eRNA2_plot, device = "pdf")

## eRNA3
eRNA3.correlation <- cor.test(df[["eRNA3"]], df[["ENSG00000189067.12"]],  method = "pearson")
eRNA3.cor.label <- paste0("correlation: ", round(eRNA3.correlation$estimate, 2))
eRNA3.pval.label <- paste0("p-value: ", signif(eRNA3.correlation$p.val, 5))

eRNA3_plot <- ggplot(df, aes(x=eRNA3, ENSG00000189067.12)) + geom_point(colour = "blue") + 
  geom_smooth(method="lm", colour = "red") + 
  annotate(x = -1.05, y = 3,  geom = "text", label = eRNA3.cor.label, size = 5) + 
  annotate(x = -1, y = 2.9,  geom = "text", label = eRNA3.pval.label, size = 5) + 
  theme_bw()

ggsave("./scripts/target_gene/chr16/eRNA3.cor.pdf", plot = eRNA3_plot, device = "pdf")

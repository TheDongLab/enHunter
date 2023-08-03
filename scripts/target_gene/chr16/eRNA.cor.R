# =============================
# 07/31/2023 Rosan Wang
# finds the correlation between the eRNA pairs in the LITAF locus
# =============================
library(data.table)
library(dplyr)
library(ggplot2)

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

##################
### make pairs ###
##################
make_pair <- function(df1, df2, colname) {
  combined <- merge(df1, df2)
  combined[, eval(colname) := combined[[2]] + combined[[3]]]
  
  return(combined[, -c(2,3)])
}

eRNA1 <- make_pair(chr16_11611980_11612400, chr16_11612780_11613560, "eRNA1")
eRNA2 <- make_pair(chr16_11613470_11613780, chr16_11613950_11614560, "eRNA2")
eRNA3 <- make_pair(chr16_11639850_11640300, chr16_11640470_11641120, "eRNA3")

eRNA_paired <- Reduce(merge, list(eRNA1, eRNA2, eRNA3))

### run pearson's correlation test on each pair ###
eRNA_pair_correlations <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(eRNA_pair_correlations) <- c("cur_eRNAs", "eRNAs", "r", "pval")

for (i in 2:length(colnames(eRNA_paired)) ) {

  for (x in 2:(length(colnames(eRNA_paired)))) {
    correlation <- cor.test(eRNA_paired[[i]], eRNA_paired[[x]], method = "pearson")
    eRNA_pair_correlations[nrow(eRNA_pair_correlations) + 1, ] <- c(colnames(eRNA_paired)[[i]], colnames(eRNA_paired)[[x]], 
                                                                   correlation$estimate, correlation$p.value)
  }
}

eRNA_pair_correlations$r <- as.numeric(eRNA_pair_correlations$r)

plot <- ggplot(eRNA_pair_correlations, aes(x =cur_eRNAs, y = eRNAs, fill = r)) + 
  geom_tile(color = "white", lwd = 1.5, linetype = 1) + 
  geom_text(aes(label = signif(r, digits = 3)), color = "white", size = 4) + 
  coord_fixed()
  
plot

ggsave("eRNA.correlations.pdf", plot = plot, device = "pdf")


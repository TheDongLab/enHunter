# =============================
# 07/31/2023 Rosan Wang
# INPUTS:
# - txt of expression values for the neighboring genes (in TPM)
# (see neighboringgenes.sh)
# - txt of eRNA expression values 
# OUTPUT: 
# - Pearson's correlation of the eRNA with each of the provided neighboring genes
# USAGE: 
# Rscript neighboring.genes.heatmap.R neighboring.genes.txt "LITAF.neighbors.csv" eRNA1_exp_table.txt ...
# =============================

library(data.table)
library(dplyr)
library(purrr)
library(ggplot2)
library(BiocManager)
library(biomaRt)

args<-commandArgs(trailingOnly=TRUE)

#neighbors <- "./input_files/eRNAs/tables/HLA_E_locus_expr.txt"
neighbors <- args[1]

#output <- "test.csv"
output <- args[2]

#args <- c(neighbors, output, "./input_files/eRNAs/tables/chr6_30497480_30497870_plus_RPM_table.txt")

### reading in the TPM values of ALL the neighboring genes in one table ###
neighbors <- fread(neighbors, header = FALSE)

# change ENSEMBL ids to SYMBOL 
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

neighbors$ensLookup <- gsub("\\.[0-9]*$", "", neighbors$V1)

annotLookup <- getBM(
  mart=mart,
  attributes=c( "ensembl_gene_id", "external_gene_name"),
  filter="ensembl_gene_id",
  values=neighbors$ensLookup,
  uniqueRows=TRUE)

annotLookup$name <- apply(annotLookup, 1, function(x){
  if (x[["external_gene_name"]] == "" | x[["external_gene_name"]] == "Y_RNA") {
    return(x[["ensembl_gene_id"]])
  } else {
    return(x[["external_gene_name"]])
  }
})

all_neighbors <- merge(neighbors, annotLookup, by.x = "ensLookup", by.y = "ensembl_gene_id", all = T)
temp <- subset(all_neighbors, select = -c(V1, ensLookup, external_gene_name)) %>% dplyr::select(name, everything())
temp[is.na(temp)] <- ""

transposed <- t(temp)
colnames(transposed) <- replace(transposed[1, ], 1, "sample")
new <- transposed[-1, ] %>% as.data.frame() 

### read in eRNA as table ###

## eRNA functions! see cor_value_selection.R ##
make_table <- function(file) {
  table <- fread(file = file) %>% melt(id = 1, variable.name = "sample" )
  name <- table[1,1][[1]]
  table <- table[, -c(1)]
  colnames(table) <- c("sample", name)
  #table[, eval(name) := log10(table[[name]])]
  #table[name] <- log10(table[[name]])
  return(table)
}

make_pair <- function(df1, df2, colname) {
  combined <- merge(df1, df2)
  combined[, eval(colname) := combined[[2]] + combined[[3]]]
  #combined[colname] <- combined[[2]] + combined[[3]]
  
  return(combined[, -c(2,3)])
}

read_eRNAs <- lapply(args[3:length(args)], make_table)

if (length(args) >= 4) {
  eRNA <- Reduce(function(x,y) make_pair(x, y, "eRNA"), read_eRNAs)
} else {
  eRNA <- read_eRNAs[[1]]
  colnames(eRNA)[2] <- "eRNA"
}

eRNA[, eRNA := log10(eRNA)]

### merge eRNA and gene tables ###
final <- merge(eRNA, new)


neighboring_genes <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(neighboring_genes) <- c('gene', 'estimate', 'p.value')


for (i in 3:ncol(final)) {
  eRNA.tmp <- final[["eRNA"]]
  gene <- final[, ..i][[1]] %>% as.numeric() %>% log10()
  
  # add eRNA and current gene being studied into separate df 
  temp <- data.frame(eRNA.tmp, gene)
  # drop rows with inf (happens with log10(0))
  temp <- temp[is.finite(rowSums(temp)),]
  
  correlation <- cor.test(temp[[1]], temp[[2]], method = "pearson")
  
  neighboring_genes[nrow(neighboring_genes) + 1, ] <- c(colnames(final)[[i]], 
                                                        correlation$estimate, 
                                                        correlation$p.value)
}

fwrite(neighboring_genes, file = paste0(output, ".csv"))

neighboring_genes$estimate <- as.double(neighboring_genes$estimate)
neighboring_genes$p.value <- as.double(neighboring_genes$p.value)

plot <- ggplot(neighboring_genes, aes(x=reorder(gene, p.value, decreasing = TRUE), 
                              y=estimate, fill = -log10(p.value))) + 
  geom_bar(stat="identity") + coord_flip() + scale_fill_gradient(high = "red", low = "blue") + 
  xlab("gene") + theme_bw()
plot
ggsave(paste0(output, ".pdf"), plot = plot, device = "pdf", width = 4, height = 9, units = "in")
#ggsave("test.pdf", plot = plot, device = "pdf", width = 5, height = 10,  units = "in")

  
library(ggplot2)

setwd("~/Documents/enHunter/scripts/gene-expression")
exp_genes <- read.table("plus.minus.genes.tsv", sep = "\t", header = TRUE)

# there are some genes which are HIGHLY expressed on the minus and plus strand
# TODO more analysis later
ggplot(exp_genes, aes(x=minus, y=plus, color=strand)) + geom_point()

# more zoomed in 
ggplot(exp_genes, aes(x=minus, y=plus, color=strand)) +
  geom_point() + xlim(0,1) + ylim(0,1) 

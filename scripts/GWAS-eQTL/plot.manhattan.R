# R script to generate manhattan plots and qqplots 
# usage: Rscript plot.manhattan.R 
library(qqman)
library(data.table)
library(dplyr)

args<-commandArgs(TRUE)

ciseQTL = args[1]
output = args[2]

df <- fread(ciseQTL) %>% select(-(beta:"t-stat"))
df$CHR <- 1

par(mar=c(1,1,1,1))

png(file = paste0(output, ".png"), units="px", width=2000, height=500)
manhattan(df, chr="CHR", bp="pos", snp="snpid", p="p-value", logp = TRUE,
          suggestiveline = F, genomewideline = F)
dev.off()

png(file = paste0(output, ".qq.png"), units = "px", width = 500, height = 500)
qq(df$`p-value`)
dev.off()
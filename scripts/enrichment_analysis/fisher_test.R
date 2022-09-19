# =============================
## calculate the p-value using a fishers exact test (one tailed)
## to determine if a tissue/trait is significantly enriched in TNEs 
## for ALL data (not GWAS) 
## usage: Rscript fisher_test.R sample control graphName
# =============================

library(data.table)
library(dplyr)
library(ggplot2)

args<-commandArgs(trailingOnly=TRUE)

sample <- fread(args[1])
control <- fread(args[2])

#sample <- fread("./data/minus/eRNA.minus.f08.CAGEenhtissue.counts")
#control <- fread("./data/all_enhancer_tissues_hg38.counts")

total <- merge(sample, control, by="V1", all.x = TRUE, all.y = FALSE)

colnames(total) <- c("Trait", "TNE", "Total")
total <- total[, TNE := as.numeric(TNE)]
total <- total[, Total := as.numeric(Total)]

# TNE = AB
total <- total[, nAB := Total - TNE]

total <- total[, AnB := sum(total[["TNE"]]) - TNE ]
total <- total[, nAnB := sum(total[["nAB"]]) - nAB ]

# total <- total[, pval := apply(total, 1, function(x) {
#   matrix <- matrix(c(as.numeric(x[["TNE"]]),
#                      as.numeric(x[["nAB"]]),
#                      as.numeric(x[["AnB"]]),
#                      as.numeric(x[["nAnB"]])), nrow = 2)
#   
#   fisher <- fisher.test(matrix, alternative='greater')
#   print(paste(x[["Trait"]], fisher[["p.value"]]))
#   fisher[["p.value"]]
# })]

total <- cbind(total, do.call("rbind", apply(total, 1, function(x) {
  matrix <- matrix(c(as.numeric(x[["TNE"]]),
                     as.numeric(x[["nAB"]]),
                     as.numeric(x[["AnB"]]),
                     as.numeric(x[["nAnB"]])), nrow = 2)
  
  fisher <- fisher.test(matrix, alternative='greater')
  list(pval = fisher[["p.value"]], OR = fisher[["estimate"]][["odds ratio"]])
  })
))


## TODO what is the right statistical test to use here? 

# FDR correction - Bonferroni Correction 
# a/N (0.05 / length(unique(total$Trait)))
# n = sample size (number of diseases/tissues)
# a = pvalue cut off 

## this seems too stringent 

N = 0.05 / length(unique(total$Trait))

signif <- subset(total, pval < N & TNE > 3)
total <- total[, pval := as.numeric(pval)]
p <- ggplot(total, aes(x=reorder(Trait, -pval), y=-log10(pval))) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  geom_hline(yintercept=-log10(N), size=.5,linetype = 2) + xlab("Trait")
p

ggsave(args[3], p)

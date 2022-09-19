# =============================
## calculate the p-value using a fishers exact test (one tailed)
## to determine if a tissue/trait is significantly enriched in TNEs for GWAS data 
## usage: Rscript eRNA.GWAS.enrichment.R sample control graphName
# =============================

# sample and control files have to be grouped by disease/trait tallied with number of occurances 
# sort -k x | bedtools groupby ....

library(data.table)
library(dplyr)
library(ggplot2)

args<-commandArgs(trailingOnly=TRUE)

sample <- fread(args[1])
control <- fread(args[2])

# this line is only needed for GWAS because sort doesn't work all the time? 
control <- control %>% group_by(V1) %>% summarise(TNE=sum(V2))

total <- merge(sample, control, by="V1", all.x = TRUE, all.y = FALSE)

colnames(total) <- c("Trait", "TNE", "Total")

# TNE = AB
total <- total[, nAB := Total - TNE]

total <- total[, AnB := sum(total[["TNE"]]) - TNE ]
total <- total[, nAnB := sum(total[["nAB"]]) - nAB ]

total <- cbind(total, do.call("rbind", apply(total, 1, function(x) {
  matrix <- matrix(c(as.numeric(x[["TNE"]]),
                     as.numeric(x[["nAB"]]),
                     as.numeric(x[["AnB"]]),
                     as.numeric(x[["nAnB"]])), nrow = 2)
  
  print(x)
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

### graphing only significant disease traits (p < 0.05)

# took the TNE > 3 criteria from https://github.com/sterding/BRAINcode/blob/63322c7b192f79fe7819035182d64153a1e5c67d/src/eRNA.TFBSencode.enrichment.R
signif <- subset(total, pval < 0.05 & TNE > 3) 
signif <- signif[, pval := as.numeric(pval)]

p <- ggplot(signif, aes(x=reorder(Trait, -pval), y=-log10(pval))) + 
  geom_bar(stat="identity") + 
  coord_flip() + geom_hline(yintercept=-log10(N), size=.5,linetype = 2)

ggsave(args[3], p)
### all disease traits (grouped together) -> should I group together some categories?
#total <- total[, pval := as.numeric(pval)]


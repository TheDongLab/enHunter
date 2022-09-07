library(data.table)
library(dplyr)
library(ggplot2)
setwd("~/Documents/college/dong_lab/code/playground/enrichment_analysis")

sample <- fread("./data/minus/eRNA.minus.f16.GWASDisease.counts.v2")
control <- fread("./data/GWAS_20220810.v1.02.counts.v2")

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

p <- ggplot(signif, aes(x=reorder(Trait, -pval), y=-log10(pval))) + 
  geom_bar(stat="identity") + 
  coord_flip() + geom_hline(yintercept=-log10(N), size=.5,linetype = 2)
p

### all disease traits (grouped together) -> should I group together some categories?
total <- total[, pval := as.numeric(pval)]


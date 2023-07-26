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

#sample <- fread("./input_files/characterization/feature.enrichment/counts/merged/class1.eRNA.f16.GWASDisease.counts")
#control <- fread("./input_files/characterization/feature.enrichment/counts/GWAS_20220810.v1.02.counts.v2")

total <- merge(sample, control, by="V1", all.x = TRUE, all.y = FALSE)

colnames(total) <- c("Trait", "TNE", "Total")
total <- total[, TNE := as.numeric(TNE)]
total <- total[, Total := as.numeric(Total)]

# TNE = AB
# a single SNP may be in 
# total <- total[, nAB := apply(total, 1, function(x)if (x[["TNE"]] > x[["Total"]]) {0} else {as.numeric(x[["Total"]]) - as.numeric(x[["TNE"]])})]

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

signif <- subset(total, pval < 0.05 & TNE > 3)

if (nrow(signif) > 5) {
  final <- signif[, pval := as.numeric(pval)]
  print('hi')
  #final <- total[, pval := as.numeric(pval)]
  
  #final <- final[1452:1472,]
} else {
  final <- total[, pval := as.numeric(pval)]
  final <- final[order(pval)]
  
  if(nrow(final) > 20) {
  
    final <- final[1:20,]
  }
}

final <- final[, OR := apply(final, 1, function(x) {if (is.infinite(x[["OR"]])) {NA} else {x[["OR"]]} })]
#final <- final[order(pval)]
#final <- final[1:10,]

p <- ggplot(final, aes(x=reorder(Trait, -pval), y=-log10(pval), size=TNE, colour=OR ) ) + 
  geom_point(stat="identity") + scale_colour_gradient(na.value = "red") +
  coord_flip()  + labs(title = args[4]) +
 # labs(title = "class 1 GWAS") +
  theme(text = element_text(size = 15))  + theme_bw() +
  geom_hline(yintercept=-log10(N), size=.5,linetype = 2) + xlab("Trait")
p

# resize gtex graphs (more narrow )
ggsave(args[3], plot = p, device = "pdf", width = 30, height = 10, units = "mm")

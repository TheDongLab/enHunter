# =============================
## calculate the p-value using a fishers exact test (one tailed)
## to determine if DE TNEs is significantly enriched in a gene 
## usage: Rscript DE_gene_enrichment.R.R sample control graphName
# =============================

library(data.table)
library(dplyr)
library(ggplot2)

args<-commandArgs(trailingOnly=TRUE)

#sample <- fread(args[1])
#control <- fread(args[2])

# number of differentially expressed TNEs inside the host gene 
sample <- fread("./input_files/hostGenes/DEeRNAs.f19.Hostgene.counts")

# number of non differentially expressed TNEs inside the host gene 
control <- fread("./input_files/hostGenes/nonDEeRNAs.f19.Hostgene.counts")

total <- merge(sample, control, by="V1", all.x = TRUE, all.y = FALSE)

colnames(total) <- c("Trait", "DE_TNE", "non_DE")
total[is.na(total$non_DE)]$non_DE <- 0

total <- total[, DE_TNE := as.numeric(DE_TNE)]
total <- total[, non_DE := as.numeric(non_DE)]

#AB: DE_TNE (DE TNEs in the gene )
#nAB: (DE TNEs not in gene )
#AnB: non_DE (non differential enhancers in gene)
#nAnB: (non differential enhancers not in gene )

total <- total[, nAB := sum(DE_TNE) - DE_TNE]
total <- total[, nAnB := sum(non_DE) - non_DE ]

total <- cbind(total, do.call("rbind", apply(total, 1, function(x) {
  matrix <- matrix(c(as.numeric(x[["DE_TNE"]]),
                     as.numeric(x[["nAB"]]),
                     as.numeric(x[["non_DE"]]),
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

signif <- subset(total, pval < 0.05 & DE_TNE > 3)

if (nrow(signif) > 5) {
  final <- signif[, pval := as.numeric(pval)]
  
  #final <- total[, pval := as.numeric(pval)]
  #final <- final[order(pval)]
  #final <- final[1452:1472,]
} else {
  final <- total[, pval := as.numeric(pval)]
  final <- final[order(pval)]
  
  if(nrow(final) > 20) {
    final <- final[1:20,]
  }
}

final <- final[, OR := apply(final, 1, function(x) {if (is.infinite(x[["OR"]])) {NA} else {x[["OR"]]} })]

p <- ggplot(final, aes(x=reorder(Trait, -pval), y=-log10(pval), size=DE_TNE, colour=OR ) ) + 
  geom_point(stat="identity") + scale_colour_gradient(na.value = "red") +
  coord_flip()  + 
  #labs(title = "class 1 GWAS") +
  labs(title = args[4]) +
  geom_hline(yintercept=-log10(N), size=.5,linetype = 2) + xlab("Trait")
p

# resize gtex graphs (more narrow )
ggsave(args[3], plot = p, device = "pdf", width = 15, height = 50)



##### ahahahhah redoing it 

library(data.table)
library(dplyr)
library(ggplot2)

df <- data.frame(name= "DE gene", 
                 AB = 103, 
                 AnB = 769, 
                 nAB = 629, 
                 nAnB = 34217)

matrix <- matrix(c(103, 629, 769, 34217), nrow = 2)

fisher <- fisher.test(matrix, alternative='greater')

df$pval = fisher[["p.value"]]
df$OR = fisher[["estimate"]][["odds ratio"]]

p <- ggplot(df, aes(x=name, y=-log10(pval), size=AB, colour=OR ) ) + 
  geom_point(stat="identity") + scale_colour_gradient(na.value = "red") +
  coord_flip()  
 # geom_hline(yintercept=-log10(N), size=.5,linetype = 2) + xlab("Trait")

p

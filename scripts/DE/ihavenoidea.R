

# gene_ids = c("ENSG00000111424.10","ENSG00000003402.19","ENSG00000093010.13")
# for (gene_i in gene_ids){
#   # gene_i = "ENSG00000141084.10"
#   gene = as.numeric(as.vector(tpm[rownames(tpm)==gene_i,]))
#   condition = covarianceTable$CONDITION
#   df = data.frame(gene,condition)
#   df$gene = log((df$gene + 1), 2)
#   pdf(file.path(output_dir, paste0(gene_i ,"_boxplot.pdf")), paper = 'USr')
#   boxplot(gene~condition,
#           main = gene_i,
#           data=df,
#           xlab="Condition",
#           ylab="Expression log2(tpm+1)"
#   )
#   dev.off()
#   print(resLFC[gene_i,"padj"])
# }

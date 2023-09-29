### making a manhattan plot ###

# Load the library
library(qqman)

par(mar=c(1,1,1,1))
png(file = "sample.manhattan.png", units="px", width=2000, height=500)
# Make the Manhattan plot on the gwasResults dataset
manhattan(gwasResults, chr="CHR", bp="BP", snp="SNP", p="P",
          suggestiveline = F, genomewideline = F)
dev.off()

png(file = "sample.qq.png", units = "px", width = 500, height = 500)
qq(gwasResults$P)
dev.off()

## see https://github.com/sterding/BRAINcode/blob/master/modules/_eQTL_manhanttenPlot.R 
#Find the largest (unadjusted) p-value for which the FDR is below the desired level. Draw the line at that value of p.
# -log10(max(df$p.value[df$FDR<=0.05])

# see https://r-graph-gallery.com/101_Manhattan_plot.html 

snps <- fread(paste0(base.dir, "/snps/chr1.snps.coordinates")) %>% select(-chr)

test <- merge(snps, tmp, by.x = "snpid", by.y = "SNP")




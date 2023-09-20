# step 1: read in single nucleotide variants 
library(motifbreakR)
library(BSgenome) 
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38)
library(stringr)

snps.bed.file <- "input_files/eRNA.LITAF.snps.bed"
# see the contents
pca.snps <- read.table(snps.bed.file, header = FALSE)$V4

snps.mb <- snps.from.rsid(rsid = pca.snps,
                          dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh38,
                          search.genome = BSgenome.Hsapiens.UCSC.hg38)

snps.mb

# step 2: find broken motifs 
library(MotifDb)

# check which organisms are available under which sources
table(mcols(MotifDb)$organism, mcols(MotifDb)$dataSource) 

# the MotifList introduced by MotifDb includes additional sets of motifs 
data(motifbreakR_motif)
motifbreakR_motif
# include ENCODE-motif, FactorBook, HOCOMOCO, and HOMER

results <- motifbreakR(snpList = snps.mb[1:12], filterp = TRUE,
                       pwmList = motifbreakR_motif,
                       threshold = 1e-4,
                       method = "default",
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::SerialParam())

x = data.frame(results)
x$motifPos <- as.character(x$motifPos)

write.table(x, file = "tfbs.motifs.snps.tsv", sep="\t", 
            col.names=T, row.names=F, quote=F)

plotMB(results = results, rsid = "rs8052797", effect = "strong")
plotMB(results = results, rsid = "rs13330449", effect = "strong")
plotMB(results = results, rsid = "rs13335254", effect = "strong")
plotMB(results = results, rsid = "rs13331560", effect = "strong")
plotMB(results = results, rsid = "rs12149000", effect = "strong")
plotMB(results = results, rsid = "rs11075004", effect = "strong")
plotMB(results = results, rsid = "rs9932684", effect = "strong")
plotMB(results = results, rsid = "rs9932278", effect = "strong")
plotMB(results = results, rsid = "rs1344533", effect = "strong")
plotMB(results = results, rsid = "rs12444261", effect = "strong") 
plotMB(results = results, rsid = "rs12445804", effect = "strong") 





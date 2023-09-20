# step 1: read in single nucleotide variants 
library(motifbreakR)
library(BSgenome) 
library(SNPlocs.Hsapiens.dbSNP155.GRCh37)
library(BSgenome.Hsapiens.UCSC.hg19)

snps.bed.file <- system.file("extdata", "test.snps.bed", package = "motifbreakR")
# see the contents
read.table(snps.bed.file, header = FALSE)

snps.mb.frombed <- snps.from.file(file = snps.bed.file,
                                  dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh37,
                                  search.genome = BSgenome.Hsapiens.UCSC.hg19,
                                  format = "bed", check.unnamed.for.rsid = TRUE)
snps.mb.frombed

 

### Here we can see which organisms are available under which sources
### in MotifDb
table(mcols(MotifDb)$organism, mcols(MotifDb)$dataSource)

data(hocomoco)
hocomoco

data(motifbreakR_motif)
motifbreakR_motif

results <- motifbreakR(snpList = snps.mb[1:12], filterp = TRUE,
                       pwmList = subset(MotifDb, 
                                        dataSource %in% c("HOCOMOCOv11-core-A", "HOCOMOCOv11-core-B", "HOCOMOCOv11-core-C")),
                       threshold = 1e-4,
                       method = "ic",
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::SerialParam())

rs1006140 <- results[results$SNP_id == "rs1006140"]
rs1006140

rs1006140 <- calculatePvalue(rs1006140, granularity = 1e-6)
rs1006140

plotMB(results = results, rsid = "rs1006140", effect = "strong", altAllele = "G")

##### a seperate sample analyis ##########
# prepare variants
load(system.file("extdata",
                 "pca.enhancer.snps.rda",
                 package = "motifbreakR")) # loads snps.mb
pca.enhancer.snps <- sample(snps.mb, 20)
# Get motifs to interrogate
data(hocomoco)
motifs <- sample(hocomoco, 50)
# run motifbreakR
results <- motifbreakR(pca.enhancer.snps,
                       motifs, threshold = 0.85,
                       method = "ic",
                       BPPARAM=BiocParallel::SerialParam())

plotMB(results = results, rsid = "rs11986220")

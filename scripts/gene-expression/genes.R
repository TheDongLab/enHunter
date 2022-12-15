# Usage: Rscript genes.R strand.minus.tab strand.plus.tab strand.genes.tab out.tsv 
# to get strand.*.tab files: cut -f1,5 *.tab 

args <- commandArgs(TRUE)

minus <- read.table(args[1], sep = "\t")
colnames(minus) <- c("gene", "minus")

plus <- read.table(args[2], sep = "\t")
colnames(plus) <- c("gene", "plus")
 
genes <- read.table(args[3], sep = "\t")
colnames(genes) <- c("gene", "strand")

df <- merge(minus, plus, by="gene")
df <- merge(df, genes, by="gene") 

write.table(df, file=args[4], quote= FALSE, sep = "\t", row.names = FALSE) 
 

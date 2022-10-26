library(ggplot2)

setwd("~/Documents/enHunter/scripts/gene-expression")
exp_genes <- read.table("plus.minus.genes.tsv", sep = "\t", header = TRUE)

# there are some genes which are HIGHLY expressed on the minus and plus strand
# TODO more analysis later
ggplot(exp_genes, aes(x=minus, y=plus, color=strand)) + geom_point()

# more zoomed in 
ggplot(exp_genes, aes(x=minus, y=plus, color=strand)) +
  geom_point() + xlim(0,1) + ylim(0,1) 

exp_transcripts <- read.table("plus.minus.transcripts.tsv", sep = "\t", header = TRUE)
transcripts <- ggplot(exp_transcripts, aes(x=minus, y=plus, color=strand)) + geom_point()
ggsave("transcripts.png", plot = transcripts)

transcripts_zoomed <- ggplot(exp_transcripts, aes(x=minus, y=plus, color=strand)) + geom_point() + xlim(0,1) + ylim(0,1) 
ggsave("zoom_transcripts.png", plot = transcripts_zoomed)

exp_transcripts$ratio <- exp_transcripts$minus/exp_transcripts$plus

# ratio < 1 AND strand = minus 
wrong_minus <- exp_transcripts[exp_transcripts$ratio < 1 & exp_transcripts$strand == "-" , ]
# ratio > 1 AND strand = plus 
wrong_plus <- exp_transcripts[exp_transcripts$ratio > 1 & exp_transcripts$strand == "+" , ]

# ratio = 1 (thats just wrong)
# exp_transcripts[exp_transcripts$ratio == 1, ]

# protein coding only 
protein_coding <- exp_transcripts[grep("protein_coding", exp_transcripts$gene),]
protein <- ggplot(protein_coding, aes(x=minus, y=plus, color=strand)) + geom_point() 
ggsave("transcript_protein.png", plot = protein)

library(ggiraph)
library(plotly)

protein <- ggplot(protein_coding) + geom_point_interactive(aes(x=minus, y=plus, color=strand, tooltip = gene)) + theme_minimal()

x <- girafe(ggobj = protein)
x <- girafe_options(x = x, opts_zoom(min = 1, max = 10))

print(x)

protein_zoomed <- ggplot(protein_coding) + geom_point_interactive(aes(x=minus, y=plus, color=strand, tooltip = gene)) + 
  theme_minimal() + xlim(0,100) + ylim(0,100) 

z <- girafe(ggobj = protein_zoomed)
z <- girafe_options(x = z)

print(z)

protein_coding$ratio <- protein_coding$minus/protein_coding$plus

# ratio < 1 AND strand = minus 
wrong_minus_protein <- protein_coding[protein_coding$ratio < 1 & protein_coding$strand == "-" , ]
# ratio > 1 AND strand = plus 
wrong_plus_protein <- protein_coding[protein_coding$ratio > 1 & protein_coding$strand == "+" , ]

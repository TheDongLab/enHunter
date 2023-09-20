library(ggplot2)

setwd("./input_files/split.strands")

### 8/2/2023: strand splitting analysis based on genes less accurate than transcripts
#exp_genes <- read.table("plus.minus.genes.tsv", sep = "\t", header = TRUE)

# there are some genes which are HIGHLY expressed on the minus and plus strand
# TODO more analysis later
#ggplot(exp_genes, aes(x=minus, y=plus, color=strand)) + geom_point()

# more zoomed in 
#ggplot(exp_genes, aes(x=minus, y=plus, color=strand)) +
  #geom_point() + xlim(0,1) + ylim(0,1) 

exp_transcripts <- read.table("plus.minus.transcripts.tsv", sep = "\t", header = TRUE)
transcripts <- ggplot(exp_transcripts, aes(x=minus, y=plus, color=strand)) + geom_point()
transcripts
ggsave("transcripts.png", plot = transcripts)

transcripts_zoomed <- ggplot(exp_transcripts, aes(x=minus, y=plus, color=strand)) + geom_point() + xlim(0,1) + ylim(0,1) 
ggsave("zoom_transcripts.png", plot = transcripts_zoomed)

exp_transcripts$ratio <- exp_transcripts$minus/exp_transcripts$plus

# ratio < 1 AND strand = minus 
wrong_minus <- exp_transcripts[exp_transcripts$ratio < 1 & exp_transcripts$strand == "-" , ]
# ratio > 1 AND strand = plus 
wrong_plus <- exp_transcripts[exp_transcripts$ratio > 1 & exp_transcripts$strand == "+" , ]

# ratio = 1 (thats just wrong)
# good thing is that there are 0 transcript with ratios of 1 (equal in plus and minus) ! 
exp_transcripts[exp_transcripts$ratio == 1, ]

# protein coding only 
protein_coding <- exp_transcripts[grep("protein_coding", exp_transcripts$gene),]
wrong_plus_protein <- protein_coding[protein_coding$ratio > 1 & protein_coding$strand == "+" , ]
wrong_minus_protein <- protein_coding[protein_coding$ratio < 1 & protein_coding$strand == "-" , ]

write.table(wrong_plus_protein, "plus.mapping.to.minus.protein.tsv", quote = F, sep = "\t", row.names = F)
write.table(wrong_minus_protein, "minus.mapping.to.plus.protein.tsv", quote = F, sep = "\t", row.names = F)


colors <- c("-" = "blue", "+" = "red")
protein <- ggplot(protein_coding, aes(x=minus, y=plus, color=strand)) +
  geom_point() + scale_color_manual(values=colors) + theme_bw()
protein
ggsave("transcript_protein.png", plot = protein, device = "png")


### same graph as protein without HBB and HBA1/2
no.globin <- protein_coding[!startsWith( protein_coding$gene, "HBB"),]
no.globin <- no.globin[!startsWith(no.globin$gene, "HBA"),]

no.globin.plot <- ggplot(no.globin, aes(x=minus, y=plus, color=strand)) +
  geom_point() + scale_color_manual(values=colors) + theme_bw()
no.globin.plot
ggsave("no.globin.transcript_protein.png", plot = no.globin.plot, device = "png")

######################interactive graphs######################################
library(ggiraph)
library(plotly)

protein <- ggplot(protein_coding) + geom_point_interactive(aes(x=minus, y=plus, color=strand, tooltip = gene)) + theme_minimal()

x <- girafe(ggobj = protein)
x <- girafe_options(x = x, opts_zoom(min = 1, max = 10))

print(x)

no.globin.plot <- ggplot(no.globin) + geom_point_interactive(aes(x=minus, y=plus, color=strand, tooltip = gene)) + theme_minimal()
a <- girafe(ggobj = no.globin.plot)
a <- girafe_options(x = a, opts_zoom(min = 1, max = 10))
print(a)

protein_zoomed <- ggplot(protein_coding) + geom_point_interactive(aes(x=minus, y=plus, color=strand, tooltip = gene)) + 
  theme_minimal() + xlim(0,100) + ylim(0,100) 

z <- girafe(ggobj = protein_zoomed)
z <- girafe_options(x = z)

print(z)

## with log10(x+1) scale 
protein_log10 <- ggplot(protein_coding) + geom_point_interactive(aes(x=log10(minus + 1), y=log10(plus + 1), color=strand, tooltip = gene)) + theme_minimal()
w <- girafe(ggobj = protein_log10)
w <- girafe_options(x = w)
print(w)

protein_coding$ratio <- protein_coding$minus/protein_coding$plus

# ratio < 1 AND strand = minus 
wrong_minus_protein <- protein_coding[protein_coding$ratio < 1 & protein_coding$strand == "-" , ]
# ratio > 1 AND strand = plus 
wrong_plus_protein <- protein_coding[protein_coding$ratio > 1 & protein_coding$strand == "+" , ]


html_graph <- function(df) {
  protein_log10 <- ggplot(df) + geom_point_interactive(aes(x=log10(minus + 1), y=log10(plus + 1), color=strand, tooltip = gene)) + theme_minimal()
  w <- girafe(ggobj = protein_log10)
  w <- girafe_options(x = w)
  print(w)
}

######## reordering data ###########
library(dplyr)
# X = minus 
# Y = plus 

# minus >= plus plot 
# plot + strand genes first, - genes second 
# move minus genes to the bottom 
protein_plus_minus <- protein_coding %>% filter(minus >= plus) %>% arrange(strand)
tail(protein_plus_minus)

html_graph(protein_plus_minus)
# X < Y plot 
# 
protein_minus_plus <- protein_coding %>% filter(minus < plus) %>% arrange(desc(strand))
tail(protein_minus_plus)
html_graph(protein_minus_plus)

final <- ggplot() + geom_point_interactive(protein_minus_plus, mapping = aes(x=log10(minus + 1), y=log10(plus + 1), color=strand, tooltip = gene)) + 
  geom_point_interactive(protein_plus_minus, mapping = aes(x=log10(minus + 1), y=log10(plus + 1), color=strand, tooltip = gene)) + theme_minimal()

w <- girafe(ggobj = final)
w <- girafe_options(x = w)
print(w)






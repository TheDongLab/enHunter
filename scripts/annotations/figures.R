### 08/02/2023 making figures for paper 
# 1. TNE length 
# 2. TSS distance 
# 3. GC content 

library(ggplot2)
library(dplyr)
library(data.table)

##### TNE length plot 
lengths <- fread("./input_files/eRNA.stranded.length.bed") 
colnames(lengths) <- c("TNE", "strand", "length")

med <- median(lengths$length, na.rm=TRUE)

length.plot <- ggplot(lengths, aes(x=length)) + 
  geom_histogram(fill="gray", binwidth=0.1)  + 
  scale_x_continuous(trans="log10", n.breaks = 28) + 
  scale_y_continuous(n.breaks= 10) + 
  theme_bw() + 
  geom_vline(xintercept=med, color="red", linetype="dashed") +
  theme(axis.line = element_line(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_blank(), 
        text = element_text(size = 20)) + xlab("TNE Length (bp)") + ylab("Count of TNEs") +
  annotate("text", x = med+560, y=5000, size = 6,
           label=paste0("median: ", round(med, digits = 2), " bp"), color="red")

length.plot
ggsave("tne.length.pdf", plot=length.plot, device="pdf")

#----------------- graphs which require eRNA characterization data -------------------------# 

char <- fread("./input_files/eRNA.characterize.feature.color.xls") 

##### TSS plot 
intergenic <- round((nrow(char[f01.dis2TSS < 0]) / nrow(char) ) * 100, digits = 2)
intergenic_txt <- paste0("intergenic (N=", nrow(char[f01.dis2TSS < 0]), "; ", intergenic, "%)")

intronic <- round((nrow(char[f01.dis2TSS > 0]) / nrow(char) ) * 100, digits = 2)
intronic_txt <- paste0("intronic (N=", nrow(char[f01.dis2TSS > 0]), "; ", intronic, "%)")

dis2TSS <- ggplot(char, aes(x=f01.dis2TSS/1000, fill=f01.dis2TSS < 0)) + 
  geom_histogram( binwidth=20) +  
  theme(axis.line = element_line(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        panel.border = element_blank(), legend.title=element_blank(), legend.position=c(0.3, 0.9), text = element_text(size = 20)) +
  xlab("Distance to the TSS (x 1000 bp)") + 
  scale_x_continuous(n.breaks = 19) + scale_fill_manual(values = c("purple", "green"), labels = c(intronic_txt, intergenic_text))
dis2TSS
ggsave("dis2TSS_minus.pdf",plot=dis2TSS, device = "pdf")

##### CpG content 

# id, norm_CpG, GC_content 

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
ggsave("dis2TSS.pdf",plot=dis2TSS, device = "pdf")

##### CpG content 

# id, norm_CpG, GC_content
rand.CpG.minus <- fread("./input_files/CpG/new.random.minus.f05.CpG.txt")
colnames(rand.CpG.minus) <- c("id", "normalized.score", "GC.content")
rand.CpG.minus$name <- "rand"
#rand.CpG.minus[, c("id", "GC.content"):=NULL] 

rand.CpG.plus <- fread("./input_files/CpG/new.random.plus.f05.Cp.txt")
colnames(rand.CpG.plus) <- c("id", "normalized.score", "GC.content")
rand.CpG.plus$name <- "rand"
#rand.CpG.plus[, c("id", "GC.content"):=NULL] 

#prom.CpG.minus is the same as prom.CpG.plus 
#(both identify promoter regions of all genes)
prom.CpG.minus <- fread("./input_files/CpG/new.promoters.minus.f05.CpG.txt")
colnames(prom.CpG.minus) <- c("id", "normalized.score", "GC.content")
prom.CpG.minus$name <- "prom"
#prom.CpG.minus[, c("id", "GC.content"):=NULL] 

#TNE.CpG <- data.frame(normalized.score = char[["f05.CpG"]], name = "TNE")

eRNA.minus.CpG <- fread("./input_files/CpG/new.eRNA.minus.f06.CpG.txt")
colnames(eRNA.minus.CpG) <- c("id", "normalized.score", "GC.content")
eRNA.minus.CpG$name <- "TNE"

eRNA.plus.CpG <- fread("./input_files/CpG/new.eRNA.plus.f05.CpG.txt")
colnames(eRNA.plus.CpG) <- c("id", "normalized.score", "GC.content")
eRNA.plus.CpG$name <- "TNE"

l <- list( eRNA.plus.CpG,eRNA.plus.CpG , prom.CpG.minus, rand.CpG.plus)
CpG <- rbindlist(l) 

# remove normalization scores of NA, this is likely due to 
# hard masked repeat regions/undefined sequences 
filtered.CpG <- CpG %>% filter(!is.na(normalized.score))

# plot 
colors <- c(TNE = "red", prom = "green", rand = "black")
gc.content <- ggplot(filtered.CpG, aes(GC.content, color = name)) + 
  geom_density(adjust = 1.5) + theme_bw() + 
  scale_color_manual(values = colors)
gc.content 

ggsave("gc.content.pdf",plot=gc.content, device = "pdf")

###############################################################################
# making pie chart
all <- fread("./input_files/eRNA.characterize.feature.color.xls")
class1.or.2 <- fread("./input_files/eRNA.characterize.feature.color.xls") %>% subset(class == 1 | class == 2)

DNaseENCODE.TNEs <- ifelse(all$f12.DNaseENCODE > 0, 1, 0) 

# Open a pdf file
pdf("DNase.ENCODE.pdf") 

## only class 1 or 2 enhancers 
DNaseENCODE <- c(sum(DNaseENCODE.TNEs), length(DNaseENCODE.TNEs) - sum(DNaseENCODE.TNEs))
pie(DNaseENCODE, labels = c("Overlap with DNase ENCODE", ""))

dev.off() 

percent.overlap.TNEs <- (sum(DNaseENCODE.TNEs) / length(DNaseENCODE.TNEs)) * 100
percent.overlap.TNEs #12.08792

###### test TNEs with overlap with enhancers found from Kim et al. (uPRO? uChRO?)
TNE.class <- all %>% select(V1, class)
colnames(TNE.class) <- c("TNE", "class")

cCRE.df <- fread("./input_files/eRNA.cCRE.intersect") 
colnames(cCRE.df) <- c("TNE", "cCRE.count")

cCRE.df.class <- merge(TNE.class, cCRE.df, by="TNE")

class1.or.2.cCRE.df <- cCRE.df.class %>% subset(class == 1 | class == 2)

cCRE.count <- ifelse(class1.or.2.cCRE.df$cCRE.count > 0, 1, 0) 
cCRE <- c(sum(cCRE.count), length(cCRE.count) - sum(cCRE.count))
pie(cCRE)

######## finding the percent of enhancers supported by features ######
all <- fread("./input_files/eRNA.characterize.feature.color.xls")

all.plot <- ggplot(all, aes(as.factor(class), fill = as.factor(class))) + geom_bar() + 
  theme(axis.line = element_line(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        panel.border = element_blank(), text = element_text(size = 20)) + ylab("Count") + 
  geom_text(position = position_dodge(width = 0.9), vjust = 0, stat="count", aes(label=after_stat(count))) +
  xlab("Class") + scale_fill_manual(values=c('#ff3333', '#ff9966', '#ffcccc'))
all.plot

ggsave("all.tnes.class.pdf",plot=all.plot, device = "pdf")

num.tne.annotation <- nrow(all[class== 2 | class== 1])
per.tne.annocation <- (num.tne.annotation / nrow(all)) * 100



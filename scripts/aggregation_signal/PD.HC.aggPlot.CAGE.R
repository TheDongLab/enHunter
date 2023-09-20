library(ggplot2)
library(dplyr)
library(stringr)

dir="/data/bioinformatics/projects/donglab/AMPPD_eRNA/binned/"

make_graph <- function(df, lim) {
  plot <- ggplot(df, aes(x = as.numeric(distance), y = avg, color = strand)) + 
    geom_line() + geom_vline(xintercept = 25.5, linetype="dashed") + 
    theme_bw() + scale_color_manual(values=c("blue", "red")) + 
    ylim(0, lim) + scale_x_continuous(n.breaks = 5)
  
  return(plot)
}

read_df <- function(plus, minus, name) {
  PLUS=read.table(plus, check.names =F, stringsAsFactors=F, header=F)[, c(1,5)]; 
  colnames(PLUS)=c("names", name)
  
  MINUS=read.table(minus, check.names =F, stringsAsFactors=F, header=F)[, c(1,5)]; 
  colnames(MINUS)=c("names", name)
  
  # add TNE and distance to mid point columns
  
  PLUS[c('chr', 'start', 'end', 'distance')] <- str_split_fixed(PLUS$names, '_', 4)
  PLUS['TNE'] <- paste(PLUS$chr, PLUS$start, PLUS$end, sep="_")
  
  MINUS[c('chr', 'start', 'end', 'distance')] <- str_split_fixed(MINUS$names, '_', 4)
  MINUS['TNE'] <- paste(MINUS$chr, MINUS$start, MINUS$end, sep="_")
  #MINUS[name] <- MINUS[name] * -1 
  
  PLUS$strand <- "+"
  MINUS$strand <- "-"
  
  df <- rbind(MINUS, PLUS)
  
  return(df)
}

features = read.table("/data/bioinformatics/projects/donglab/AMPPD_eRNA/eRNA.characterize.feature.color.xls", 
                      check.names =F, stringsAsFactors=F, header=T, row.names = NULL)

colnames(features)[colnames(features) == "row.names"] ="TNE"
features = features[, c("TNE", "class", "features", "strand")]


group_dist <- function(input, col){
  input %>% group_by(distance, strand) %>% summarise(avg = mean(get(col), trim=0.01))
}

############## 1. combined PD + HC CAGE signal from blood ###################
# - combined.ctss
binfile = paste0(dir, "combined.ctss")
df <- read_df(paste0(binfile, ".plus.500bp.tab"), paste0(binfile, ".minus.500bp.tab"), "combined.ctss")

# add class filters 
combined = merge(df, features, by = c("strand", "TNE") )

CAGE <- combined %>% filter(class == 1)
CAGE <- CAGE %>% group_by(distance, strand) %>% summarise(avg = mean(combined.ctss, trim=0.01))

pdf("bidir.CAGE.pdf")
make_graph(CAGE, 0.05)
dev.off()

############## 2. RNAseq signal from blood ###################
RNAseqplusfile = paste0(dir, "combined.mean.normalized.random.samplesN200.plus.500bp.tab")
RNAseqminusfile = paste0(dir, "combined.mean.normalized.random.samplesN200.minus.500bp.tab")

RNAseq <- read_df(RNAseqplusfile, RNAseqminusfile, "combined.mean")
RNAseq.class <- merge(RNAseq, features, by = c("strand", "TNE"))

RNAseq.plot <- RNAseq.class %>% filter(class == 1) %>% group_dist(., "combined.mean")

pdf("bidir.RNAseq.pdf")
make_graph(RNAseq.plot, 0.2)
dev.off()

############## 2. Total CTSS CAGE ###################
totalplusfile = paste0(dir, "CAGE.FANTOM5.total.fwd.500bp.tab")
totalminusfile = paste0(dir, "CAGE.FANTOM5.total.rev.minus.500bp.tab")

total <- read_df(totalplusfile, totalminusfile, "total.ctss")
total.class <- merge(total, features, c("strand", "TNE"))

total.plot <- total.class %>% filter(class == 1) %>% group_dist(., "total.ctss")

pdf("bidir.total.ctss.pdf")
make_graph(total.plot, 0.6)
dev.off()


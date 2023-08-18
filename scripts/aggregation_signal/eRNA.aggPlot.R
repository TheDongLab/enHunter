####### aggregation plot #######
####### AMPPD eRNA #############

## run on eristwo


dir="/data/bioinformatics/projects/donglab/AMPPD_eRNA/binned/"
marks=c("combined.mean.normalized.random.samplesN200.minus.bigwig", 
        "combined.mean.normalized.random.samplesN200.plus.bigwig", 
        "PhyloP", "TFBS", "DNase",
        "RoadmapBlood.hg19.H3K27acMerged.bigwig", 
        "RoadmapBlood.hg19.H3K4me1Merged.bigwig", 
        "RoadmapBlood.hg19.H3K4me3Merged.bigwig", 
        "CAGE.FANTOM5.total.rev.bigwig", 
        "CAGE.FANTOM5.total.fwd.bigwig")


PLUS=data.frame()
MINUS=data.frame()
for(i in 1:length(marks)) {
  binfile = paste0(dir, marks[i])
  plusfile = paste0(binfile, ".plus.tab")
  minusfile = paste0(binfile, ".minus.tab")
  message(paste("Loading:",binfile))
  
  plus0=read.table(plusfile, check.names =F, stringsAsFactors=F, header=F)[, c(1,5)]; 
  colnames(plus0)=c("names", marks[i])
  
  minus0=read.table(minusfile, check.names =F, stringsAsFactors=F, header=F)[, c(1,5)]; 
  colnames(minus0)=c("names", marks[i])
  
  #rownames(x0)=x0[,1]; x0=x0[,-1];
  if(ncol(PLUS)>0) {
    PLUS=merge(PLUS,plus0);
    MINUS=merge(MINUS,minus0);
  } else { 
    PLUS=plus0;
    MINUS=minus0;
  }
}

dnase <- paste0(dir, "DNase", ".plus.tab")
plus_dnase <- read.table(dnase, check.names =F, stringsAsFactors=F, header=F)[, c(1,5)]

# add TNE and distance to mid point columns
library(stringr)
PLUS[c('chr', 'start', 'end', 'distance')] <- str_split_fixed(PLUS$names, '_', 4)
PLUS['TNE'] <- paste(PLUS$chr, PLUS$start, PLUS$end, sep="_")


MINUS[c('chr', 'start', 'end', 'distance')] <- str_split_fixed(MINUS$names, '_', 4)
MINUS['TNE'] <- paste(MINUS$chr, MINUS$start, MINUS$end, sep="_")


PLUS$strand <- "+"
# drop CAGE and RNAseq minus cols 
PLUS <- PLUS[, !names(PLUS) %in% c("CAGE.FANTOM5.total.rev.bigwig",
                                   "combined.mean.normalized.random.samplesN200.minus.bigwig")]

# rename columns
colnames(PLUS)[colnames(PLUS) == "CAGE.FANTOM5.total.fwd.bigwig"] ="CAGE"
colnames(PLUS)[colnames(PLUS) == "combined.mean.normalized.random.samplesN200.plus.bigwig"] ="RNAseq"


MINUS$strand <- "-"
# drop CAGE and RNAseq plus cols 
MINUS <- MINUS[, !names(MINUS) %in% c("CAGE.FANTOM5.total.fwd.bigwig",
                                      "combined.mean.normalized.random.samplesN200.plus.bigwig")]

# rename columns
colnames(MINUS)[colnames(MINUS) == "CAGE.FANTOM5.total.rev.bigwig"] ="CAGE"
colnames(MINUS)[colnames(MINUS) == "combined.mean.normalized.random.samplesN200.minus.bigwig"] ="RNAseq"


# combine both plus and minus into one df 
df <- rbind(MINUS, PLUS)


#############################################
####### read in classification table ########
#############################################

features = read.table("/data/bioinformatics/projects/donglab/AMPPD_eRNA/eRNA.characterize.feature.color.xls", 
                      check.names =F, stringsAsFactors=F, header=T, row.names = NULL)

colnames(features)[colnames(features) == "row.names"] ="TNE"
features = features[, c("TNE", "class", "features", "strand")]

X = merge(df, features, by = c("strand", "TNE") )
X$distance <- as.numeric(X$distance)

rm(PLUS)
rm(plus0)
rm(MINUS)
rm(minus0)
rm(df)

### start plotting ###
library(dplyr)
library(ggplot2)

pdf("aggplot.pdf")
allplots <- list()
for (i in colnames(X)[4:11]) {
  print(i)
  test <- X %>% group_by(class, distance) %>% summarise(avg = mean(get(i), trim =0.01))
  all <- X %>% group_by(distance) %>% summarise(avg = mean(get(i), trim=0.01))
  all$class <- "none"
  
  final <- rbind(as.data.frame(test), as.data.frame(all))
  plot <- ggplot(final, aes(x = distance, y = avg, color = as.factor(class))) + 
    geom_line() + geom_vline(xintercept = 10.5, linetype="dashed") + 
    scale_color_manual(values=c('#ff3333', '#ff9966', '#ffcccc', 'black')) + 
    theme_bw() + ylab(i)
  allplots[[i]] <- plot
  print(plot)
  
}
dev.off() 









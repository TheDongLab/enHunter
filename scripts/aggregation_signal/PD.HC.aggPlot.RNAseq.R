
# files
#"combined.mean.normalized.control.samplesN200.minus.bigwig",
#"combined.mean.normalized.PD.samplesN200.minus.bigwig", 
#"combined.mean.normalized.PD.samplesN200.plus.bigwig", 
#"combined.mean.normalized.control.samplesN200.plus.bigwig"

dir="/data/bioinformatics/projects/donglab/AMPPD_eRNA/binned/"
marks=c("combined.mean.normalized.control.samplesN200.bigwig",
        "combined.mean.normalized.PD.samplesN200.bigwig")

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

# add TNE and distance to mid point columns
library(stringr)
PLUS[c('chr', 'start', 'end', 'distance')] <- str_split_fixed(PLUS$names, '_', 4)
PLUS['TNE'] <- paste(PLUS$chr, PLUS$start, PLUS$end, sep="_")


MINUS[c('chr', 'start', 'end', 'distance')] <- str_split_fixed(MINUS$names, '_', 4)
MINUS['TNE'] <- paste(MINUS$chr, MINUS$start, MINUS$end, sep="_")


PLUS$strand <- "+"
MINUS$strand <- "-"

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
class1.2.only <- X %>% filter(class != "3")

rm(PLUS)
rm(plus0)
rm(MINUS)
rm(minus0)
rm(df)

### convert to long data frame 
library(reshape2)
class1.2.only <- select(class1.2.only, -chr, -start, -end, -names)

long.df <- melt(class1.2.only, id = c("strand", "TNE", "class", "features", "distance"))
long.df <- mutate(long.df, variable= ifelse(grepl("control", variable, fixed = T ), "HC", "PD"))

### start plotting ###
library(dplyr)
library(ggplot2)
grouped.df <- long.df %>% group_by(distance, variable) %>% summarise(avg = mean(value, trim =0.01))

pdf("PD.HC.aggplot.RNAseq.pdf")
ggplot(grouped.df, aes(x = distance, y = avg, color = variable)) + 
  geom_line() + geom_vline(xintercept = 10.5, linetype="dashed") + 
  theme_bw()

dev.off() 


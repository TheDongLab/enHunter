library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)

#how many TNES have each feature 
# one feature, two feature, ...
# TNE supported by x number features 
# what are the features that support TNE presence?
setwd("~/Documents/college/dong_lab/code/playground/annotations")

# class 1 : DHS site and at least one of five additional external features 
# class 2 : external features but not DHS evidence 
# class 3 : no supporting external features 

###### MINUS STRAND ######
eRNA_minus <- fread("eRNA.minus.characterize.xls")

#### class 1
class1_eRNA_minus <- eRNA_minus[class==1]

print(paste0("Percentage of class 1 eRNAs : ", nrow(class1_eRNA_minus) / nrow(eRNA_minus))) # 29%

#### class 2
class2_eRNA_minus <- eRNA_minus[class==2]

print(paste0("Percentage of class 2 eRNAs : ", nrow(class2_eRNA_minus) / nrow(eRNA_minus))) #20.8%

# enhancers defined by one or more genomic or epigenetic features 
print(paste0("Percentage of class 2 eRNAs : ", 
             (nrow(class1_eRNA_minus) + nrow(class2_eRNA_minus))/ nrow(eRNA_minus))) #49.8%

####### PLUS STRAND ######
eRNA_plus <- fread("eRNA.plus.characterize.xls")

#### class 1
class1_eRNA_plus <- eRNA_plus[class==1]

print(paste0("Percentage of class 1 eRNAs : ", nrow(class1_eRNA_plus) / nrow(eRNA_plus))) # 29%

#### class 2
class2_eRNA_plus <- eRNA_plus[class==2]

print(paste0("Percentage of class 2 eRNAs : ", nrow(class2_eRNA_plus) / nrow(eRNA_plus))) #20.8%

# enhancers defined by one or more genomic or epigenetic features 
print(paste0("Percentage of feature supported eRNAs : ", 
             (nrow(class1_eRNA_plus) + nrow(class2_eRNA_plus))/ nrow(eRNA_plus))) #49.9%

###### graphs #######
# class 1 vs 2 vs 3 enhancers (as defined in Dong et al. 2018)
total_num <- nrow(eRNA_plus) + nrow(eRNA_minus)

class_df <- data.frame(classes = c("class1", "class2", "class3"), 
                       num_TNE = c(nrow(class1_eRNA_plus) + nrow(class1_eRNA_minus), 
                                   nrow(class2_eRNA_plus) + nrow(class2_eRNA_minus), 
                                   nrow(eRNA_plus[class==3]) + nrow(eRNA_minus[class==3])))

class_df <- class_df %>% mutate(percentage = paste0(round((num_TNE / total_num) * 100, digits=2), "%"))

class_graph <- ggplot(class_df, aes(x=classes, y=num_TNE, label=percentage)) + 
  geom_bar(stat="identity") + 
  theme(axis.line = element_line(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        panel.border = element_blank()) + ylab("Number of Plus and Minus TNEs") + 
  xlab("Class Classification") + geom_label()

class_graph

ggsave("class_graph.png", class_graph)

# combining plus and minus data frames together 
eRNA_plus <- eRNA_plus[, strand:="+"]
nrow(eRNA_plus) # 110599
eRNA_minus <- eRNA_minus[, strand:="-"]
nrow(eRNA_minus) # 111368
eRNA <- rbind(eRNA_plus, eRNA_minus)
nrow(eRNA) # 221967

# 110599 + 111368 = 221967 
# no rows are lost by merging 

## feature 1: distance to TSS 
#TODO this graph doesn't represent the data very well 
# NOTE: this can not be log transformed because the negative values will be lost 
dis2TSS <- ggplot(eRNA, aes(x=f01.dis2TSS)) + 
  geom_histogram(fill="gray", binwidth=900) +  
  theme(axis.line = element_line(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        panel.border = element_blank())
dis2TSS
ggsave("dis2TSS.png", dis2TSS)

## feature 2: RPKM 
RPKM <- ggplot(eRNA, aes(x=f02.RPKM)) + 
  geom_histogram(fill="gray", binwidth=0.01) +  
  theme(axis.line = element_line(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        panel.border = element_blank()) + scale_x_continuous(trans="log10")
RPKM
ggsave("RPKM.png", RPKM)

## feature 3: RPM 
RPM <- ggplot(eRNA, aes(x=f03.RPM)) + 
  geom_histogram(fill="gray", binwidth=0.02) +  
  theme(axis.line = element_line(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        panel.border = element_blank()) + scale_x_log10()
RPM
ggsave("RPM.png", RPM)

## feature 13: phyloP
phyloP <- ggplot(eRNA, aes(x=f13.phyloP)) + 
  geom_histogram(fill="gray", binwidth=0.05) +  
  theme(axis.line = element_line(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        panel.border = element_blank()) 
phyloP
ggsave("phyloP.png", phyloP)

## feature 5: CpG score ------ IN PROGRESS -----------------

#NOTE idk the merge function in eRNA.merge didn't work 
#all(is.na(eRNA[["f05.CpG"]])) # every CpG value is NA... 

all_cpgs <- list()
for (i in list("eRNA", "promoters", "random") ) {
  minus <- paste0(i, ".minus.f05.CpG.txt")
  plus <- paste0(i, ".plus.f05.CpG.txt")
  
  cpg_minus <- fread(minus) 
  cpg_minus <- cpg_minus[, strand:="-"]
  
  cpg_plus <- fread(plus)
  cpg_plus <- cpg_plus[, strand:="+"]
  
  df <- rbind(cpg_minus, cpg_plus)
  
  all_cpgs[[i]] <- df
}

cpg_erna <- rbind(cpg_erna_minus, cpg_erna_plus)

## feature annotations broken down by classes 
# features = TFBS, P300, CAGEbloodonlyenhancer, CAGEenhancer.1n2, 
# chromHMM_blood, VISTA, DNaseENCODE, DNaseROADMAP, HCNE
# NOTE: "f08.CAGEenhancer.1n2", "f10.VISTA", "f12.DNase.ENCODE" weren't used in classification  

# select columns
features <- c("f06.TFBS", "f07.P300", 
              "f08.CAGEbloodonlyenhancer", 
              "f08.CAGEenhancer.1n2", "f09.chromHMM_blood", "f10.VISTA", "f12.DNaseENCODE", 
              "f12.DNaseROADMAP", "f15.HCNE", "strand")

# make values binary 
eRNA_features <- eRNA %>% select(features) %>% 
  mutate(f06.TFBS=ifelse(f06.TFBS>5, 1, 0)) %>% 
  mutate(f10.VISTA=ifelse(is.na(f10.VISTA),0, 1))

# supported by at least one occurrence of the feature = 1 
# supported by no features = 0
eRNA_features[eRNA_features>0]=1

eRNA_features$class <- eRNA$class

names(eRNA_features) <- gsub("f[0-9]+.", "", names(eRNA_features))

#transform data from wide to long 
long_eRNA <- eRNA_features %>% pivot_longer(cols=TFBS:HCNE, names_to="Feature", values_to="Value")

#calculate the values in each group 
grouped_long_eRNA <- long_eRNA %>% group_by(class, strand, Feature) %>% summarize(count = sum(Value))
grouped_long_eRNA$class <- as.factor(grouped_long_eRNA$class)

# fill=class 
# x=strand
# facet_grid( ~ Feature)
feature_annot <- ggplot(grouped_long_eRNA, aes(x=strand, y=count, fill=class)) + 
  geom_bar(stat="identity", position = "stack") + facet_grid( ~ Feature) + theme_bw()

feature_annot

ggsave("feature_annot.png", plot= feature_annot)

## look at the TNE length distribution of different classes 
dir <- "~/Documents/college/dong_lab/code/playground"
plus.length <- fread(paste0(dir, "/plus.length.bed"))
minus.length <- fread(paste0(dir, "/minus.length.bed"))

length.df <- rbind( merge(eRNA_plus, plus.length, by="V1"), 
                    merge(eRNA_minus, minus.length, by="V1") ) %>% 
  select(class, V2) %>% 
  mutate(class = as.factor(class))

class_lengths <- ggplot(length.df, aes(x=V2, fill=class)) + 
  geom_histogram(binwidth=0.05) + facet_grid(~class) + 
  scale_x_log10() + xlab("TNE Length")

class_lengths
ggsave("class_lengths.png", plot= class_lengths)

## bidirectional orientation of highly supported TNEs 
wd <- "~/Documents/college/dong_lab/code/playground/peaks/attempt-3"
minus.plus <- fread(paste0(wd, "/inputs/minus.plus.bed")) %>% select(V4, V14, V21)

# are there TNEs on the plus strand that have the exact same coordinates as on the minus strand? 
test <- merge( eRNA, minus.plus, by.x="V1", by.y="V4") %>% mutate(class=as.factor(class))

minus_plus_bidir_class <- ggplot(test, aes(x=V21, fill=class)) + 
  geom_histogram(bins = 200) + facet_grid(~class) +
  scale_x_continuous(limits=c(-1000,1000)) + xlab("Distance to Nearest Peak") + 
  ggtitle("Minus-Plus Orientation Peak Distance") 

bidir_class 
ggsave("minus_plus_bidir_class.png", plot= minus_plus_bidir_class)

## look at number of features supporting TNE elements
#f06.TFBS>=5 to be considered one feature 

# similar to the code for generating the eRNA_features table
feature_ct <- eRNA_features
feature_ct$ID <- eRNA$V1

# data is already complied in eRNA_features table
feature_ct$num <- rowSums(feature_ct[, c("TFBS", "P300", "CAGEbloodonlyenhancer", 
                                                 "CAGEenhancer.1n2", "chromHMM_blood", "VISTA", "DNaseENCODE", 
                                                 "DNaseROADMAP", "HCNE")])

# select for the TNE ids with the most number of features supporting it 
max_num_feature_ids <- feature_ct %>% filter(num == max(feature_ct$num) ) %>% select(ID, strand)
# maximum number of features is 8

# continuing data cleaning ... 
feature_ct$num <- as.factor(feature_ct$num)
#feature_ct <- feature_ct %>% group_by(num) %>% summarise(count=n()) %>% mutate(per = count/sum(count))

TNE_feature_count <- ggplot(feature_ct, aes(x=num)) + geom_bar() + xlab("Number of Features")
TNE_feature_count

ggsave("TNE_feature_count.png", TNE_feature_count)

### number of GWAS SNPS 

eRNA_GWAS <- ggplot(eRNA, aes(x=f16.GWAS)) + 
  geom_bar() + facet_grid(~class) + scale_y_continuous(trans="log10") + 
  xlab("Number of GWAS Snps") +   ggtitle("Distribution of GWAS snps localized to TNE") 

ggsave("TNE_GWAS.png", eRNA_GWAS)

### number of eQTL SNPs

eRNA_gtexDapg <- ggplot(eRNA, aes(x=f18.eSNP.gtexDapg)) + 
  geom_bar() + facet_grid(~class) + scale_y_continuous(trans="log10") + 
  xlab("Number of gtexDapg Snps") +   ggtitle("Distribution of gtexDapg snps localized to TNE") 
eRNA_gtexDapg 

max(eRNA$f18.eSNP.gtexDapg) # 2248

ggsave("TNE_gtexDapg.png", eRNA_gtexDapg)

eRNA_gtexCaviar <- ggplot(eRNA, aes(x=f18.eSNP.gtexCaviar)) + 
  geom_bar() + facet_grid(~class) + scale_y_continuous(trans="log10") + 
  xlab("Number of gtexDapg Snps") +   ggtitle("Distribution of gtexCaviar snps localized to TNE") 
eRNA_gtexCaviar

max(eRNA$f18.eSNP.gtexCaviar) # 107

ggsave("TNE_gtexCaviar.png", eRNA_gtexCaviar)

eRNA_pval <- ggplot(eRNA, aes(x=f18.eSNP.pval)) + 
  geom_bar() + facet_grid(~class) + scale_y_continuous(trans="log10") + 
  xlab("Number of gtexDapg Snps") +   ggtitle("Distribution of pval snps localized to TNE") 
eRNA_pval

max(eRNA$f18.eSNP.pval) # 18384

ggsave("TNE_pval.png", eRNA_pval)


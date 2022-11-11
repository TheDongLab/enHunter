library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)

#how many TNES have each feature 
# one feature, two feature, ...
# TNE supported by x number features 
# what are the features that support TNE presence?

# class 1 : DHS site and at least one of five additional external features 
# class 2 : external features but not DHS evidence 
# class 3 : no supporting external features 

wd <- "./input_files/characterization/"

outdir <- "./scripts/annotations/"

###### MINUS STRAND ######
eRNA_minus <- fread("./input_files/characterization/eRNA.minus.characterize.xls")

#### class 1
class1_eRNA_minus <- eRNA_minus[class==1]

print(paste0("Percentage of class 1 eRNAs : ", nrow(class1_eRNA_minus) / nrow(eRNA_minus))) # 27.8%

#### class 2
class2_eRNA_minus <- eRNA_minus[class==2]

print(paste0("Percentage of class 2 eRNAs : ", nrow(class2_eRNA_minus) / nrow(eRNA_minus))) #19.8%

# enhancers defined by one or more genomic or epigenetic features 
print(paste0("Percentage of class 2 and class 1 eRNAs : ", 
             (nrow(class1_eRNA_minus) + nrow(class2_eRNA_minus))/ nrow(eRNA_minus))) #47.6

####### PLUS STRAND ######
eRNA_plus <- fread("./input_files/characterization/eRNA.plus.characterize.xls")

#### class 1
class1_eRNA_plus <- eRNA_plus[class==1]

print(paste0("Percentage of class 1 eRNAs : ", nrow(class1_eRNA_plus) / nrow(eRNA_plus))) # 21%

#### class 2
class2_eRNA_plus <- eRNA_plus[class==2]

print(paste0("Percentage of class 2 eRNAs : ", nrow(class2_eRNA_plus) / nrow(eRNA_plus))) #15.9%

# enhancers defined by one or more genomic or epigenetic features 
print(paste0("Percentage of feature supported eRNAs : ", 
             (nrow(class1_eRNA_plus) + nrow(class2_eRNA_plus))/ nrow(eRNA_plus))) #37.1%

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

ggsave("./scripts/annotations/class_graph.png", class_graph)

## class annotations split by strand 
class_df_split <- data.frame(class = c("class1", "class2", "class3", "class1", "class2", "class3"), 
                             num_TNE = c(nrow(class1_eRNA_minus), nrow(class2_eRNA_minus), nrow(eRNA_minus[class==3]), 
                                         nrow(class1_eRNA_plus), nrow(class2_eRNA_plus), nrow(eRNA_plus[class==3])), 
                             total_TNE = c(nrow(eRNA_minus), nrow(eRNA_minus), nrow(eRNA_minus), nrow(eRNA_plus), nrow(eRNA_plus), nrow(eRNA_plus)),
                             strand = c("-", "-", "-", "+", "+", "+"))

class_df_split <- class_df_split %>% mutate(percentage = paste0(round((num_TNE / total_TNE) * 100, digits=2), "%"))

class_graph_split <- ggplot(class_df_split, aes(x=strand, y=num_TNE, label=percentage, fill=strand, alpha=class)) + 
  geom_bar(stat="identity", position = position_dodge()) + 
  theme(axis.line = element_line(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        panel.border = element_blank(), text = element_text(size = 20)) + ylab("TNE Count") + 
  geom_text(position = position_dodge(width = 0.9), vjust = 0) + 
  xlab("") + scale_y_continuous(breaks = seq(0, 70000, by= 5000)) + scale_alpha_manual(values = c(.25, .5, .75, 1)) + 
  scale_fill_manual(values = c("blue", "red"))

class_graph_split
ggsave("./scripts/annotations/class_graph_split.png", class_graph_split)


# combining plus and minus data frames together 
eRNA_plus <- eRNA_plus[, strand:="+"]
nrow(eRNA_plus) # 66782
eRNA_minus <- eRNA_minus[, strand:="-"]
nrow(eRNA_minus) # 42815
eRNA <- rbind(eRNA_plus, eRNA_minus)
nrow(eRNA) # 109597

# 66782 + 42815 = 109597
# no rows are lost by merging 

## feature 1: distance to TSS 
# NO TNEs have a distance of 0 to TSS

intergenic_minus <- round((nrow(eRNA_minus[f01.dis2TSS < 0]) / nrow(eRNA_minus) ) * 100, digits = 2)
intergenic_txt_minus <- paste0("intergenic (N=", nrow(eRNA_minus[f01.dis2TSS < 0]), "; ", intergenic_minus, "%)")

intronic_minus <- round((nrow(eRNA_minus[f01.dis2TSS > 0]) / nrow(eRNA_minus) ) * 100, digits = 2)
intronic_txt_minus <- paste0("intronic (N=", nrow(eRNA_minus[f01.dis2TSS > 0]), "; ", intronic_minus, "%)")

dis2TSS_minus <- ggplot(eRNA_minus, aes(x=f01.dis2TSS/1000, fill=f01.dis2TSS < 0)) + 
  geom_histogram( binwidth=20) +  
  theme(axis.line = element_line(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        panel.border = element_blank(), legend.title=element_blank(), legend.position=c(0.3, 0.9), text = element_text(size = 20)) +
  xlab("Distance to the TSS (x 1000 bp)") + 
  scale_x_continuous(n.breaks = 6) + scale_fill_manual(values = c("purple", "green"), labels = c(intronic_txt_minus, intergenic_txt_minus))
dis2TSS_minus
ggsave("./scripts/annotations/dis2TSS_minus.png", dis2TSS_minus)


intergenic_plus <- round((nrow(eRNA_plus[f01.dis2TSS < 0]) / nrow(eRNA_plus) ) * 100, digits = 2)
intergenic_txt_plus <- paste0("intergenic (N=", nrow(eRNA_plus[f01.dis2TSS < 0]), "; ", intergenic_plus, "%)")

intronic_plus <- round((nrow(eRNA_plus[f01.dis2TSS > 0]) / nrow(eRNA_plus) ) * 100, digits = 2)
intronic_txt_plus <- paste0("intronic (N=", nrow(eRNA_plus[f01.dis2TSS > 0]), "; ", intronic_plus, "%)")

dis2TSS_plus <- ggplot(eRNA_plus, aes(x=f01.dis2TSS/1000, fill=f01.dis2TSS < 0)) + 
  geom_histogram( binwidth=20) +  
  theme(axis.line = element_line(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        panel.border = element_blank(), legend.title=element_blank(), legend.position=c(0.2, 0.7), text = element_text(size = 20)) + xlab("Distance to the TSS (x 1000 bp)") + 
  scale_x_continuous(n.breaks = 6) + scale_fill_manual(values = c("purple", "green"), labels = c(intronic_txt_plus, intergenic_txt_plus))
dis2TSS_plus
ggsave("./scripts/annotations/dis2TSS_plus.png", dis2TSS_plus)

## feature 2: RPKM 

minus.RPKM.med <- median(eRNA_minus[, f02.RPKM], na.rm = TRUE)

RPKM_minus <- ggplot(eRNA_minus, aes(x=f02.RPKM)) + 
  geom_histogram(fill="gray", binwidth=0.01) +  
  theme(axis.line = element_line(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        panel.border = element_blank(), text = element_text(size = 20)) + 
  geom_vline(xintercept=minus.RPKM.med, color="red", linetype="dashed") +
  scale_x_continuous(trans="log10") + xlab("RPKM") + ylab("Count of TNEs") + 
  annotate("text", x = minus.RPKM.med+240, y=8000, size = 6,
           label=paste0("median: ", round(minus.RPKM.med, digits = 2)), color="red")
RPKM_minus
ggsave("./scripts/annotations/RPKM_minus.png", RPKM_minus)

plus.RPKM.med <- median(eRNA_plus[, f02.RPKM], na.rm = TRUE)

RPKM_plus <- ggplot(eRNA_plus, aes(x=f02.RPKM)) + 
  geom_histogram(fill="gray", binwidth=0.01) +  
  theme(axis.line = element_line(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        panel.border = element_blank(), text = element_text(size = 20)) + 
  geom_vline(xintercept=plus.RPKM.med, color="red", linetype="dashed") +
  scale_x_continuous(trans="log10") + xlab("RPKM") + ylab("Count of TNEs") + 
  annotate("text", x = plus.RPKM.med+240, y=8000, size = 6,
           label=paste0("median: ", round(plus.RPKM.med, digits = 2)), color="red")
RPKM_plus
ggsave("./scripts/annotations/RPKM_plus.png", RPKM_plus)


## feature 3: RPM

minus.RPM.med <- median(eRNA_minus[, f03.RPM], na.rm = T)

RPM_minus <- ggplot(eRNA_minus, aes(x=f03.RPM)) + 
  geom_histogram(fill="gray", binwidth=0.02) +  
  theme(axis.line = element_line(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        panel.border = element_blank(), text = element_text(size = 20)) + 
  scale_x_log10() + geom_vline(xintercept=minus.RPM.med, color="red", linetype="dashed") +
  xlab("RPM") + ylab("Count of TNEs") + 
  annotate("text", x = minus.RPM.med+1, y=5000, size = 6,
           label=paste0("median: ", round(minus.RPM.med, digits = 2)), color="red")
RPM_minus
ggsave("./scripts/annotations/RPM_minus.png", RPM_minus)

plus.RPM.med <- median(eRNA_plus[, f03.RPM], na.rm = T)

RPM_plus <- ggplot(eRNA_plus, aes(x=f03.RPM)) + 
  geom_histogram(fill="gray", binwidth=0.02) +  
  theme(axis.line = element_line(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        panel.border = element_blank(), text = element_text(size = 20)) + 
  scale_x_log10() +
  geom_vline(xintercept=plus.RPM.med, color="red", linetype="dashed") +
  xlab("RPM") + ylab("Count of TNEs") + 
  annotate("text", x = plus.RPM.med+1, y=5000, size = 6,
           label=paste0("median: ", round(plus.RPM.med, digits = 2)), color="red")
RPM_plus
ggsave("./scripts/annotations/RPM_plus.png", RPM_plus)


## feature 13: phyloP
phyloP <- ggplot(eRNA, aes(x=f13.phyloP)) + 
  geom_histogram(fill="gray", binwidth=0.05) +  
  theme(axis.line = element_line(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        panel.border = element_blank()) 
phyloP
ggsave("./scripts/annotations/phyloP.png", phyloP)

## feature 5: CpG score ------ IN PROGRESS -----------------

#NOTE idk the merge function in eRNA.merge didn't work 
#all(is.na(eRNA[["f05.CpG"]])) # every CpG value is NA... 

# all_cpgs <- list()
# for (i in list("eRNA", "promoters", "random") ) {
#   minus <- paste0(i, ".minus.f05.CpG.txt")
#   plus <- paste0(i, ".plus.f05.CpG.txt")
#   
#   cpg_minus <- fread(minus) 
#   cpg_minus <- cpg_minus[, strand:="-"]
#   
#   cpg_plus <- fread(plus)
#   cpg_plus <- cpg_plus[, strand:="+"]
#   
#   df <- rbind(cpg_minus, cpg_plus)
#   
#   all_cpgs[[i]] <- df
# }
# 
# cpg_erna <- rbind(cpg_erna_minus, cpg_erna_plus)

## feature annotations broken down by classes 
# features = TFBS, P300, CAGEbloodenhancer, CAGEenhancer.1n2, 
# chromHMM_blood, VISTA, DNaseENCODE, DNaseROADMAP, HCNE
# NOTE: "f08.CAGEenhancer.1n2", "f10.VISTA", "f12.DNase.ENCODE" weren't used in classification  

# select columns
features <- c("f06.TFBS", "f07.P300", 
              "f08.CAGEbloodenhancer", 
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
  geom_bar(stat="identity", position = "stack") + facet_wrap( ~ Feature, scales = "free_y") +
  theme_bw() + theme( text = element_text(size = 15))

feature_annot

ggsave("./scripts/annotations/feature_annot.png", plot= feature_annot)

## look at the TNE length distribution of different classes 
dir <- "./input_files/TNE"
plus.length <- fread(paste0(dir, "/plus.length.bed"))
minus.length <- fread(paste0(dir, "/minus.length.bed"))

length.df <- rbind( merge(eRNA_plus, plus.length, by="V1"), 
                    merge(eRNA_minus, minus.length, by="V1") ) %>% 
  select(class, V2) %>% 
  mutate(class = as.factor(class))

class_lengths <- ggplot(length.df, aes(x=V2, fill=class)) + 
  geom_histogram(binwidth=0.05) + facet_wrap(~class, scales = "free_y") + 
  scale_x_log10() + xlab("TNE Length") 

class_lengths
ggsave("./scripts/annotations/plus_minus_class_lengths.png", plot= class_lengths)

## bidirectional orientation of highly supported TNEs 
wd <- "./input_files/closest/"
minus.plus <- fread(paste0(wd, "minus.plus.bed")) %>% select(V4, V14, V21)

# are there TNEs on the plus strand that have the exact same coordinates as on the minus strand? 
test <- merge( eRNA, minus.plus, by.x="V1", by.y="V4") %>% mutate(class=as.factor(class))

minus_plus_bidir_class <- ggplot(test, aes(x=V21, fill=class)) + 
  geom_histogram(bins = 200) + facet_grid(~class, scales = "free_y") +
  scale_x_continuous(limits=c(-1000,1000)) + xlab("Distance to Nearest Peak") + 
  ggtitle("Minus-Plus Orientation Peak Distance") 

minus_plus_bidir_class

ggsave("./scripts/annotations/minus_plus_bidir_class.png", plot= minus_plus_bidir_class)

## look at number of features supporting TNE elements
#f06.TFBS>=5 to be considered one feature 

# similar to the code for generating the eRNA_features table
feature_ct <- eRNA_features
feature_ct$ID <- eRNA$V1
feature_ct$vals <- apply(feature_ct, 1, function(x){
  x <- x[c("TFBS", "P300", "CAGEbloodenhancer", 
             "CAGEenhancer.1n2", "chromHMM_blood", "VISTA", "DNaseENCODE", 
             "DNaseROADMAP", "HCNE")]
  paste0(unlist(  unlist(names(x[x==1]))), collapse=" ,")

  })

# data is already complied in eRNA_features table
feature_ct$num <- rowSums(feature_ct[, c("TFBS", "P300", "CAGEbloodenhancer", 
                                                 "CAGEenhancer.1n2", "chromHMM_blood", "VISTA", "DNaseENCODE", 
                                                 "DNaseROADMAP", "HCNE")])

# select for the TNE ids with the most number of features supporting it 
max_num_feature_ids <- feature_ct %>% filter(num == max(feature_ct$num) ) %>% 
  select(ID, strand, class, vals) %>% mutate(num_features= max(feature_ct$num))
# maximum number of features is 8

testing <- feature_ct %>% filter(num == max(feature_ct$num) - 1) %>% 
  select(ID, strand, class, vals) %>% mutate(num_features= max(feature_ct$num) - 1)

max_num_feature_ids <- rbind(max_num_feature_ids, testing)

# add RPM and RPKM info 
max_num_feature_ids <- merge(max_num_feature_ids, eRNA[, c("V1", "strand", "f03.RPM", "f02.RPKM")], 
                             by.x = c("ID", "strand"), by.y = c("V1", "strand")) 

write.table(max_num_feature_ids , file = "top_TNEs.tsv", row.names = FALSE, sep = "\t")

# continuing data cleaning ... 
feature_ct$num <- as.factor(feature_ct$num)
counts <- feature_ct %>% group_by(num) %>% summarise(count = n())

feature_ct <- merge(feature_ct, counts)

#feature_ct <- feature_ct %>% group_by(num) %>% summarise(count=n()) %>% mutate(per = count/sum(count))

TNE_feature_count <- ggplot(feature_ct, aes(x= num)) + 
  geom_bar() + xlab("Number of Features") + 
  theme(axis.line = element_line(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        panel.border = element_blank(), text = element_text(size = 20)) + 
  scale_y_continuous(n.breaks=14)
TNE_feature_count

ggsave("./scripts/annotations/TNE_feature_count.png", TNE_feature_count)

### number of GWAS SNPS 

eRNA_GWAS <- ggplot(eRNA, aes(x=f16.GWAS)) + 
  geom_bar() + facet_grid(~class) + scale_y_continuous(trans="log10") + 
  xlab("Number of GWAS Snps") +   ggtitle("Distribution of GWAS snps localized to TNE") 

eRNA_GWAS

ggsave("./scripts/annotations/TNE_GWAS.png", eRNA_GWAS)

### number of eQTL SNPs

eRNA_gtexDapg <- ggplot(eRNA, aes(x=f18.eSNP.gtexDapg)) + 
  geom_bar() + facet_grid(~class) + scale_y_continuous(trans="log10") + 
  xlab("Number of gtexDapg Snps") +   ggtitle("Distribution of gtexDapg snps localized to TNE") 
eRNA_gtexDapg 

max(eRNA$f18.eSNP.gtexDapg) # 2248

ggsave("./scripts/annotations/TNE_gtexDapg.png", eRNA_gtexDapg)

eRNA_gtexCaviar <- ggplot(eRNA, aes(x=f18.eSNP.gtexCaviar)) + 
  geom_bar() + facet_grid(~class) + scale_y_continuous(trans="log10") + 
  xlab("Number of gtexDapg Snps") +   ggtitle("Distribution of gtexCaviar snps localized to TNE") 
eRNA_gtexCaviar

max(eRNA$f18.eSNP.gtexCaviar) # 170

ggsave("./scripts/annotations/TNE_gtexCaviar.png", eRNA_gtexCaviar)

eRNA_pval <- ggplot(eRNA, aes(x=f18.eSNP)) + 
  geom_bar() + facet_grid(~class) + scale_y_continuous(trans="log10") + 
  xlab("Number of gtexDapg Snps") +   ggtitle("Distribution of pval snps localized to TNE") 
eRNA_pval

max(eRNA$f18.eSNP) # 18384

ggsave("./scripts/annotations/TNE_pval.png", eRNA_pval)


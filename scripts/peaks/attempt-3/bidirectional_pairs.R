# selecting for bidirectional TNE pairs based on minus.plus.bed and plus.minus.bed 
# (outputs from narrowPeak -> bedtools closest )

# 1. both TNEs are in minus to plus orientation (negative distance)
# 2. summit distance is < 500 bp 
# 3. one of the pairs is class 1 or 2 

##### 1. both TNEs are in minus to plus orientation (negative distance)

# plus.minus.bed 
system("mkdir scripts/peaks/attempt-3/pairs")
system("awk '$21<0 {print $0}' input_files/closest/plus.minus.bed > scripts/peaks/attempt-3/pairs/plus.minus.step1.bed")
system("wc -l scripts/peaks/attempt-3/pairs/plus.minus.step1.bed")
# there are 37613 pairs in minus to plus orientation for plus.minus.bed 

# minus.plus.bed
system("awk '$21<0 {print $0}' input_files/closest/minus.plus.bed > scripts/peaks/attempt-3/pairs/minus.plus.step1.bed")
system("wc -l scripts/peaks/attempt-3/pairs/minus.plus.step1.bed")
# there are 24078 pairs in minus to plus orientation for minus.plus.bed 

##### 2. summit distance is < 500 bp 

# plus.minus.bed 
system("awk '$21>-500 {print $0}' scripts/peaks/attempt-3/pairs/plus.minus.step1.bed > scripts/peaks/attempt-3/pairs/plus.minus.step2.bed")
system("wc -l scripts/peaks/attempt-3/pairs/plus.minus.step2.bed")
# there are 341 pairs < 500 bp apart for plus.minus orientation 

# minus.plus.bed 
system("awk '$21>-500 {print $0}' scripts/peaks/attempt-3/pairs/minus.plus.step1.bed > scripts/peaks/attempt-3/pairs/minus.plus.step2.bed")
system("wc -l scripts/peaks/attempt-3/pairs/minus.plus.step2.bed")
# there are 339 pairs < 500 bp apart for minus.plus.bed 

system("awk '$21>-700 {print $0}' scripts/peaks/attempt-3/pairs/minus.plus.step1.bed > scripts/peaks/attempt-3/pairs/minus.plus.step2-700bp.bed")
system("wc -l scripts/peaks/attempt-3/pairs/minus.plus.step2-700bp.bed")

system("awk '$21>-700 {print $0}' scripts/peaks/attempt-3/pairs/plus.minus.step1.bed > scripts/peaks/attempt-3/pairs/plus.minus.step2-700bp.bed")
system("wc -l scripts/peaks/attempt-3/pairs/plus.minus.step2-700bp.bed")


## merge minus.plus.step2.bed and plus.minus.step2.bed, checking for duplicate pairs 
## reorganize minus.plus.step2.bed into plus minus TNE listing 
system("paste <(cut -f11-20 scripts/peaks/attempt-3/pairs/minus.plus.step2.bed) <(cut -f1-10 scripts/peaks/attempt-3/pairs/minus.plus.step2.bed) <(cut -f21 scripts/peaks/attempt-3/pairs/minus.plus.step2.bed) > scripts/peaks/attempt-3/pairs/plus.minus.step2.tmp")
system("cat scripts/peaks/attempt-3/pairs/plus.minus.step2.tmp scripts/peaks/attempt-3/pairs/plus.minus.step2.bed | sort | uniq > scripts/peaks/attempt-3/pairs/merged.bed")
system("wc -l scripts/peaks/attempt-3/pairs/merged.bed")
# there are 372 pairs in merged.bed 

system("cut -f4,6,14,16,21 scripts/peaks/attempt-3/pairs/merged.bed > scripts/peaks/attempt-3/pairs/bidirectional_pairs.txt")

## 700 bp 
system("paste <(cut -f11-20 scripts/peaks/attempt-3/pairs/minus.plus.step2-700bp.bed) <(cut -f1-10 scripts/peaks/attempt-3/pairs/minus.plus.step2-700bp.bed) <(cut -f21 scripts/peaks/attempt-3/pairs/minus.plus.step2-700bp.bed) > scripts/peaks/attempt-3/pairs/plus.minus.step2-700.tmp")
system("cat scripts/peaks/attempt-3/pairs/plus.minus.step2-700.tmp scripts/peaks/attempt-3/pairs/plus.minus.step2-700bp.bed | sort | uniq > scripts/peaks/attempt-3/pairs/merged-700.bed")
system("wc -l scripts/peaks/attempt-3/pairs/merged-700.bed")

system("cut -f4,6,14,16,21 scripts/peaks/attempt-3/pairs/merged-700.bed > scripts/peaks/attempt-3/pairs/bidirectional_pairs_700.txt")

## no distance filter 
system("paste <(cut -f11-20 input_files/closest/minus.plus.bed) <(cut -f1-10 input_files/closest/minus.plus.bed) <(cut -f21 input_files/closest/minus.plus.bed) > scripts/peaks/attempt-3/pairs/plus.minus.all.tmp")
system("cat scripts/peaks/attempt-3/pairs/plus.minus.all.tmp input_files/closest/plus.minus.bed | sort | uniq > scripts/peaks/attempt-3/pairs/merged_all.bed")
system("wc -l scripts/peaks/attempt-3/pairs/merged_all.bed")

system("cut -f4,6,14,16,21 scripts/peaks/attempt-3/pairs/merged_all.bed > scripts/peaks/attempt-3/pairs/bidirectional_pairs_all.txt")

##### 3. class 1 or class 2 TNEs
library(data.table)
library(dplyr)

eRNA.class <- fread("input_files/characterization/eRNA.characterize.feature.color.xls") %>% select(V1, strand, class)
#bidir_pairs <- fread("scripts/peaks/attempt-3/pairs/bidirectional_pairs.txt", header = FALSE)
#bidir_pairs <- fread("scripts/peaks/attempt-3/pairs/bidirectional_pairs_700.txt", header = FALSE)
bidir_pairs <- fread("scripts/peaks/attempt-3/pairs/bidirectional_pairs_all.txt", header = FALSE)

colnames(bidir_pairs) <- c("TNE_1", "strand_1", "TNE_2", "strand_2", "distance")

TNE1_bidir <- merge(bidir_pairs, eRNA.class, by.x = c("TNE_1", "strand_1"), by.y = c("V1", "strand")) %>% rename(TNE_1_class = class)

final <- merge(TNE1_bidir, eRNA.class, by.x = c("TNE_2", "strand_2"), by.y = c("V1", "strand")) %>% rename(TNE_2_class = class)
# nrow(final) = 372

# 228 pairs were filtered out due to class limitations 

# no class 3 class 3 pairs 
final_bidir_pairs <- final %>% filter( !(TNE_1_class == 3 & TNE_2_class == 3 ) )
# nrow(final_bidir_pairs) = 144
#no_classIII_III_pairs.distance.density.plot.pdf

# no class 3 TNEs 
no_3_bidir_pairs <- final %>% filter( !(TNE_1_class == 3 | TNE_2_class == 3)  )
#classI_or_II.distance.density.plot.pdf

nrow(no_3_bidir_pairs %>% filter(distance > -1000  & distance < 1000))
test <- no_3_bidir_pairs %>% filter(distance > -1000  & distance < 1000)

class1_pair <- final %>% filter( TNE_1_class == 1 | TNE_2_class == 1  )
#single_classI_pairs.distance.density.plot.pdf

#write.table(final_bidir_pairs, file = "scripts/peaks/attempt-3/pairs/class_bidirectional_pairs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)
#write.table(final_bidir_pairs, file = "scripts/peaks/attempt-3/pairs/class_bidirectional_pairs_700.tsv", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(final_bidir_pairs, file = "scripts/peaks/attempt-3/pairs/class_bidirectional_pairs_all.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

system("cat <(cut -f1 class_bidirectional_pairs.tsv | sed '1d') <(cut -f3 class_bidirectional_pairs.tsv | sed '1d') | wc -l  > all_bidirectional_class_TNEs.txt")

#system("cat <(cut -f1 class_bidirectional_pairs_700.tsv | sed '1d') <(cut -f3 class_bidirectional_pairs_700.tsv | sed '1d') | wc -l  > all_bidirectional_class_TNEs_700.txt")


###### graphing the bidirectional distance of pairs 
library(ggplot2)

final_pairs <- fread("scripts/peaks/attempt-3/pairs/class_bidirectional_pairs_all.tsv", sep = "\t")

#no_3_bidir_pairs <- no_3_bidir_pairs %>% filter(TNE_1_class == 1 & TNE_2_class == 1)

dist_density_minus.plus <- ggplot(final_bidir_pairs, aes(x=distance)) + 
  geom_histogram(color="black", fill="white", binwidth = 30) + 
  scale_x_continuous(limits=c(-1000,1000), n.breaks = 10) + 
  theme(axis.line = element_line(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        panel.border = element_blank(), text = element_text(size = 20))  +
  xlab("Distance to Nearest Peak") + 
  ggtitle("Bidirectional Pairs Peak Distance")

dist_density_minus.plus

ggsave("single_classI_pairs.distance.density.plot.pdf",
       plot=dist_density_minus.plus, device="pdf")

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

## merge minus.plus.step2.bed and plus.minus.step2.bed, checking for duplicate pairs 
## reorganize minus.plus.step2.bed into plus minus TNE listing 
system("paste <(cut -f11-20 scripts/peaks/attempt-3/pairs/minus.plus.step2.bed) <(cut -f1-10 scripts/peaks/attempt-3/pairs/minus.plus.step2.bed) <(cut -f21 scripts/peaks/attempt-3/pairs/minus.plus.step2.bed) > scripts/peaks/attempt-3/pairs/plus.minus.step2.tmp")
system("cat scripts/peaks/attempt-3/pairs/plus.minus.step2.tmp scripts/peaks/attempt-3/pairs/plus.minus.step2.bed | sort | uniq > scripts/peaks/attempt-3/pairs/merged.bed")
system("wc -l scripts/peaks/attempt-3/pairs/merged.bed")
# there are 372 pairs in merged.bed 

system("cut -f4,6,14,16,21 scripts/peaks/attempt-3/pairs/merged.bed > scripts/peaks/attempt-3/pairs/bidirectional_pairs.txt")


##### 3. class 1 or class 2 TNEs
library(data.table)
library(dplyr)

eRNA.class <- fread("input_files/characterization/eRNA.characterize.xls") %>% select(V1, strand, class)
bidir_pairs <- fread("scripts/peaks/attempt-3/pairs/bidirectional_pairs.txt", header = FALSE)
colnames(bidir_pairs) <- c("TNE_1", "strand_1", "TNE_2", "strand_2", "distance")

TNE1_bidir <- merge(bidir_pairs, eRNA.class, by.x = c("TNE_1", "strand_1"), by.y = c("V1", "strand")) %>% rename(TNE_1_class = class)

final <- merge(TNE1_bidir, eRNA.class, by.x = c("TNE_2", "strand_2"), by.y = c("V1", "strand")) %>% rename(TNE_2_class = class)
# nrow(final) = 372

# 228 pairs were filtered out due to class limitations 

final_bidir_pairs <- final %>% filter( !(TNE_1_class == 3 & TNE_2_class == 3) )
# nrow(final_bidir_pairs) = 144

write.table(final_bidir_pairs, file = "scripts/peaks/attempt-3/pairs/class_bidirectional_pairs.tsv",
            quote = FALSE, sep = "\t", row.names = FALSE)

# usage:  Rscript bidirectional_pairs.v2.R input_files/closest/plus.minus.bed input_files/closest/minus.plus.bed 900

##### 1. both TNEs are in minus to plus orientation (negative distance)
##### 2. summit distance is < x bp 

args<-commandArgs(trailingOnly=TRUE)

plus.minus <- args[1]
minus.plus <- args[2]
dist <- args[3]

system(paste0("mkdir scripts/peaks/bidirectional/pairs/",dist,"bp") )
system(paste0("cd scripts/peaks/bidirectional/pairs/",dist,"bp") )

##### 1. both TNEs are in minus to plus orientation (negative distance)

# plus.minus.bed 
#input_files/closest/plus.minus.bed
system(paste0("awk '$21<0 {print $0}' ", plus.minus, " > scripts/peaks/bidirectional/pairs/plus.minus.step1.bed") )
system("wc -l scripts/peaks/bidirectional/pairs/plus.minus.step1.bed")
# there are 37613 pairs in minus to plus orientation for plus.minus.bed 

# input_files/closest/minus.plus.bed
system(paste0("awk '$21<0 {print $0}' ", minus.plus, " > scripts/peaks/bidirectional/pairs/minus.plus.step1.bed") )
system("wc -l scripts/peaks/bidirectional/pairs/minus.plus.step1.bed")
# there are 24078 pairs in minus to plus orientation for minus.plus.bed 

##### 2. summit distance is < 500 bp 

# plus.minus.bed 
system(paste0("awk '$21>-",dist," {print $0}' scripts/peaks/bidirectional/pairs/plus.minus.step1.bed > scripts/peaks/bidirectional/pairs/plus.minus.step2.bed") )
system("wc -l scripts/peaks/bidirectional/pairs/plus.minus.step2.bed")

# minus.plus.bed 
system(paste0("awk '$21>-", dist, " {print $0}' scripts/peaks/bidirectional/pairs/minus.plus.step1.bed > scripts/peaks/bidirectional/pairs/minus.plus.step2.bed") )
system("wc -l scripts/peaks/bidirectional/pairs/minus.plus.step2.bed")

## merge minus.plus.step2.bed and plus.minus.step2.bed, checking for duplicate pairs 
## reorganize minus.plus.step2.bed into plus minus TNE listing
system(paste0("cd scripts/peaks/bidirectional/pairs/",dist,"bp") )
system("paste <(cut -f11-20 minus.plus.step2.bed) <(cut -f1-10 minus.plus.step2.bed) <(cut -f21 minus.plus.step2.bed) > plus.minus.step2.tmp")
system("cat plus.minus.step2.tmp plus.minus.step2.bed | sort | uniq > merged.bed")
system("wc -l merged.bed")

system("cut -f4,6,14,16,21 merged.bed > bidirectional_pairs.txt")

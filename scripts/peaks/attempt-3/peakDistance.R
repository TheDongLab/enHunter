library(data.table)
library(ggplot2)
library(dplyr)

### this is different than attempt-2 as it contains the narrowPeak distances of TNEs 
# with 8000 sample family significance in step 6 

## from an atlas of active enhancers across human cell types and tissues 
#Enhancer-associated reverse and forward strand transcription initiation events were, on average, separated by 180 bp 

###### with the minus plus orientation 
minus.plus <- fread("input_files/closest/minus.plus.bed", sep = "\t")

minus.plus.dist <- ggplot(minus.plus, aes(x=V21)) + 
  geom_histogram(color="black", fill="white") + 
  scale_y_continuous(trans="log2")

minus.plus.dist

dist_minus.plus <- ggplot(minus.plus, aes(x=V21, y=..density..)) + 
  geom_histogram(color="black", fill="white", binwidth = 30) + 
  geom_density() + 
  xlab("Distance to Nearest Peak") + 
  ggtitle("Minus-Plus Orientation Peak Distance") 

dist_minus.plus

dist_density_minus.plus <- ggplot(minus.plus, aes(x=V21, y=..density..)) + 
  geom_histogram(color="black", fill="white", binwidth = 5) + 
  scale_x_continuous(limits=c(-100,100)) + 
  geom_density() + 
  xlab("Distance to Nearest Peak") + 
  ggtitle("Minus-Plus Orientation Peak Distance") 

dist_density_minus.plus

ggsave("minus.plus.distance.density.plot.pdf",
       plot=dist_density_minus.plus, device="pdf")
ggsave("minus.plus.distance.density.plot.png",
       plot=dist_density_minus.plus, device="png")


min_dis <- min(minus.plus$V21) # -5511680
max_dis <- max(minus.plus$V21) # 1764470

# no distances are na 
any(is.na(minus.plus$V21))

# only neg distances
neg.minus.plus <- minus.plus[minus.plus$V21 < 0]
nrow(neg.minus.plus) # 24078
length(minus.plus$V21) # 42818

print(paste0("percentage of TNEs where - is before + ", nrow(neg.minus.plus) / length(minus.plus$V21) ))
# ~ 56%

minus.plus.dist <- ggplot(minus.plus) + geom_histogram(aes(x=V21), color="black", fill="white", binwidth = 10) + 
  scale_x_continuous(limits=c(-200,200)) # limit minimum x value to see distribution of values at higher res

minus.plus.dist

median(minus.plus)
####### with the plus minus orientation 

plus.minus <- fread("./inputs/plus.minus.bed", sep = "\t")

plus.minus.dist <- ggplot(plus.minus, aes(x=V21)) + 
  geom_histogram(color="black", fill="white") + 
  scale_y_continuous(trans="log2")

plus.minus.dist

dist_density_plus.minus <- ggplot(plus.minus, aes(x=V21, y=..density..)) + 
  geom_histogram(color="black", fill="white", binwidth = 30) + 
  scale_x_continuous(limits=c(-1000,1000)) + 
  geom_density() + 
  xlab("Distance to Nearest Peak") + ggtitle("Plus-Minus Orientation Peak Distance")

dist_density_plus.minus 

ggsave("plus.minus.distance.density.plot.pdf",
       plot=dist_density_plus.minus, device="pdf")

ggsave("plus.minus.distance.density.plot.png",
       plot=dist_density_plus.minus, device="png")

min_dis <- min(plus.minus$V21) # -1912881
max_dis <- max(plus.minus$V21) # 2815664

# no distances are na 
any(is.na(plus.minus$V21))

# only neg distances

neg.plus.minus <- plus.minus[plus.minus$V21 < 0, ]
print(nrow(neg.plus.minus)) # 34068
length(plus.minus$V21) # 110647

print(paste0("percentage of TNEs where - is before + ", nrow(neg.plus.minus) / length(plus.minus$V21) ))
# ~ 30%

minus.plus.dist <- ggplot() + geom_histogram(aes(x=neg.plus.minus), color="black", fill="white", binwidth = 10) + 
  scale_x_continuous(limits=c(-500,10)) # limit minimum x value to see distribution of values at higher res

length(plus.minus$V21[plus.minus$V21 == 0]) #175

plus.minus %>%
  ggplot(aes(x=V16, y=V21)) +
  geom_jitter(color="black", size=0.4, alpha=0.9) + 
  geom_boxplot() + scale_y_continuous(limits=c(-1000,10))

#+
#geom_jitter(color="black", size=0.4, alpha=0.9) +
  #theme(
    #legend.position="none",
    #plot.title = element_text(size=11)
  #) +
  #ggtitle("A boxplot with jitter") +
  #xlab("")

filtered.plus.minus  <- neg.plus.minus[neg.plus.minus$V21 > -400]
nrow(filtered.plus.minus) # 19,363

filtered.minus.plus  <- neg.minus.plus[neg.minus.plus$V21 > -400]
nrow(filtered.minus.plus) # 19,353


# there are 19,365 which are in reverse forward orientation with a distance of < 400 bp apart

# reformat narrowPeak format into bed format 
bed.filtered.plus.minus <- data.table(chrom=character(), 
                                      start=numeric(), 
                                      end=numeric(), 
                                      name=character(), 
                                      score=numeric(), 
                                      strand=character())

### this was only run on the pairs found on minus strand by using plus strand as a reference 
# aka plus.minus orientation 
for(i in 1:nrow(filtered.plus.minus)) {
  row <- filtered.plus.minus[i, ]
  first_bed <- data.table(chrom=row$V1, start=row$V2, end=row$V3, name=row$V4, score=row$V5, strand=row$V6)
  sec_bed <- data.table(chrom=row$V11, start=row$V12, end=row$V13, name=row$V14, score=row$V15, strand=row$V16)
  
  bed.filtered.plus.minus <- rbindlist(list(bed.filtered.plus.minus, first_bed))
  bed.filtered.plus.minus <- rbindlist(list(bed.filtered.plus.minus, sec_bed))
  
  #bed.filtered.plus.minus[nrow(bed.filtered.plus.minus) + 1, ] <- first_bed
 #bed.filtered.plus.minus[nrow(bed.filtered.plus.minus) + 1, ] <- sec_bed
}

# bed file with coordinates as start and end of peak (only 1 bp long)
write.table(bed.filtered.plus.minus, file="filtered.plus.minus.peaks", quote=FALSE, 
            row.names=FALSE, col.names=FALSE, sep="\t")

bed.format.plus.minus.filtered <- apply(bed.filtered.plus.minus[1:100,], 1, function(i) {
  row <- i
  vals <- unlist(strsplit(row[["name"]], "_"))
  
  list(row[["chrom"]], vals[2], vals[3], row[["name"]],row[["score"]], row[["strand"]])
})

final.bed.plus.minus.filtered <- rbindlist(bed.format.plus.minus.filtered)

# bed file with the coordinates as the start and end of the TNE 
write.table(final.bed.plus.minus.filtered, file="filtered.plus.minus.bed", quote=FALSE, 
            row.names=FALSE, col.names=FALSE, sep="\t")

#compare the similarity of values between plus.minus and minus.plus
# compile a full list of bidirectional pairs
init_sample_num <- ( nrow(minus.plus) + nrow(plus.minus) ) * 0.5
print(init_sample_num)

total_bi_pairs <- merge(filtered.plus.minus, filtered.minus.plus, by.x="V4", by.y="V14")
nrow(total_bi_pairs) #14172
ncol(total_bi_pairs) #41

nrow(total_bi_pairs) / init_sample_num # 0.1276383

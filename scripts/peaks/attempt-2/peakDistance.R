library(data.table)
library(ggplot2)
library(dplyr)
setwd("~/Documents/college/dong_lab/code/playground/peaks/attempt-2")

# with the minus plus orientation 
minus.plus <- fread("minus.plus.bed", sep = "\t")

minus.plus.dist <- ggplot(minus.plus, aes(x=V21)) + 
  geom_histogram(color="black", fill="white", binwidth = 500) + 
  scale_y_continuous(trans="log2")

#TODO create a filter function to select for V21 values which are from -1000 to 1000 
minus.plus.density <- density(minus.plus$V21)
#returns the index 
y_max <- max(minus.plus.density$y)
index <- which(density(minus.plus$V21)$y == y_max)

ggplot(minus.plus, aes(x=V21, y=..density..)) + 
  geom_histogram(color="black", fill="white", binwidth = 30) + 
  scale_x_continuous(limits=c(-1000,1000)) + 
  geom_density() + 
  xlab("Distance to Nearest Peak") + 
  ggtitle("Minus-Plus Orientation Peak Distance") 
  #geom_vline(xintercept = minus.plus.density$x[index])


min_dis <- min(minus.plus$V21) # -2453205
max_dis <- max(minus.plus$V21) # 2230126

# no distances are na 
any(is.na(minus.plus$V21))

# only neg distances
neg_list <- minus.plus$V21 < 0

neg.minus.plus <- minus.plus[minus.plus$V21 < 0]
nrow(neg.minus.plus) # 34465 
length(minus.plus$V21) # 111440

print(paste0("percentage of TNEs where - is before + ", nrow(neg.minus.plus) / length(minus.plus$V21) ))
# ~ 30%

minus.plus.dist <- ggplot() + geom_histogram(aes(x=neg.minus.plus), color="black", fill="white") + 
  scale_x_continuous(limits=c(-1000,10)) # limit minimum x value to see distribution of values at higher res

# with the plus minus orientation 

plus.minus <- fread("plus.minus.bed", sep = "\t")

plus.minus.dist <- ggplot(plus.minus, aes(x=V21)) + 
  geom_histogram(color="black", fill="white") + 
  scale_y_continuous(trans="log2")

ggplot(plus.minus, aes(x=V21, y=..density..)) + 
  geom_histogram(color="black", fill="white", binwidth = 30) + 
  scale_x_continuous(limits=c(-1000,1000)) + 
  geom_density() + 
  xlab("Distance to Nearest Peak") + ggtitle("Plus-Minus Orientation Peak Distance")


plus.minus.dist

min_dis <- min(plus.minus$V21) # -1912881
max_dis <- max(plus.minus$V21) # 2815664

# no distances are na 
any(is.na(plus.minus$V21))

# only neg distances

neg.plus.minus <- plus.minus[plus.minus$V21 < 0, ]
print(nrow(neg.plus.minus)) # 34073
length(plus.minus$V21) # 110671

print(paste0("percentage of TNEs where - is before + ", nrow(neg.plus.minus) / length(plus.minus$V21) ))
# ~ 30%

minus.plus.dist <- ggplot() + geom_histogram(aes(x=neg.plus.minus), color="black", fill="white") + 
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
nrow(filtered.plus.minus) # 19,365

filtered.minus.plus  <- neg.minus.plus[neg.minus.plus$V21 > -400]
nrow(filtered.minus.plus) # 19,355


# there are 19,365 which are in reverse forward orientation with a distance of < 400 bp apart

###########reformat narrowPeak format into bed format 
narrowPeakToBed <- function(df, filename) {
  bed.filtered.df <- data.table(chrom=character(), 
                                        start=numeric(), 
                                        end=numeric(), 
                                        name=character(), 
                                        score=numeric(), 
                                        strand=character())
  
  for(i in 1:nrow(df)) {
    row <- df[i, ]
    first_bed <- data.table(chrom=row$V1, start=row$V2, end=row$V3, name=row$V4, score=row$V5, strand=row$V6)
    sec_bed <- data.table(chrom=row$V11, start=row$V12, end=row$V13, name=row$V14, score=row$V15, strand=row$V16)
    
    bed.filtered.df <- rbindlist(list(bed.filtered.df, first_bed))
    bed.filtered.df <- rbindlist(list(bed.filtered.df, sec_bed))
  }
  
  # this step is to get the original TNE coordinates 
  bed.format.df.filtered <- apply(bed.filtered.df[1:100,], 1, function(i) {
    row <- i
    vals <- unlist(strsplit(row[["name"]], "_"))
    
    list(row[["chrom"]], vals[2], vals[3], row[["name"]],row[["score"]], row[["strand"]])
  })
  
  
  final.bed.df.filtered <- rbindlist(bed.format.df.filtered)
  
  write.table(final.bed.df.filtered, file=filename, quote=FALSE, 
              row.names=FALSE, col.names=FALSE, sep="\t")
  
}


############compare the similarity of values between plus.minus and minus.plus
# compile a full list of filtered bidirectional pairs

total_sample_num <- nrow(filtered.minus.plus) + nrow(filtered.plus.minus)
print(total_sample_num)

setcolorder(filtered.plus.minus, 
            c(paste0("V", 11:20), paste0("V", 1:10), "V21"))

setnames(filtered.plus.minus, colnames(filtered.plus.minus), paste0("V", 1:21))


# All identified + and - pairs ! 
total_bi_pairs <- merge(filtered.plus.minus, filtered.minus.plus, all=T, sort=FALSE)
total_bi_pairs_test <- merge(filtered.plus.minus, filtered.minus.plus, sort=FALSE)
nrow(total_bi_pairs) #24740
ncol(total_bi_pairs)

nrow(total_bi_pairs) / init_sample_num # 0.2227715

num_shared_pairs <- total_sample_num - nrow(total_bi_pairs) #13980

#average of total number of filtered plus minus and minus plus pairs 
avg_bi_pair_filtered_num <- (nrow(filtered.plus.minus) + nrow(filtered.minus.plus)) * 0.5  
num_shared_pairs / avg_bi_pair_filtered_num # 0.7221074


# okay so out of 100,000 + TNE candidate regions, 22% are bidirectionally supported...
# does that make sense?
# I feel like the number should be higher like 40 or 30 loll 

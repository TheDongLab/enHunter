library(data.table)
library(dplyr)

setwd("~/Documents/college/dong_lab/code/enHunter/")

###### with the minus plus orientation 
minus.plus <- fread("input_files/closest/minus.plus.bed", sep = "\t")

####### with the plus minus orientation 
plus.minus <- fread("input_files/closest/plus.minus.bed", sep = "\t")


####### find the mode of the distances 
# what is the most common distance between peaks? 

# Create the function.
getmode <- function(v) {
  uniqv <- unique(v)
  #match returns a vector of position of (first) matches of its first argument in its second 
  #tabulate takes the integer-valued vector bin and counts the number of times each integer occurs in it 
  #which.max returns the index of the element with the maximal value in a vector 
  
  list(num = uniqv[which.max(tabulate(match(v, uniqv)))],
       rep_times = max(tabulate(match(v, uniqv))) )
}

plus.minus.mode <- getmode(plus.minus$V21)
print(plus.minus.mode)
# plus.minus.mode$rep_times = 1373

minus.plus.mode <- getmode(minus.plus$V21)
print(minus.plus.mode)
# minus.plus.mode$rep_times = 859

# 0 is the most common distance between peaks ! 
# what is the percentage of peak pairs with this distance? 
nrow(minus.plus) #42818
# 82 times 
nrow(plus.minus) #66785
# 82 times 

print(paste0("percentage of peak pairs with 0 distance between peaks for minus plus pairs: ", minus.plus.mode$rep_times/nrow(minus.plus)))
#0.00191508244196366
print(paste0("percentage of peak pairs with 0 distance between peaks for plus minus pairs: ", plus.minus.mode$rep_times/nrow(plus.minus)))
#0.00122782061840234

# only 0.1% of all peaks are the same between pairs  

####### filter for peak distances with 0 
zero_dist_plus_minus <- plus.minus[plus.minus$V21==0]
# nrow = 82

zero_dist_minus_plus <- minus.plus[minus.plus$V21==0]
# nrow = 82

write.table(zero_dist_plus_minus, file="zero_dist_plus_minus.tsv", sep="\t", col.names=FALSE, row.names=F, quote=F)
write.table(zero_dist_minus_plus, file="zero_dist_minus_plus.tsv", sep="\t", col.names=FALSE, row.names=F, quote=F)



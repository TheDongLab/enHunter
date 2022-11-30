library(data.table)
library(dplyr)
library(ggplot2)

wd <- "/Users/rosanwang/Documents/college/dong_lab/code/enHunter"
setwd(wd)

##### output from bedtools closest 
minus.plus <- fread("./input_files/closest/minus.plus.bed") %>% select(V1, V2, V3, V4, V11, V12, V13, V14, V21)
colnames(minus.plus) <- c("chr-minus", "start-minus", "end-minus", "minus", "chr-plus", "start-plus", "end-plus","plus", "distance")

plus.minus <- fread("./input_files/closest/plus.minus.bed") %>% select(V1, V2, V3, V4, V11, V12, V13, V14, V21)
colnames(plus.minus) <- c("chr-plus", "start-plus", "end-plus", "plus", "chr-minus", "start-minus", "end-minus", "minus", "distance")

##### TNE id and class information 
eRNA_minus <- fread("./input_files/characterization/eRNA.minus.characterize.xls") %>% select(V1, class)
eRNA_plus <- fread("./input_files/characterization/eRNA.plus.characterize.xls") %>% select(V1, class)

##### merging TNE id and class information 
class_minus.plus <- merge(minus.plus, eRNA_minus, by.x = "minus", by.y = "V1") %>% 
  rename(minus.class = class ) %>% merge(. , eRNA_plus, by.x = "plus", by.y = "V1") %>% 
  rename(plus.class = class ) %>% mutate(class_pairing = paste(minus.class, plus.class, sep = "_"))

class_plus.minus <- merge(plus.minus, eRNA_plus, by.x = "plus", by.y = "V1") %>% 
  rename(plus.class = class ) %>% merge(. , eRNA_minus, by.x = "minus", by.y = "V1") %>% 
  rename(minus.class = class ) %>% mutate(class_pairing = paste(plus.class, minus.class, sep = "_"))

##### making the bi directional pairing distance graphs 

minus_plus_bidir_class <- ggplot(class_minus.plus , aes(x=distance, fill=class_pairing)) + 
  geom_histogram(bins = 60) + facet_wrap(~class_pairing, scales = "free_y") +
  scale_x_continuous(limits=c(-1000,1000), n.breaks = 10) + xlab("Distance to Nearest Peak") + 
  theme(axis.line = element_line(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        text = element_text(size = 20),
        panel.border = element_blank()) +
  ggtitle("Minus-Plus Orientation Peak Distance") 

minus_plus_bidir_class

ggsave("minus_plus_bidir_class.pdf",
       plot=minus_plus_bidir_class, device="pdf")
ggsave("minus_plus_bidir_class.png", plot= minus_plus_bidir_class)

plus_minus_bidir_class <- ggplot(class_plus.minus , aes(x=distance, fill=class_pairing)) + 
  geom_histogram(bins = 60) + facet_wrap(~class_pairing, scales = "free_y") +
  scale_x_continuous(limits=c(-1000,1000), 
                     n.breaks = 10) + xlab("Distance to Nearest Peak") + 
  theme(axis.line = element_line(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        text = element_text(size = 20),
        panel.border = element_blank()) +
  ggtitle("Plus-Minus Orientation Peak Distance") 

plus_minus_bidir_class

ggsave("plus_minus_bidir_class.pdf",
       plot=plus_minus_bidir_class, device="pdf")
ggsave("plus_minus_bidir_class.png", plot= plus_minus_bidir_class)

# Warning messages:
# 1: Removed 14442 rows containing non-finite values (stat_bin). 
# this happens when the x axis is limited and cuts off values 

### combining both plus_minus and minus_plus class graphs 
# first combine the data frame 
temp_class_plus.minus <- class_plus.minus %>% select(plus, minus, distance, plus.class, minus.class)
temp_class_minus.plus <- class_minus.plus %>% select(minus, plus, distance, minus.class, plus.class)

class_combined_bidirectional <- full_join(temp_class_plus.minus, temp_class_minus.plus, 
                                          by =c("plus", "minus", "distance", "minus.class", "plus.class")) %>% 
  mutate( class = apply(., 1, function(x) {
    minus.class = x["minus.class"]
    plus.class = x["plus.class"]
    if( (minus.class == "1" & plus.class == "3") ||  (minus.class == "3" & plus.class == "1") ) {
      "3 and 1"
    } else if ((minus.class == "1" & plus.class == "2") |  (minus.class == "2" & plus.class == "1")) {
      "1 and 2"
    } else if ((minus.class == "3" & plus.class == "2") |  (minus.class == "2" & plus.class == "3")) {
      "2 and 3"
    } else if (minus.class == plus.class) {
      paste0(minus.class, " and ", plus.class)
    } else {
      "Failed"
    }
    
  }))

combined_bidir_class <- ggplot(class_combined_bidirectional, aes(x=distance, fill=class)) + 
  geom_histogram(bins = 60) + facet_wrap(~class, scales = "free_y") +
  scale_x_continuous(limits=c(-1000,1000), 
                     n.breaks = 10) + xlab("Distance to Nearest Peak") + 
  theme(axis.line = element_line(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        text = element_text(size = 20),
        panel.border = element_blank()) +
  ggtitle("Peak Distance") 

combined_bidir_class

ggsave("combined_bidir_class.pdf",
       plot=combined_bidir_class, device="pdf")

##### selecting only class 1 - class 1 pairs 
class1_plus.minus <- class_plus.minus[class_pairing == "1_1"] %>% mutate(strand1 = "+") %>% mutate(strand2 = "-")
class1_plus.minus_test <- class1_plus.minus[, c("strand1", "plus", "chr-plus", "start-plus", "end-plus", "plus.class",
                                           "strand2", "minus", "chr-minus", "start-minus", "end-minus", "minus.class", "distance")]

# plus                        minus                      distance  minus.class   plus.class  class_pairing 
# chr10_100524807_100524999   chr10_100524764_100524913    -40         1             1           1_1

class1_minus.plus <- class_minus.plus[class_pairing == "1_1"] %>% mutate(strand1 = "-") %>% mutate(strand2 = "+")
class1_minus.plus_test <- class1_minus.plus[, c("strand1", "minus", "chr-minus", "start-minus", "end-minus", "minus.class",
                                                "strand2", "plus", "chr-plus", "start-plus", "end-plus", "plus.class","distance")]

write.table(class1_plus.minus_test, file="class1_plus_minus_test.txt", sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(class1_minus.plus_test, file="class1_minus_plus_test.txt", sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

ggplot(class1_plus.minus, aes(x=distance)) + 
  geom_histogram(bins = 200) + 
  scale_x_continuous(limits=c(-1000,1000)) + 
  xlab("Distance to Nearest Peak") + 
  ggtitle("Plus-Minus Orientation Peak Distance") 

ggplot(class1_minus.plus, aes(x=distance)) + 
  geom_histogram(bins = 200) + 
  scale_x_continuous(limits=c(-1000,1000)) + xlab("Distance to Nearest Peak") + 
  ggtitle("Minus-Plus Orientation Peak Distance") 


### checking class1_list.sh in R 
#26326 + 26170 = 52496 number of non unique items 
all_minus <- unique(c(unlist(class1_minus.plus_test %>% select(minus)) , unlist(class1_plus.minus_test %>% select(minus)))) 
# number of unique items = 27266 (matches the shell script)

all_plus <- unique(c(unlist(class1_minus.plus_test %>% select(plus)) , unlist(class1_plus.minus_test %>% select(plus)))) 
# number of unique items = 27403 (matches shell script )

##### out of the birdirectional pairs with a peak distance of 0, what are their class break downs?
plus.minus.class.zero <- class_plus.minus[distance == 0]
minus.plus.class.zero <- class_minus.plus[distance == 0]

class_zero <- ggplot(plus.minus.class.zero, aes(x=class_pairing)) +
  geom_bar() + geom_text(aes(label=..count..), stat="count", vjust=0, colour="red") + 
  xlab("Class Pairing") + theme(axis.line = element_line(), panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(), panel.background = element_blank()) +
  ggtitle("Class Pairings of Bidirectional Pairs with 0 bp Distance") 

class_zero

ggsave("class_zero.png", plot=class_zero)

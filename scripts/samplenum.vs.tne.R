library(ggplot2)
library(dplyr)
library(ggtext)
library(ggrepel)

setwd("~/Documents/college/dong_lab/code/playground/peaks/attempt-2")

## number of TNE graph 
TNE_num <- read.csv("sample.num.tnes.csv", header = TRUE, sep = ",")

TNE_num_graph <- ggplot(TNE_num, aes(x= SampleNum, y=TNENum)) + geom_point() +
  expand_limits(y=0)


## similarity ratio graph 
similarity_ratio <- read.csv("similarity.ratios.csv", header = TRUE, sep = ",")
similarity_ratio$index <- 1:nrow(similarity_ratio)

#
similarity_graph <- ggplot(similarity_ratio, aes(x=factor(Ratio, Ratio), y=Similarity)) + geom_point() +
  expand_limits(y=0) + xlab("Sample Sizes Compared \n Number of TNEs Determined") + ylab("Similarity Ratio") + 
  geom_smooth(aes(x=index, y=Similarity), method = "loess", inherit.aes = FALSE) + 
  ggtitle("TNE Region Similarity for Each Sample Size") + theme(plot.title = element_text(hjust = 0.5))

##------------------------------------------------------------------------------
## THIS IS THE CODE YOU ACTUALLY USED 
similarity_ratio$label <- lapply(similarity_ratio$Ratio, function(x) {
  vals <- unlist(strsplit(x, " | ", fixed=TRUE))
  first_tne <- filter(TNE_num, SampleNum == vals[[1]]) %>% select(TNENum)
  first <- first_tne[1,1]
  
  sec_tne <- filter(TNE_num, SampleNum == vals[[2]]) %>% select(TNENum)
  sec <- sec_tne[1,1]
  print(sec)
  
  label <- paste0(vals[[1]], ":", first," TNEs\n",
                  vals[[2]], ":", sec, " TNEs")
  label
})


#,400,500,600,800,1000
similarity_ratio <- similarity_ratio[1:8,]
similarity_graph_adjusted_axis <- ggplot(
  similarity_ratio, 
  aes(x=c(2,10,20,40,60,80,100,200), y=Similarity, 
      label = label
      )) + 
  geom_point(color = "red") +
  expand_limits(y= c(0, 1)) + 
  xlab("Sample Size Ratio") + 
  ylab("Similarity Ratio") + 
  stat_smooth( method = "loess", se = FALSE, span = 1.3) + 
  scale_x_continuous(breaks=c(2,10,20,40,60,80,100,200), labels=similarity_ratio$Ratio) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(face = "bold"), 
        axis.text.y = element_text(face = "bold")
        ) + 
 geom_label_repel(nudge_y = 0.1, size = 2.75, fill = "white",  box.padding = 0, 
                   max.overlaps = Inf)
#angle=30,
#vjust = -0.005,
# ggtitle("TNE Region Similarity for Each Sample Size")


#subsampled similarity graph 
subsampled_similarity_ratio <- read.csv("subset_similarity.csv", header = TRUE, sep = ",")
subsampled_similarity_ratio$index <- 1:nrow(subsampled_similarity_ratio)

subsampled_tnenum <- read.csv("subsampled.num.tnes.csv", header = TRUE, sep = ",")

subsampled_similarity_ratio$label <- lapply(subsampled_similarity_ratio$Ratio, function(x) {
  vals <- unlist(strsplit(x, " | ", fixed=TRUE))
  first_tne <- filter(subsampled_tnenum, SampleNum == vals[[1]]) %>% select(TNENum)
  first <- first_tne[1,1]
  
  sec_tne <- filter(subsampled_tnenum, SampleNum == vals[[2]]) %>% select(TNENum)
  sec <- sec_tne[1,1]
  print(sec)
  
  label <- paste0(vals[[1]], ":", first," TNEs\n",
                  vals[[2]], ":", sec, " TNEs")
  label
})


subsampled_similarity_graph <- ggplot(
  subsampled_similarity_ratio, 
  aes(x=factor(Ratio, Ratio), y=Similarity, label = label)) + 
  geom_point(color = "red") +
  stat_smooth(aes(x=index, y=Similarity), method = "loess", se = FALSE, span = 1.15) + 
  expand_limits(y= c(0, 1)) + xlab("Subsampling Ratio") + ylab("Similarity Ratio") + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),  
        axis.text.x = element_text(face = "bold"), 
        axis.text.y = element_text(face = "bold")) + 
  geom_label_repel(nudge_y = 0.1, size = 2.75, fill = "white",  box.padding = 0, 
                   max.overlaps = Inf)

#ggtitle("TNE Region Similarity for Each Subsampled Percentage") 


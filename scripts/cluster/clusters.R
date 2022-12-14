### graph to plot the number of clusters computed based on the cut off distance 
## numbers were generated using bedtools cluster 

library(ggplot2)

distance <- c(300, 1000, 3000, 5000, 6000, 
              7000, 8000, 9000, 12500, 25000, 50000, 100000, 200000)
clusters <- c(72915, 47369, 31550, 25342, 23315,
              21578, 20110, 18985, 16070, 11232, 7940, 5417, 3362)

df <- data.frame(distance, clusters)

ggplot(df, aes(x=distance, y=clusters)) + 
  geom_point() + theme_minimal()


## taking distance = 100,000 and cluster = 5,417



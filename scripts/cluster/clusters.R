library(ggplot2)

distance <- c(300, 1000, 3000, 5000, 6000, 
              7000, 8000, 9000, 12500, 50000, 100000)
clusters <- c(72915, 47369, 31550, 25342, 23315,
              21578, 20110, 18985, 16070, 7940, 5417)

df <- data.frame(distance, clusters)

ggplot(df, aes(x=distance, y=clusters)) + geom_point()

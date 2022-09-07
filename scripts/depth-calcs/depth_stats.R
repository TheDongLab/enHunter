library(ggplot2)

depths <- read.table("depths.txt")
depths$index <- 1:nrow(depths)
mean <- mean(depths[[1]])
median <- median(depths[[1]])
hist(depths[[1]])

ggplot(depths, aes(x=index, y=V1)) + geom_point()


full_depths <- read.table("full-sizeddepths.txt")
full_mean <- mean(full_depths[[1]])
full_median <- median(full_depths[[1]])
hist(full_depths[[1]])

filtered_depths <- read.table("filtered-depth.txt")
filtered_mean <- mean(filtered_depths[[1]])
filtered_median <- median(filtered_depths[[1]])
hist(filtered_depths[[1]])

filtered_depths_0x4 <- read.table("filtered_depths0x4.txt")
filtered_mean_0x4 <- mean(filtered_depths_0x4[[1]])
filtered_median_0x4 <- median(filtered_depths_0x4[[1]])
hist(filtered_depths_0x4[[1]])

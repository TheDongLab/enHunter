library(data.table)
library(ggplot2)
library(dplyr)

# bidirectional pairs < 900 bp apart 
bidir.pairs <- fread("input_files/closest/bidirectional_pairs_class.csv")

dist_density_minus.plus <- ggplot(bidir.pairs, aes(x=distance)) + 
  geom_histogram(color="black", fill="white", binwidth = 25) + 
  scale_x_continuous(n.breaks = 10) + 
  theme(axis.line = element_line(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        panel.border = element_blank(), text = element_text(size = 20))  +
  xlab("Distance to Nearest Peak") + 
  ggtitle("Bidirectional Pairs Peak Distance")

dist_density_minus.plus


class.1.or.2.pairs <- bidir.pairs %>% filter((class1 == 1 | class1 == 2) & (class2 == 1 | class2 == 2))
class.1.or.2.pairs.med <- median(class.1.or.2.pairs$distance, na.rm=TRUE)

class.1.or.2.pairs.bidir <- ggplot(class.1.or.2.pairs, aes(x=distance)) + 
  geom_histogram(color="black", fill="lightgray", binwidth = 25) + 
  scale_x_continuous(n.breaks = 10) + 
  scale_y_continuous(n.breaks = 10) +
  theme(axis.line = element_line(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        panel.border = element_blank(), text = element_text(size = 20))  +
  xlab("Distance to Nearest Peak") + 
  ggtitle("Bidirectional Pairs Peak Distance") +
  geom_vline(xintercept=class.1.or.2.pairs.med, color="red", linetype="dashed", linewidth = 1)

class.1.or.2.pairs.bidir

ggsave( "class.1.or.2.bidir.dist.900bp.pdf",plot = class.1.or.2.pairs.bidir, device = "pdf")




################ all bidirectional pairs ##########
all.pairs <- fread("all_bidirectional_pairs_class.csv")
all.class.1.or.2.pairs <- all.pairs  %>% filter((class1 == 1 | class1 == 2) & (class2 == 1 | class2 == 2))

all.class.1.or.2.pairs.bidir <- ggplot(all.class.1.or.2.pairs, aes(x=distance)) + 
  geom_histogram(color="black", fill="lightgray", binwidth = 20) + 
  scale_x_continuous(n.breaks = 10, limits = c(-1000, 1000)) + 
  scale_y_continuous(n.breaks = 10) +
  theme(axis.line = element_line(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        panel.border = element_blank(), text = element_text(size = 20))  +
  xlab("Distance to Nearest Peak") + 
  ggtitle("Bidirectional Pairs Peak Distance") +
  geom_vline(xintercept=0, color="red", linetype="dashed", linewidth = 1) +
  geom_vline(xintercept=-900, color="red", linetype="dashed", linewidth = 1)

all.class.1.or.2.pairs.bidir

ggsave("all.class.1.or.2.pairs.bidir.pdf", plot = all.class.1.or.2.pairs.bidir, device = "pdf")

all.class.1.or.2.pairs[-900<distance & distance < 0]


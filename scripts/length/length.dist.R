library(data.table)
library(ggplot2)

wd <- "~/Documents/college/dong_lab/code/playground"
setwd(wd)

### calculating the length of each TNE with awk 
system("awk '{OFS=\"\t\"; print $4,$3-$2}' plus.bed > plus.length.bed")

### generating distribution of TNE lengths 
plus <- fread("plus.length.bed", sep = "\t")

plus.mean <- mean(plus$V2, na.rm=TRUE)

plusLength <- ggplot(plus, aes(x=V2)) + 
  geom_histogram(fill="gray", binwidth=0.02)  + 
  scale_x_continuous(trans="log10") + theme_bw() + 
  geom_vline(xintercept=plus.mean, color="red", linetype="dashed") + 
  theme(axis.line = element_line(), panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(), 
        panel.border = element_blank()) + xlab("TNE Length") + ylab("Count") + 
  annotate("text", x=plus.mean+240, y=8000, 
           label=paste0("mean: ", round(plus.mean, digits = 2), " bp"), color="red")

plusLength
ggsave("./length/plusLength.pdf", plot=plusLength, device="pdf")
ggsave("./length/plusLength.png", plot=plusLength, device="png")

system("awk '{OFS=\"\t\"; print $4,$3-$2}' minus.bed > minus.length.bed")
minus <- fread("minus.length.bed", sep = "\t")

minus.mean <- mean(minus$V2, na.rm=TRUE)

minusLength <- ggplot(minus, aes(x=V2)) + 
  geom_histogram(fill="gray", binwidth=0.02)  + 
  scale_x_continuous(trans="log10") + theme_bw() + 
  geom_vline(xintercept=minus.mean, color="red", linetype="dashed") +
  theme(axis.line = element_line(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_blank()) + xlab("TNE Length") + ylab("Count") +
  annotate("text", x =minus.mean+240, y=8000, 
           label=paste0("mean: ", round(minus.mean, digits = 2), " bp"), color="red")

minusLength
ggsave("./length/minusLength.pdf", plot=minusLength, device="pdf")
ggsave("./length/minusLength.png", plot=minusLength, device="png")

#testing to see how the distribution looks when not in log 
ggplot(minus, aes(x=V2)) + 
  geom_histogram(fill="gray", binwidth=10)  + 
  scale_x_continuous(limits=c(0, 1000)) + theme_bw() + 
  geom_vline(xintercept=minus.mean, color="red")

### generating distribution of enhancer lengths from fantom5
system(paste0("awk '{OFS=\"\t\"; print $4,$3-$2}' ", wd, 
"/length/permissive_enhancers_hg38_blood_only.bed > ", wd, "/length/permissive_enhancers_hg38_blood_only_length.bed"))

fantom5 <- fread("./length/permissive_enhancers_hg38_blood_only_length.bed")

fantom5.mean <- mean(fantom5$V2, na.rm=TRUE)

fantom5Length <- ggplot(fantom5 , aes(x=V2)) + 
  geom_histogram(fill="gray", binwidth=0.02)  + 
  scale_x_continuous(trans="log10") + theme_bw() + 
  geom_vline(xintercept=fantom5.mean, color="red", linetype="dashed") + 
  theme(axis.line = element_line(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_blank()) + xlab("TNE Length")+ ylab("Count") +
  annotate("text", x=fantom5.mean+300, y=2500, 
           label=paste0("mean: ", round(fantom5.mean, digits = 2), " bp"), color="red")

fantom5Length

ggsave("./length/fantom5Length.pdf", plot=fantom5Length, device="pdf")
ggsave("./length/fantom5Length.png", plot=fantom5Length, device="png")

#testing to see how the distribution looks when not in log 
ggplot(fantom5, aes(x=V2)) + 
  geom_histogram(fill="gray", binwidth=10)  + 
  scale_x_continuous(limits=c(0, 1000)) + theme_bw() + 
  geom_vline(xintercept=fantom5.mean, color="red")

### generating distribution of enhancer lengths from Roadmap
system(paste0("awk '{OFS=\"\t\"; print $4,$3-$2}' ", wd, 
              "/length/15_coreMarks_hg38lift_segments.E6E7E12.bed > ", wd, "/length/15_coreMarks_hg38lift_segments.E6E7E12_length.bed"))

roadmap <- fread("./length/15_coreMarks_hg38lift_segments.E6E7E12_length.bed")
roadmap.mean <- mean(roadmap$V2, na.rm=TRUE)

roadmapLength <- ggplot(roadmap , aes(x=V2)) + 
  geom_histogram(fill="gray", binwidth=0.08)  + 
  scale_x_continuous(trans="log10") + theme_bw() + 
  geom_vline(xintercept=roadmap.mean, color="red", linetype="dashed") + 
  theme(axis.line = element_line(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_blank()) + xlab("TNE Length")+ ylab("Count") +
  annotate("text", x=roadmap.mean+2000, y=400000, 
           label=paste0("mean: ", round(roadmap.mean, digits = 2), " bp"), color="red")

roadmapLength

ggsave("./length/roadmapLength.pdf", plot=roadmapLength, device="pdf")
ggsave("./length/roadmapLength.png", plot=roadmapLength, device="png")


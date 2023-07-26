#########################################
# referenced from: https://github.com/sterding/BRAINcode/blob/7b4c5e816ff1cf9af86041326b71cf3f3e2e4bf6/src/eRNA.classification.R#L58 
#########################################

library(pheatmap)

features = read.table("./input_files/eRNA.characterize.feature.color.stranded.xls")
df=subset(features, select=c(f06.TFBS, f07.P300, f08.CAGEbloodenhancer, f09.chromHMM_blood, f12.DNaseROADMAP, f15.HCNE))

df$f06.TFBS=ifelse(df$f06.TFBS>=5,1,0)

df[df>0]=1 #if the value is greater than 0, then convert to 1

library(plyr)
df0=count(df, vars = colnames(df)) #get frequencies
mat_data=df0[,1:6]

# clustering 
library(ade4)
library(cluster)

# matrix of distance between each of the rows
row_distance = daisy(mat_data, metric = "gower", type=list(asymm=c(1:6)), weights = c(1,1,1,1,5,1))
# cluster rows based on their distance
row_cluster = hclust(row_distance, method = "single") # single linkage method 
plot(row_cluster)
dev.off()

# cor(mat_data) = correlation of the different features (TFBS, P300, ...)
col_distance = as.dist((1 - cor(mat_data))/2)
col_cluster = hclust(col_distance, method = "complete") # complete linkage or farthest neighbor

gr.row <- cutree(row_cluster, 7) # form 6 groups (6 features)
gr.col <- cutree(col_cluster, 4) # form 3 groups (3 classes)

# changed to 43, as that is the number of rows of the number of unique frequencies 
row_cluster <- rev(reorder(as.dendrogram(row_cluster), 43:1, min))  
col_cluster <- rev(reorder(as.dendrogram(col_cluster), 1:6))


library(gplots)
library(RColorBrewer)

col1=c("blue","purple","red")
# specify how many rows are class 1, 2, or 3
gr.row[order.dendrogram(row_cluster)]=c(1,rep(2,22), rep(3,20))

col2 <- brewer.pal(4, "Pastel1")
gr.col[order.dendrogram(col_cluster)]=c(1,1,2,2,3,4)

pdf("eRNA.clustering.class.pdf", width=4, height=7)
heatmap.2(as.matrix(mat_data),
          margins=c(9,9), keysize=0.1, sepcolor=NA,sepwidth=c(0,0),
          lwid=c(1.5,0.25,4), lhei=c(0.1,0.01,0.5),
          #lmat=rbind(c(5,0,4),c(3,1,2)),
          lmat=rbind(c(6,0,5),c(0,0,2),c(4,1,3)),
          #cellnote = mat_data,  # same data set for cell labels
          main = NA, # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          col=colorRampPalette(c("gray", "black"))(n = 2),       # use on color palette defined earlier
          cexCol=1,
          cexRow=1,
          labRow = df0$freq, #[rownames(mat_data)],
          labCol = sub("f[0-9]+.(.*)","\\1",colnames(mat_data)),
          #breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="both",     # only draw a row dendrogram
          Colv = as.dendrogram(col_cluster),
          Rowv = as.dendrogram(row_cluster),
          RowSideColors=col1[gr.row],
          ColSideColors=col2[gr.col],
          key.title=NA, key.ylab=NA, key.xlab = NA, key.par=list(mgp=c(3,0.3,0))
)

dev.off()

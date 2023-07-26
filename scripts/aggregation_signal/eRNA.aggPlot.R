## script to analysis enhancer regions
## taken from https://github.com/sterding/BRAINcode/blob/master/src/eRNA.aggPlot.R 

setwd("~/eRNAseq/HCILB_SNDA")

#============================================================
# draw aggregation plot (to run in R console)
#============================================================

# the following inputs for this script are taken from amppd.eRNA.aggregation.sh 
# H3K4me3.Roadmap.brainMerged.bigwig -> Roadmap/RoadmapBrainH3K4me3Merged.bigwig
# H3K4me1.Roadmap.brainMerged.bigwig -> Roadmap/RoadmapBrainH3K4me1Merged.bigwig
# H3K27ac.Roadmap.brainMerged.bigwig -> Roadmap/RoadmapBrainH3K27acMerged.bigwig
# DNase/merged.DNase.pval.signal

#from the get.regulatoryMarks.data.sh script 
#TFBS/TFBS.ENCODE.all.count
#CAGE/CAGE.FANTOM5.total.rev
#CAGE/CAGE.FANTOM5.total.fwd
#Conservation/Conservation/Conservation.phyloP100way

mark="eRNA"
marks=c("RNAseq/RNAseq.BRAINCODE.HCILB_SNDA",
        "DNase/merged.DNase.pval.signal",  # all Roadmap brain samples
        "Histone/H3k27ac.Roadmap.brainMerged", # all Roadmap brain samples
        "Histone/H3k4me1.Roadmap.brainMerged", # all Roadmap brain samples
        "Histone/H3k4me3.Roadmap.brainMerged", # all Roadmap brain samples
        "TFBS/TFBS.ENCODE.all.count",
        "CAGE/CAGE.FANTOM5.total.rev",
        "CAGE/CAGE.FANTOM5.total.fwd",
        "Conservation/Conservation.phyloP46way")

## read the binned data
X=data.frame()
for(i in 1:length(marks))
{
  binfile = paste0("~/neurogen/external_download/externalData/", marks[i], ".bigwig.",mark,".1kbp.summit.or.mid.1kbp.100bins")
  message(paste("Loading:",binfile))
  
  x0=read.table(binfile, check.names =F, stringsAsFactors=F, header=F); 
  rownames(x0)=x0[,1]; x0=x0[,-1];
  if(ncol(X)>0) {X=cbind(X,x0);} else { X=x0;}
}

## read the classification of eRNA candidate
features = read.table("eRNA.characterize.xls", check.names =F, stringsAsFactors=F, header=T)
library('dplyr')
df=features %>% select(f12.DNaseROADMAP, f09.chromHMM_brain, f07.P300, f06.TFBS, f08.CAGEenhancer, f15.HCNE) 
df$f06.TFBS=ifelse(df$f06.TFBS>=5,1,0)
df[df>0]=1;
df$class=ifelse(apply(df,1,sum)==0,1,ifelse(df$f12.DNaseROADMAP==0,ifelse(apply(df,1,sum)==1,2,ifelse(apply(df,1,sum)==2,3,4)),5))
df$subclass=apply(df[,1:6],1,paste0,collapse="")

MARKS=c("RNAseq","DNase","H3K27ac","H3K4me1","H3K4me3","nTFBS","CAGE(-)","CAGE(+)","phyloP")
cols=c("red", 'darkgreen','darkorange4','darkorange4','blue','purple','magenta', 'darkmagenta','darkblue')

## to show aggregation of each class in separate rows
## -------------------------------------
# png("aggregation.eRNAs.class.all.png",width=1000, height=800)
# par(mfrow=c(6,8), mar=c(.5,0,0,0), oma=c(2,1,1,1))
# layout(matrix(seq(6*length(marks)),nrow=6,ncol=length(marks), byrow=1),widths=rep(1,length(marks)),heights=rep(1,6))
# ymin=c(0.03,0.68,0.68,0,0.2,0.2,0,0.05)
# ymax=c(0.1,1.45,1.45,8.5,2.2,2.3,10,0.3)
# for(k in c(0:5))
# {
#   if(k==0) index=rownames(X) else { index=rownames(subset(df,class==k))}
#   for(i in 1:length(marks))
#   {
#     x=X[index,(100*(i-1)+1):(i*100)]
#     plot(apply(x, 2, mean, trim =0.01), type='l',col=cols[i], ylim=c(ymin[i],ymax[i]),xlab="",ylab="",xaxt="n",las=2, cex.axis=1,tck=0.04, frame.plot=T, mgp=c(3,-2,0))
#     V=ncol(x)/2+.5
#     abline(v=V, lty=2, lwd=.5)
#   }
# }
# dev.off()

# figure 2a (Braincode paper)
message("
        ## to show aggregation of each class in one row  (Fig. 3d)
        ## -------------------------------------
        ")
pdf("aggregation.eRNAs.class.four.pdf",width=18, height=2, paper = 'usr')
par(mfrow=c(1,9), mar=c(.5,0,2,0), oma=c(2,1,1,1), adj=c(0,0))
layout(matrix(seq(1*length(marks)),nrow=1,ncol=length(marks), byrow=1),widths=rep(1,length(marks)),heights=1)
ymin=c(0.03, 18,  1, 1,  1, 0.0, 0.8, 0.8, 0.03)
ymax=c(0.12, 180,10,13,13, 2.5, 1.4, 1.4, 0.2)
for(i in 1:length(marks))
{
  # all
  x=X[,(100*(i-1)+1):(i*100)]
  plot(apply(x, 2, mean, trim =0.01), main=MARKS[i], type='l',col='black', ylim=c(ymin[i],ymax[i]),xlab="",ylab="",xaxt="n",las=2, cex.axis=1,tck=0.04, frame.plot=T, mgp=c(3,-2,0))
  
  # class III (no mark)
  x=X[rownames(subset(df, class==1)),(100*(i-1)+1):(i*100)]
  points(x=1:100, y=apply(x, 2, mean, trim =0.01), type='l',col='#ffcccc')
  
  # class II (no DNase, but at least one other mark)
  x=X[rownames(subset(df, class>1 & class<5)),(100*(i-1)+1):(i*100)]
  points(x=1:100, y=apply(x, 2, mean, trim =0.01), type='l',col='#ff9966')
  
  # class I (DNase)
  x=X[rownames(subset(df, class==5)),(100*(i-1)+1):(i*100)]
  points(x=1:100, y=apply(x, 2, mean, trim =0.01), type='l',col='#ff3333')
  
  V=ncol(x)/2+.5
  abline(v=V, lty=2, lwd=.5)
}

dev.off()



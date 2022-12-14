#!/usr/bin/env Rscript 
## Rscript to calculate the consistence of TNE expression across the samples
## Usage: Rscript TNE_caller.consistency.R $SAMPLE_GROUP $sampleIDs $sampleRDBGs $sampleMeanRPM
## Author: Xianjun Dong / Rosan Wang 
## Date: Nov-16-2022

args<-commandArgs(TRUE)

SAMPLE_GROUP = args[1]
SampleIDs = args[2]
SampleRDBGs = args[3]
SamplemeanRPM = args[4]

#write the meanRPM and rdbg file paths to a file which I then read in 
#TODO: write the files of 
#1. sample IDs 
#2. meanRPM files 
#3. rdbg files 

sampleID = read.table(SampleIDs, header = F, col.names = c('ID'), stringsAsFactors = F)
rdbgs = read.table(SampleRDBGs, header = F, col.names = c('RDBG'), stringsAsFactors = F)
meanRPM = read.table(SamplemeanRPM, header = F, col.names = c('RPM'), stringsAsFactors = F)

EXP=data.frame(); PV=data.frame(); 

for(i in 1:nrow(sampleID)){
  cur_sampleName <- sampleID$ID[[i]]

  print(cur_sampleName)

  cur_meanRPM <- meanRPM$RPM[[i]]
  cur_rdbg <- rdbgs$RDBG[[i]]

  # read background
  df=read.table(cur_rdbg, header=F)[,2] # mean RPM (mean0 from bigWigAverageOverBed)
  Fn=ecdf(df)
  
  # read expression
  expression=read.table(cur_meanRPM, header=F)
  pvalue=as.numeric(format(1-Fn(expression[,2]), digits=3));

  # merge
  if(ncol(EXP)==0) { EXP=expression; expression[,2]=pvalue; PV=expression; }
  else {EXP=cbind(EXP, expression[,2]); PV=cbind(PV, pvalue); }
}

colnames(EXP)=c("locus", sampleID$ID); colnames(PV)=c("locus", sampleID$ID); 
write.table(EXP, "eRNA.tmp5.meanRPM.xls", col.names=T, row.names=F, sep="\t", quote=F)
write.table(PV,  "eRNA.tmp5.pvalues.xls", col.names=T, row.names=F, sep="\t", quote=F)

## binomial test for the significant HTNE (p<0.05)
N=nrow(sampleID)
binomial.pvalues = sapply(rowSums(PV[,-1]<=0.05), function(x) binom.test(x,N,0.05,'greater')$p.value)
# using the Holm-Bonferroni method (AKA step-down Bonferroni)  to correct for multiple test
p.adjusted = cbind(binomial.pvalues=binomial.pvalues, p.adjusted.HB = p.adjust(binomial.pvalues, method = "holm"), p.adjusted.bonferroni=p.adjust(binomial.pvalues, method = "bonferroni"),p.adjusted.FDR=p.adjust(binomial.pvalues, method = "fdr"))
rownames(p.adjusted) = PV[,1]

write.table(p.adjusted,  "eRNA.tmp5.pvalues.adjusted.xls", col.names=NA, row.names=T, sep="\t", quote=F)
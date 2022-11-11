# ===============================
# Usage: Rscript ./colorBedByFeatures.R eRNA.characterization.xls minus 
# ===============================

# given a list of TNEs and the characterization table, colors the TNEs based on the number of features supporting them 
# colors TNEs based on the following features: f06.TFBS, f07.P300, f08.CAGEenhancer, f09.chromHMM_brain, f12.DNaseROADMAP. f15.HCNE

#args <- commandArgs(TRUE)

#features <- read.table(args[1])
features <- read.table("./input_files/characterization/eRNA.minus.characterize.xls")

df=subset(features, select=c(f06.TFBS, f07.P300, f08.CAGEbloodenhancer, f09.chromHMM_blood, f12.DNaseROADMAP, f15.HCNE))
df$f06.TFBS=ifelse(df$f06.TFBS>=5,1,0)
df[df>0]=1;
features$features=apply(df,1,sum)

#strand=args[2]
strand <- "minus"

find_rbg <- function(x) {
  
  class <- x[["features"]]
  
  rbg <- "0,0,0" 
  
  if(strand == "plus") {
    if(class== "0") {
      # supported by only 1 feature = light red
      rbg = "255,204,204"
    } else if (class== "1") {
      rbg = "255,153,153"
    } else if (class== "2") {
      rbg = "255,102,102"
    } else if (class=="3") {
      rbg ="255,51,51"
    } else if (class=="4") {
      rbg = "255,0,0"
    } else if(class=="5") {
      rbg = "204,0,0"
    } else if(class=="6") {
      rbg = "153,0,0"
    }
    
  }  else if(strand=="minus") {
    if(class=="0") {
      # supported by only 1 feature = light red
      rbg = "204,229,255"
    } else if (class=="1") {
      rbg = "153,204,255"
    } else if (class=="2") {
      rbg = "102,178,255"
    } else if (class=="3") {
      rbg = "51,153,255"
    } else if (class=="4") {
      rbg = "0,128,255"
    } else if(class=="5") {
      rbg = "0,102,204"
    } else if(class=="6") {
      rbg = "0,76,153"
    }
  }
  return(rbg)
}

features$rbg <- apply(features,1,find_rbg)
features$strand <- "-"

#write.table(features, paste0("eRNA.minus.characterize", ".xls"), sep="\t", quote = F, col.names = NA, row.names = T)

final_df <- subset(features, select=c(rbg))


# TO CREATE BED FILES FOR UCSC GENOME BROWSER VISUALIZATION
write.table(final_df, paste0("rbg_", strand, ".tsv"), sep="\t", 
            quote = F, col.names = F, row.names = T)

# sort -k4,4 ucsc.plus.narrowPeak.bed > ucsc.plus.sorted.narrowPeak.bed
# sort -k4,4 rbg_plus.tsv > sorted_rbg_plus.tsv
# join -1 4 ucsc.plus.sorted.narrowPeak.bed sorted_rbg_plus.tsv | awk 'OFS="\t" {print $2, $3, $4, $1, $5, $6, $7, $8, $9}' > test.bed

# join -1 4 ucsc.minus.sorted.narrowPeak.bed sorted_rbg_minus.tsv | awk 'OFS="\t" {print $2, $3, $4, $1, $5, $6, $7, $8, $9}' > test.bed

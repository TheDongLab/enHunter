#!/bin/bash

### input: cell type (Mon, Mac0, Mac1, Mac2... etc)

#### functionality ####
### selects for significant interactions for the given cell type (CHiCAGO score >= 5)
### puts significant interactions in interact file format 
### liftover coordinates from hg19 to hg38 
### creates a bigInteract file 

###### remember to activate conda before running ##########
#conda activate /PHShome/rw552/condaenvs/ucsc 

mkdir $1 
cd $1 

inputdir=/data/bioinformatics/external_data/externalData/Hi-C

awk -v t="Mon" 'BEGIN {OFS="\t"}
{ if ( NR == 1 ) {
  split($0, array); 
  col=0;
  for(i=12;i<=28;i++)if(array[i]==t){col=i};
  
 } else if ($col >= 5) {
  
   if (t == "Mon") color = "#F52D24";
    else if (t == "Mac0") color = "#BB312A";
    else if (t == "Mac1") color = "#F29809";
    else if (t == "Mac2") color = "#FFC761";
    else if (t == "Neu") color = "#EFCB09";
    else if (t == "MK") color = "#D0E009";
    else if (t == "EP") color = "#C5D84A"; 
    else if (t == "Ery") color = "#01B850"; 
    else if (t == "FoeT") color = "#01A689"; 
    else if (t == "nCD4") color = "#02A3A3"; 
    else if (t == "tCD4") color = "#038Ab0"; 
    else if (t == "aCD4") color = "#7AD2F7"; 
    else if (t == "naCD4") color = "#096EBE";
    else if (t == "nCD8") color = "#8D00DE";
    else if (t == "tCD8") color = "#CB96C4"; 
    else if (t == "nB") color = "#C40093"; 
    else if (t == "tB") color = "#E10064"; 
    else color = "ERROR"
    
    if ($1 != $6) {
      start = $2; 
      end = $3;
    } else {
      if($2 < $7) {
        start=$2;
        end=$8;
      } else {
        start=$7; 
        end=$3;
      }
    }
    
    print "chr"$1, start, end, "ID"NR, 0, $12, t, color, "chr"$1, $2, $3, $4, ".", "chr"$6, $7, $8, $9, "."
  }
}' $inputdir/PCHiC_peak_matrix_cutoff5.txt > $1.test.hg19.interact

## liftover the hg19 coordinates to hg38 
## this must be done in portions as the liftOver tool only handles bed files
module load liftover/1.0

cut -f 1,2,3,4 $1.test.hg19.interact | liftOver stdin $inputdir/hg19ToHg38.over.chain.gz $1.cor.hg38.bed unMappedPCHiC
awk 'BEGIN{OFS="\t"; FS="\t"} { print $9, $10, $11, $4, $13}' $1.test.hg19.interact | liftOver stdin $inputdir/hg19ToHg38.over.chain.gz $1.start.hg38.bed unMappedPCHiCstart
awk 'BEGIN{OFS="\t"} { print $14, $15, $16, $4, $18 }' $1.test.hg19.interact | liftOver stdin $inputdir/hg19ToHg38.over.chain.gz $1.end.hg38.bed unMappedPCHiCend
cut -f 5,6,7,8,4 $1.test.hg19.interact  > $1.color.tmp

### merge based on field 4 from above 
join -1 4 -2 1 <(LC_ALL=C sort -k 4,4b $1.cor.hg38.bed) <(LC_ALL=C sort -k 1,1b $1.color.tmp) | LC_ALL=C sort -k 1,1b | join -1 1 -2 4 - <(LC_ALL=C sort -k 4,4b $1.start.hg38.bed) | join -1 1 -2 4 - <(LC_ALL=C sort -k 4,4b $1.end.hg38.bed) | awk 'BEGIN{OFS="\t"} {print $2, $3, $4, $1, $5, $6, $7, $8, $9, $10, $11, $1, $12, $13, $14, $15, $1, $16}' | LC_ALL=C sort -k1,1 -k2,2n > $1.hg38.interact

# add the header
# sed  -i '1i track type=interact name=interact_PCHiC description=An_interact_file interactDirectional=true maxHeightPixels=200:100:50 visibility=full' $1.hg38.intersect


#curl https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes > hg38.chrom.sizes
dir=/PHShome/rw552/Documents

## to install bedToBigBed, openssl must be version 1.0 
# bedToBigBed -as=$dir/interact.as -type=bed5+13 mon.hg38.sorted.intersect $dir/hg38.chrom.sizes mon.hg38.inter.bb
bedToBigBed -as=$dir/interact.as -type=bed5+13 $1.hg38.interact $dir/hg38.chrom.sizes $1.hg38.inter.bb








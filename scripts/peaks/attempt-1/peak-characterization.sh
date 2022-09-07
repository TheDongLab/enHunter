#!/bin/bash

#peak files

head -n5 full-size/peaks.n200.minus.narrowPeak > n5peaks.n200.minus.narrowPeak
head -n5 full-size/peaks.n200.plus.narrowPeak > n5peaks.n200.plus.narrowPeak

head -n5 full-size/minus.n200.eRNA.bed > n5minus.n200.eRNA.bed
head -n5 full-size/plus.n200.eRNA.bed > n5plus.n200.eRNA.bed

minusPeak=n5peaks.n200.minus.narrowPeak
plusPeak=n5peaks.n200.plus.narrowPeak

#input bed file 
sampleBed=n5plus.n200.eRNA.bed

#bed file being tested against
compBed=n5minus.n200.eRNA.bed 

#should be either "minus-plus" or "plus-minus"
#basically dictates the strand direction of sample-comp
orientation=plus-minus

# the -D ref parameter doesn't do anything 
# I just have it there to get a gauge of what the distance between peaks is 
bedtools closest -a $sampleBed -b $compBed -D ref > $orientation.closest.bed

while read curChrom curStart curEnd cur closeChrom closeStart closeEnd close dist
do
  if [[ "$orientation" = "plus-minus" ]]
  then 
      curPeak=$plusPeak
      closePeak=$minusPeak
  else
      curPeak=$minusPeak
      closePeak=$plusPeak
  fi 

closeCor=`grep -w $close $closePeak | awk '{OFS="\t"; start=$10+$2} END{ print start}'`
curCor=`grep -w $cur $curPeak | awk '{OFS="\t"; start=$10+$2} END{ print start}'`

  if [[ "$orientation" == "plus-minus" ]]
  then 
    dist="$((closeCor - curCor))"
  else
    dist="$((curCor - closeCor))"
        #echo "$closeCor"
  fi 
  
  echo -e "${cur}\t${dist}" >> $orientation.bed
  
done < $orientation.closest.bed

#Rscript --vanilla peak_distances.R $orientation $orientation.closest.bed $minusPeak $plusPeak
# peak_distances.R "minus-plus" minus.plus.bed smaller.peaks.n200.minus.narrowPeak smaller.peaks.n200.plus.narrowPeak 

#output: $orientation.bed
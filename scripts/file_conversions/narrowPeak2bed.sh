#!/bin/bash
# usage: narrowPeak2bed.sh peaks.n200.minus.narrowPeak
# note: this takes in the NEW narrowPeak format
# where start and end coordinates represent the peak 

input=$1

awk '{
  OFS="\t";
  split($4,a,"_"); 
  start=a[2];
  end=a[3];
  print $1, start, end, $4, $5, $6, $2, $3
}' $input > ucsc.narrowPeak.bed
#!/bin/bash
# usage: closest2bed.sh minus.plus.bed

input=$1

mkdir tmp 

cut -f 1-10 $input > ./tmp/first.bed
cut -f 11-20 $input > ./tmp/sec.bed

awk '{
  OFS="\t";
  split($4,a,"_"); 
  start=a[2];
  end=a[3];
  print $1, start, end, $4, $5, $6, $2, $3
}' ./tmp/first.bed > ucsc.first.bed

awk '{
  OFS="\t";
  split($4,a,"_"); 
  start=a[2];
  end=a[3];
  print $1, start, end, $4, $5, $6, $2, $3
}' ./tmp/sec.bed > ucsc.sec.bed

cat ucsc.first.bed ucsc.sec.bed > closest.bed

rm ucsc.sec.bed
rm ucsc.first.bed

rm -r tmp 
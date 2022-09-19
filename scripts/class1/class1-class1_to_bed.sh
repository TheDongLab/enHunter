#!/bin/bash
# usage: ./class1-class1_to_bed.sh class1_minus_plus.txt
# usage: ./class1-class1_to_bed.sh class1_plus_minus.txt

# plus                 minus               chr-minus start-minus end-minus chr-plus start-plus end-plus
#chr10_100524807_100524999	chr10_100524764_100524913	chr10	100524839	100524840	chr10	100524879	100524880	-40	1	1_1	-	+
#chr10_100538460_100538587	chr10_100537628_100537735	chr10	100537672	100537673	chr10	100538509	100538510	-837	1	1_1	-	+
#chr10_100538460_100538587	chr10_100538495_100538650	chr10	100538592	100538593	chr10	100538509	100538510	83	1	1_1	-	+

# NOTE this script was writting based on the output of an older version of class1_class1_pairs.R 
# with the output file class1_minus_plus.txtand class1_plus_minus.txt
input=$1

# strand 1 
 awk '{OFS="\t";
  split($1,a,"_"); start1=a[2];end1=a[3];
  split($2,b,"_"); start2=b[2];end2=b[3];
  print $6, start1, end1, $1, 0, $14, $7, $8"\n"$3, start2, end2, $2, 0, $13, $4, $5}' $input > ucsc_class1TNEpairs.bed


#!/bin/bash

# script to determine the location of the peak transcription regions and distance
# in regards to transcription peaks on the opposite strand 

## running analysis 
#NOTE: I actually ran this on Terra because bigwig file is very large 
# see https://console.cloud.google.com/storage/browser/fc-secure-4300ea7f-8e77-4f20-a20b-38417286eaaf/338622bb-9e05-4d8b-9bb0-6f7faf91f6a7?authuser=rosan2027@gmail.com 
# see https://console.cloud.google.com/storage/browser/fc-secure-4300ea7f-8e77-4f20-a20b-38417286eaaf/57661aea-bb0e-417e-9e57-1ded396f4c25?authuser=rosan2027@gmail.com
#this script is here for documentation purposes
plusbed=""
minusbed=""

plusbw=""
minusbw=""

name=""

# run this TWICE (once for plus, once for minus)
# run strands seperately
# narrowPeak files are 0- coordinate based 
cat $plusbed | while read chr start end name rest
do
  l=$(expr $end - $start)
  bigWigSummary ~{input_bw} $chr $start $end $l | awk -vstrand=~{strand} -vchr=$chr -vstart=$start -vend=$end -vname=$name '{OFS="\t"; for(i=1;i<=NF;i++) if($i!="n/a" && $i>max) {imax=i;max=$i}}END{print chr, start+imax-1, start+imax, name, 0, strand, max, -1, -1, 0}'
done > $name.narrowPeak


# bedtools closest commmand 
# run after the narrowPeak file for each strand is generated 

# note: there closest TNE region on the opposite strand differs based on the order you compare your strands in 
# direction was always set to the plus strand (ex. if minus TNE is before plus TNE, distance is negative)
bedtools closest -a peaks.n200.minus.narrowPeak -b peaks.n200.plus.narrowPeak -t all -D b > minus.plus.bed # looking for neg 
bedtools closest -a peaks.n200.plus.narrowPeak -b peaks.n200.minus.narrowPeak -t all -D a > plus.minus.bed # looking for neg 

## checking outputs 
# tested the outputs of bedtools closest manually for the first 5 regions (not shown here)
# now testing based on a previous script I wrote (peak-characterization.sh)

#peak-characterization.sh 
# NOTE: narrowPeak input files for peak-characterization are in a different format 
# than the output of the narrowPeak format generated above 

# narrowPeak for the peak-characterization file gives the TNE coordinates as bed coordinates 

./peak-characterization.sh

#OUTPUTS: 

#input files were subset to the first five lines of each file 
# closest TNE analysis was done only with the first 5 lines 
# plus-minus.bed 
#chr1_10000_10259	87027 note: these genomic coordinates are only for the plus strand 
#chr1_96648_97228	375
#chr1_97665_97853	445
#chr1_98139_98291	29
#chr1_99751_100244	491

#minus-plus.bed 
#chr1_96712_97331	375 note: these genomic coordinates are only for the minus strand 
#chr1_98173_98337	29
#chr1_99993_100535	491
#chr1_105069_105226	5268
#chr1_107595_107738	7808

#first 5 lines of the full size closest TNE analysis
head -n5 plus.minus.bed
#chr1	10082	10083	chr1_10000_10259	0	+	6.3527	-1	-1	0	chr1	97109	97110	chr1_96712_97331	0	-	0.299824	-1	-1	0	87027
#chr1	96734	96735	chr1_96648_97228	0	+	0.443508	-1	-1	0	chr1	97109	97110	chr1_96712_97331	0	-	0.299824	-1	-1	0	375
#chr1	97807	97808	chr1_97665_97853	0	+	0.137655	-1	-1	0	chr1	98252	98253	chr1_98173_98337	0	-	0.196891	-1	-1	0	445
#chr1	98223	98224	chr1_98139_98291	0	+	0.204462	-1	-1	0	chr1	98252	98253	chr1_98173_98337	0	-	0.196891	-1	-1	0	29
#chr1	99838	99839	chr1_99751_100244	0	+	0.159993	-1	-1	0	chr1	100329	100330	chr1_99993_100535	0	-	0.118081	-1	-1	0	491

head -n5 minus.plus.bed 
#chr1	97109	97110	chr1_96712_97331	0	-	0.299824	-1	-1	0	chr1	96734	96735	chr1_96648_97228	0	+	0.443508	-1	-1	0	375
#chr1	98252	98253	chr1_98173_98337	0	-	0.196891	-1	-1	0	chr1	98223	98224	chr1_98139_98291	0	+	0.204462	-1	-1	0	29
#chr1	100329	100330	chr1_99993_100535	0	-	0.118081	-1	-1	0	chr1	99838	99839	chr1_99751_100244	0	+	0.159993	-1	-1	0	491

#NOTE these lines differ to the output because they paired to values outside of the first five lines in the plus strand file 
#chr1	105106	105107	chr1_105069_105226	0	-	0.130886	-1	-1	0	chr1	105096	105097	chr1_105014_105220	0	+	0.146664	-1	-1	0	10
#chr1	107646	107647	chr1_107595_107738	0	-	0.13321	-1	-1	0	chr1	107491	107492	chr1_107426_107639	0	+	0.177124	-1	-1	0	155


## how many equal distance peaks? 

#determine number of tnes based on narrowPeak file
wc -l peaks.n200.minus.narrowPeak
# 111390 peaks.n200.minus.narrowPeak

wc -l peaks.n200.plus.narrowPeak
# 110623 peaks.n200.plus.narrowPeak

# determine number of pairings 
wc -l minus.plus.bed
# 111440 minus.plus.bed

wc -l plus.minus.bed 
# 110671 plus.minus.bed

#111440 - 111390 = 50 
# 50 minus peak instances are multimapped 

#110671 - 110623 = 48 
# 48 plus peak instances are multimapped 






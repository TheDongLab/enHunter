#!/bin/bash

./scripts/peaks/bidirectional/bidirectional_pairs.sh /Users/rosanwang/Documents/donglab/projects/eRNA/enHunter/input_files/closest/plus.minus.bed /Users/rosanwang/Documents/donglab/projects/eRNA/enHunter/input_files/closest/minus.plus.bed 900

Rscript scripts/peaks/bidirectional/bidirectional_classes.R input_files/eRNA.characterize.feature.color.xls scripts/peaks/bidirectional/pairs/900/bidirectional_pairs.txt 



### making a list of ALL unique bidirectional pairs 
plusminus=/Users/rosanwang/Documents/donglab/projects/eRNA/enHunter/input_files/closest/plus.minus.bed
minusplus=/Users/rosanwang/Documents/donglab/projects/eRNA/enHunter/input_files/closest/minus.plus.bed

paste <(cut -f11-20 $minusplus) <(cut -f1-10 $minusplus) <(cut -f21 $minusplus) > tmp
cat tmp $plusminus | sort | uniq | cut -f4,6,14,16,21 > all.bidir.pairs.merged.bed

Rscript scripts/peaks/bidirectional/bidirectional_classes.R input_files/eRNA.characterize.feature.color.xls all.bidir.pairs.merged.bed
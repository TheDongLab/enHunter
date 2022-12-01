#!/bin/bash 
## create a list of TNEs for all classes (class I , class II, class III)
## note: plus and minus strand TNEs are not in seperate files 


awk '$24==1 {print $1}' ./input_files/characterization/eRNA.characterize.feature.color.xls > class1_TNEs.txt 
awk '$24==2 {print $1}' ./input_files/characterization/eRNA.characterize.feature.color.xls > class2_TNEs.txt
awk '$24==3 {print $1}' ./input_files/characterization/eRNA.characterize.feature.color.xls > class3_TNEs.txt
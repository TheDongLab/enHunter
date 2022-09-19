#!/bin/bash

#https://unix.stackexchange.com/questions/110645/select-lines-from-text-file-which-have-ids-listed-in-another-file 

#class1 minus pairs 

# number of minus strand TNEs in class1_minus_plus_test.txt
cut -f 2 class1_minus_plus_test.txt | uniq > minus_plus_minusTNEs.txt
# number of minus strand TNEs in class1_plus_minus_test.txt
cut -f 8 class1_plus_minus_test.txt | uniq > plus_minus_minusTNEs.txt

cat minus_plus_minusTNEs.txt plus_minus_minusTNEs.txt | sort | uniq > class_1_minus.txt

rm minus_plus_minusTNEs.txt
rm plus_minus_minusTNEs.txt

## DOES NOT WORK
# cut -f 8 class1_plus_minus_test.txt | uniq | cat stdin <(cut -f 2 class1_minus_plus_test.txt | uniq ) | sort | uniq > class_1_minus_test.txt

#class1 plus pairs 

# number of plus strand TNEs in class1_minus_plus_test.txt
cut -f 8 class1_minus_plus_test.txt | uniq > minus_plus_plusTNEs.txt
# number of plus strand TNEs in class1_plus_minus_test.txt
cut -f 2 class1_plus_minus_test.txt | uniq > plus_minus_plusTNEs.txt

cat minus_plus_plusTNEs.txt plus_minus_plusTNEs.txt | sort | uniq > class_1_plus.txt

rm minus_plus_plusTNEs.txt
rm plus_minus_plusTNEs.txt

#cut -f 1 class1_plus_minus.txt | uniq | cat <(cut -f 2 class1_minus_plus.txt | uniq ) | sort | uniq > class_1_plus.txt

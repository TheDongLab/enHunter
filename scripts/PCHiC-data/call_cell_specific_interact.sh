#!/bin/bash

#Mon     #F52D24
#Mac0    #BB312A
#Mac1    #F29809
#Mac2    #FFC761
#Neu     #EFCB09
#MK      #D0E009
#EP      #C5D84A
#Ery     #01B850
#FoeT    #01A689
#nCD4    #02A3A3
#tCD4    #038Ab0
#aCD4    #7AD2F7
#naCD4   #096EBE
#nCD8    #8D00DE
#tCD8    #CB96C4
#nB      #C40093
#tB      #E10064

for cell in Mon Mac0 Mac1 Mac2 Neu MK EP Ery FoeT nCD4 tCD4 aCD4 naCD4 nCD8 tCD8 nB tB 
do 
  ./cell_specific_interact.sh $cell
done 


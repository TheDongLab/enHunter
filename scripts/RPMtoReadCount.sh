# converting the eRNA.meanRPM.xls output from the enHunter pipeline to raw read counts 

awk 'OFS="\t" {if (NR==1) print $0; else {
split($1, a, "_"); 
len=a[3]-a[2]; 

for(i=2; i <=NF; i++) {
  x = $i * len
  ival = int($i * len)  # integer part, int() truncates

   # see if fractional part
   if (ival == x) {
     $i = ival
   } else if (x < 0) {
     $i = "NEGATIVE"
    }
   else {
   # assums that x is NEVER negative! 
       fraction = x - ival
    if (fraction >= .5) {
        $i= ival + 1
     }
     else {
      $i = ival
    }
   }
   }
  print $0
}
 }' eRNA.meanRPM.xls > round_test.xls

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

#### puts ALL cell type interactions into one file 
#### colored based on the cell type with the highest score 
awk 'OFS="\t" {
  if(NR==1){
    split($0, array); 
    print "track type=interact name=interact_PCHiC description=An_interact_file interactDirectional=true maxHeightPixels=200:100:50 visibility=full"
  } else {
    m=$12;
    col=12; 
    for(i=12;i<=28;i++)if($i>m){m=$i;col=i}; 
    
    if ($1 != $6) {
      start = $2; 
      end = $3;
    } else {
      if($2 < $7) {
        start=$2;
        end=$8;
      } else {
        start=$7; 
        end=$3;
      }
    }

    t=array[col];
    
    if (t == "Mon") color = "#F52D24";
    else if (t == "Mac0") color = "#BB312A";
    else if (t == "Mac1") color = "#F29809";
    else if (t == "Mac2") color = "#FFC761";
    else if (t == "Neu") color = "#EFCB09";
    else if (t == "MK") color = "#D0E009";
    else if (t == "EP") color = "#C5D84A"; 
    else if (t == "Ery") color = "#01B850"; 
    else if (t == "FoeT") color = "#01A689"; 
    else if (t == "nCD4") color = "#02A3A3"; 
    else if (t == "tCD4") color = "#038Ab0"; 
    else if (t == "aCD4") color = "#7AD2F7"; 
    else if (t == "naCD4") color = "#096EBE";
    else if (t == "nCD8") color = "#8D00DE";
    else if (t == "tCD8") color = "#CB96C4"; 
    else if (t == "nB") color = "#C40093"; 
    else if (t == "tB") color = "#E10064"; 
    else color = "ERROR"
    
    print "chr"$1, start, end, "chr"$1"/"start"/"end"/"t, 0, m, t, color, "chr"$1, $2, $3, $4, ".", "chr"$6, $7, $8, $9, "."
  }
}' PCHiC_peak_matrix_cutoff5.txt > test.interact

#### seperating PCHiC interactions based on cell type 





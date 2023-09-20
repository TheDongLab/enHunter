#!/usr/bin/awk -f
# usage:  q
# inputs: 
# headers = columns to select for 
# file = matrix file with columns for selection

# taken from https://unix.stackexchange.com/questions/529060/how-to-select-columns-in-a-file-by-the-column-name-in-a-txt-file-in-unix 

BEGIN {
  OFS="\t";
  ORS="\n";
}
NR==FNR {
    # loop through first file
    groups[++numGroups] = $1 
    next
}
FNR==1 {
    # second file
    for (i=7; i<=NF; i++) {
        # save each column into array
        f[$i] = i
    }
}
{
    printf "%s%s", $2, OFS
    for (groupNr=1; groupNr<=numGroups; groupNr++) {
        group = groups[groupNr]
        if (f[group] != "") printf "%s%s", $(f[group]), (groupNr<numGroups ? OFS : "")
        if (groupNr==numGroups) printf "%s", ORS
    }
}


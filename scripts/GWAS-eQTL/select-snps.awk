#!/usr/bin/awk -f

NR==FNR{a[$1]++; next} FNR==1&&NR!=FNR{print $0} $1 in a
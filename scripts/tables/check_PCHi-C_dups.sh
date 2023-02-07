# used to check number of duplicate rows for TNEs to answer the question: 
# are the duplicate TNE rows coming from multiple data entries or join errors?

module load bedtools/2.20.1 
cd /data/bioinformatics/projects/donglab/AMPPD_eRNA/output/

#intersect=$(cat <(groupBy -i plus/eRNA.PCHiC/eRNA.plus.f22.PCHiCPromoters.bait.score.tmp -g 1 -c 2 -o count)  <(groupBy -i minus/eRNA.PCHiC/eRNA.minus.f22.PCHiCPromoters.bait.score.tmp -g 1 -c 2 -o count))
#target_table=$(groupBy -i target_table/TNE.PCHiC.xls -g 1 -c 2 -o count)

cat <(cat <(sort -k 1,1 plus/eRNA.PCHiC/eRNA.plus.f22.PCHiCPromoters.bait.score.tmp | groupBy -g 1 -c 2 -o count)  <(sort -k 1,1 minus/eRNA.PCHiC/eRNA.minus.f22.PCHiCPromoters.bait.score.tmp | groupBy -g 1 -c 2 -o count)) <(sort -k1,1 target_table/TNE.PCHiC.xls | groupBy -g 1 -c 2 -o count) | sort | uniq -u | wc -l 
# this should onyl return 1 due to the header 

### these two TNEs are repeated because they have the same coordinates for the plus and minus strand! 
# chr3_112462880_112463020        6
# chr6_42262460_42262610  2

# so grouped TNE for plus and TNE for minus produce the same row, which is then filtered out with uniq -u 
# but TNE target table row is kept 


### checking that there are no duplicate baitID and oeID pairs
cut -f 2,3 eRNA.minus.f22.PCHiCPromoters.bait.score.tmp | sort -k 1,1 | uniq | wc -l 
wc -l  eRNA.minus.f22.PCHiCPromoters.bait.score.tmp

## hahahaha this doesn't work because one oe Fragment can over lap multiple TNEs, causing baitID and oeID to be repeated 
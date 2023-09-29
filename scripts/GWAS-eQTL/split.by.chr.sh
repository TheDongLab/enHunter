#!/bin/bash
# used to split the inputs of matrixEQTL by chromosome for faster processing 

for i in {1..22} "X" "Y"
do 
  echo "chr$i"
  [ -e "chr$i".gene.eRNA.coordinates ] || cat corheader <(awk -v v=chr$i '$2 == v' gene.eRNA.coordinates) > "chr$i".gene.eRNA.coordinates
  # [ -e "chr$i".eRNA.gene.Matrix.reordered.txt ] || cat header <(grep -f <(cut -f 1 "chr$i".gene.eRNA.coordinates) -F eRNA.gene.Matrix.reordered.txt) > "chr$i".eRNA.gene.Matrix.reordered.txt
done 

## did not include sex chromosomes for eQTL 
for i in {2..22}
do 
  echo "chr$i"
  #[ -e chr$i.snps.coordinates ] || cat snps.header <(awk -v v=chr$i '$2 == v' snps.coordinates) > chr$i.snps.coordinates
  #[ -e chr$i.snps.tmp ] || cut -f 1 chr$i.snps.coordinates | sed '1d' > chr$i.snps.tmp
  [ -e chr$i.eRNA.SNP_Matrix.revised.txt ] || bsub -q normal -n 1 -M 1000 "awk -f select-snps.awk chr$i.snps.tmp eRNA.SNP_Matrix.revised.txt > chr$i.eRNA.SNP_Matrix.revised.txt"
done 

rm *.tmp

# running matrixEQTL 
for i in {2..22}
do 
  echo "running matrix eQTL for chr$i"
  [ -e chr$i.all.cis.eQTL.txt ] || bsub -q bigmem -n 1 -M 35000 "Rscript matrixeQTL.inputs.R snps/chr$i.eRNA.SNP_Matrix.revised.txt snps/chr$i.snps.coordinates genes/chr$i.eRNA.gene.Matrix.reordered.txt genes/chr$i.gene.eRNA.coordinates covariate.reordered.txt chr$i.cis.eQTL"
done 

bsub -q bigmem -n 1 -M 35000 "Rscript matrixeQTL.inputs.R snps/chr1.eRNA.SNP_Matrix.revised.txt snps/chr1.snps.coordinates genes/chr1.eRNA.gene.Matrix.reordered.txt genes/chr1.gene.eRNA.coordinates covariate.reordered.txt chr1.cis.eQTL"


# splitting eRNA and genes into seperate files 
for i in {2..22}
do 
  echo "seperating chr$i into eRNA and genes "
  bsub -q normal -n 1 -M 2000 awk 'NR==1{print $0} $2 ~ /^chr/' chr$i.cis.eQTL.txt | sort -k 1 > eRNA.chr$i.cis.eQTL.sorted.txt 
  bsub -q normal -n 1 -M 2000 awk 'NR==1{print $0} $2 ~ /^ENSG/' chr$i.cis.eQTL.txt | sort -k 1 > gene.chr$i.cis.eQTL.sorted.txt 
done 

bsub -q normal -n 1 -M 2000 "sort -k 1 gene.chr1.cis.eQTL.txt > gene.chr1.cis.eQTL.sorted.txt"
bsub -q normal -n 1 -M 2000 "sort -k 1 eRNA.chr1.cis.eQTL.txt > eRNA.chr1.cis.eQTL.sorted.txt"

# join each cis eQTL file with the SNP location 
# test batch with chr1
join gene.chr1.cis.eQTL.sorted.txt snps/chr1.snps.coordinates

#bsub -q normal -n 1 -M 1000 "awk -f select-snps.awk chr1.snps.tmp eRNA.SNP_Matrix.revised.txt > chr1.eRNA.SNP_Matrix.revised.tmp"
cat header chr$i.eRNA.SNP_Matrix.revised.tmp > chr$i.eRNA.SNP_Matrix.revised.txt
### GREP IS NOT WORKING@@@!!!!!!!
awk 'NR==FNR{a[$1]++; next} $1 in a' chr1.snps.tmp eRNA.SNP_Matrix.revised.txt > chr1.eRNA.SNP_Matrix.revised.tmp





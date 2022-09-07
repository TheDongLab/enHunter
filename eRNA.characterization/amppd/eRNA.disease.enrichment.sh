TNE_file=
total_file=/data/bioinformatics/external_data/externalData/GWAS_catalog_20220810/GWAS_20220810.v1.02.bed

## groupby doesn't work on newer versions of bedtools 

module load bedtools/2.23.0
sort -k 5 snp_gtexDapg_sorted.bed | bedtools groupby -g 5 -c 5 -o count > snp_gtexDapg.counts

sort -k 5 snp_gtexDapg_sorted.bed | bedtools groupby -g 5 -c 5 -o count > snp_gtexDapg.counts

sort -k 5 all_signif_varient_gene_pairs_sorted.bed | bedtools groupby -g 5 -c 5 -o count > all_signif_varient_gene_pairs.counts

sort -k 6 all_enhancer_tissues_hg38.bed | bedtools groupby -g 6 -c 6 -o count > all_enhancer_tissues_hg38.counts

sort -k 5 GWAS_20220810.v1.02.bed | bedtools groupby -g 5 -c 5 -o count > GWAS_20220810.v1.02.counts
sort -k 1 GWAS_20220810.v1.02.counts | bedtools groupby -g 1 -c 2 -o sum > GWAS_20220810.v1.02.counts.v2
sort -k 1 GWAS_20220810.v1.02.counts.v2 | bedtools groupby -g 1 -c 2 -o sum > GWAS_20220810.v1.02.counts.v3
# to check the groupby 
cut -f 2 all_signif_varient_gene_pairs.counts | paste -sd+ | bc

cut -f 2 GWAS_20220810.v1.02.counts | paste -sd+ | bc


######## minus strand 
sort -k 9 eRNA.minus.f16.GWASDisease.txt | bedtools groupby -g 9 -c 9 -o count > eRNA.minus.f16.GWASDisease.counts
cut -f 2 eRNA.minus.f16.GWASDisease.counts | paste -sd+ | bc

sort -k 1 eRNA.minus.f16.GWASDisease.counts | bedtools groupby -g 1 -c 2 -o sum > eRNA.minus.f16.GWASDisease.counts.v2

sort -k 9 eRNA.minus.f18.eSNP.gtexCaviarDisease.txt | bedtools groupby -g 9 -c 9 -o count > eRNA.minus.f18.eSNP.gtexCaviarDisease.counts
cut -f 2 eRNA.minus.f18.eSNP.gtexCaviarDisease.counts | paste -sd+ | bc

sort -k 9 eRNA.minus.f18.eSNP.gtexDapgDisease.txt | bedtools groupby -g 9 -c 9 -o count > eRNA.minus.f18.eSNP.gtexDapgDisease.counts
cut -f 2 eRNA.minus.f18.eSNP.gtexDapgDisease.counts | paste -sd+ | bc

sort -k 9 eRNA.minus.f18.eSNP.pvalDisease.txt | bedtools groupby -g 9 -c 9 -o count > eRNA.minus.f18.eSNP.pvalDisease.counts
cut -f 2 eRNA.minus.f18.eSNP.pvalDisease.counts | paste -sd+ | bc

######### plus strand 
module load bedtools/2.23.0
sort -k 9 ../eRNA.plus.f16.GWASDisease.txt | bedtools groupby -g 9 -c 9 -o count > eRNA.plus.f16.GWASDisease.counts

sort -k 9 ../eRNA.plus.f18.eSNP.gtexCaviarDisease.txt | bedtools groupby -g 9 -c 9 -o count > eRNA.plus.f18.eSNP.gtexCaviarDisease.counts

sort -k 9 ../eRNA.plus.f18.eSNP.gtexDapgDisease.txt | bedtools groupby -g 9 -c 9 -o count > eRNA.plus.f18.eSNP.gtexDapgDisease.counts

sort -k 9 ../eRNA.plus.f18.eSNP.pvalDisease.txt | bedtools groupby -g 9 -c 9 -o count > eRNA.plus.f18.eSNP.pvalDisease.counts

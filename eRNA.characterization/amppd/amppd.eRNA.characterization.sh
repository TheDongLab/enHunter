## script to characterize eRNAs with the following features
# usage: $pipeline_path/src/eRNA.characterize.sh minus
#
#Features	            Data type	    Description
#========             =========     ==================================
#dis2TSS	            integer	      distance between the middle points of HiTNE and the nearest TSS, in bp. If HiTNE is located intronic, it's plus; otherwise it's minus. 
#RPKM	                float	        normalized expression level, calculated in the same way as RPKM.
#RPM	                float	        reads density at the summit position, normalized to total mapped reads in million. 
#readsCount	          integer	      raw reads count mapped to the HiTNE
#normCpG	            float	        normalized CpG score 
#nTFBS	              integer	      number of distinct TFs bound to the HiTNE, based on ENCODE ChIPseq cluster (wgEncodeRegTfbsClusteredV3)
#P300	                boolean	      if the P300 binding site found in the HiTNE
#enhancer_CAGE	      boolean	      if overlap with any CAGE-defined permissive enhancers
#enhancer_histone	    boolean	      if overlap with any histone marks-defined enhancers (chromHMM states of E6|E7|E12 from substantial nigro)
#enhancer_VISTA	      boolean	      if overlap with any tested enhancers (positive enhancers from VISTA enhancer database)
#bDNase	              boolean	      if overlap with DNase cluster from ENCODE (wgEncodeRegDnaseClustered V2)
#DNaseROADMAP         boolean       if overlap with DNase narrow peak from Roadmap brain samples 
#conservation	        float	        mean phastCons score for the HiTNE region
#bHCNE	              boolean	      if overlapping with any HCNEs (HCNE_hg38_danRer7_70pc_50col from Ancora)

#bidirectional_trans  str,int,int   closest peak on opposite strand, distance away, and directionality score  

#GWAS	                integer	      number of GWAS SNPs in HiTNE
#bGWAS	              boolean	      if any GWAS SNPs in HiTNE
#eSNP	                integer	      number of eQTL SNPs in HiTNE
#beSNP	              boolean	      if any eQTL SNPs in HiTNE

#nHostgene	          integer	      number of HiTNEs in the host gene. 0 for intergenic HiTNE
#lenHostgene	        integer	      length of host gene. 0 for intergenic HiTNE

#HiC                  integer       number of PIR (promoter interacting regions) in TNE 
#bHiC                 boolean       if any PIR 

pipeline_path=/data/bioinformatics/projects/donglab/AMPPD_eRNA
source $pipeline_path/config.txt

STRAND=$1
echo "$STRAND"
inputBG=$pipeline_path/inputs/$STRAND/combined.mean.normalized.random.samplesN200.$STRAND.bigwig

EXTERNAL_FEATURE=/data/bioinformatics/external_data/externalData
inputbed=$pipeline_path/inputs/$STRAND/eRNA.bed

[ -e eRNA.$STRAND.random.bed ] || bedtools random -n 100000 -l 400 -seed 1234 -g $pipeline_path/inputs/ChromInfo.txt | intersectBed -a stdin -b <(cat $pipeline_path/inputs/toExclude.bed $inputbed | cut -f1-3 | sortBed | mergeBed -i -) -v > eRNA.$STRAND.random.bed

# ====================================
## dis2TSS (distance btw middle of HiTNE and the nearest TSS)
# ====================================
echo "RUNNING ---- dis2TSS"
## NOTE gencode.v37.annotation.genes.bed is the same as gencode.v37.annotation.genes.gtf.bed
# fgrep -w gene $GENOME/Annotation/Genes/gencode.v37.annotation.gtf | sed 's/[;"]//g'  | awk '{OFS="\t"; print $1, $4-1, $5, $18"___"$10"___"$14, 0, $7}' > $GENOME/Annotation/Genes/gencode.v37.annotation.genes.bed
# have to deal with intronic and intergenic separately
# intronic ones (if located in two genes' intron, just randomly pick the first hit in the file.)
TMPFILE=`mktemp /tmp/example.XXXXXXXXXX`
# $11 == 0 means TNE is located within a gene -> + value 
awk '{OFS="\t"; mid=int(($3+$2)/2); print $1, mid, mid+1,$4}' $inputbed | closestBed -a - -b <(sortBed -i $GENOME/Annotation/Genes/gencode.v37.annotation.genes.bed) -d -t first | awk '$11==0' | awk '{OFS="\t"; tss=($10=="+")?$6:$7; d=tss-$2; if(d<0) d=-d; print $4, d;}' > $TMPFILE
# intergenic ones 
# $11 != 0 means TNE is located not within a gene -> code converts to - value 
awk '{OFS="\t"; mid=int(($3+$2)/2); print $1, mid, mid+1,$4}' $inputbed | closestBed -a - -b <(sortBed -i $GENOME/Annotation/Genes/gencode.v37.annotation.genes.bed) -d -t first | awk '$11!=0' | cut -f1-4 | sort -k1,1 -k2,2n -u | closestBed -a - -b <(awk '{OFS="\t"; tss=($6=="+")?$2:($3-1);  print $1, tss, tss+1, $4, $3-$2, $6}' $GENOME/Annotation/Genes/gencode.v37.annotation.genes.bed | sortBed) -D b -t first | awk '{OFS="\t"; print $4,($11<0)?$11:-$11;}' >> $TMPFILE

sort -k1,1 $TMPFILE > eRNA.$STRAND.f01.dis2TSS.txt

# ====================================
## RPKM (Note: this RPKM is different from the normal RPKM, they might have a factor of read length difference)
# ====================================
# mean0: average over bases with non-covered bases counting as zeroes
echo "RUNNING ---- RPKM"
bigWigAverageOverBed $inputBG $inputbed stdout | cut -f1,5 | sort -k1,1 | awk '{OFS="\t"; print $1, $2*1000+0}' > eRNA.$STRAND.f02.RPKM.txt

# ====================================
### RPM 
# ====================================
echo "RUNNING ---- RPM"
bigWigAverageOverBed $inputBG $inputbed stdout -minMax | cut -f1,8 | sort -k1,1  > eRNA.$STRAND.f03.RPM.txt

# ====================================
# normalized CpG score
# ====================================
echo "RUNNING ---- CpG score"
bedtools getfasta -name -tab -fi $GENOME/Sequence/WholeGenomeFasta/genome.fa -bed $inputbed -fo eRNA.$STRAND.seq.tab
$pipeline_path/bin/getNormalizedCpGscore.awk eRNA.$STRAND.seq.tab | sort -k1,1  > eRNA.$STRAND.f05.CpG.txt
#textHistogram -col=2 -real -binSize=0.02 -maxBinCount=50 -minVal=0 eRNA.f05.CpG.tab

bedtools getfasta -name -tab -fi $GENOME/Sequence/WholeGenomeFasta/genome.fa -bed eRNA.$STRAND.random.bed -fo random.$STRAND.seq.tab
$pipeline_path/bin/getNormalizedCpGscore.awk random.$STRAND.seq.tab | sort -k1,1 > random.$STRAND.f05.Cp.txt
#textHistogram -col=2 -real -binSize=0.02 -maxBinCount=50 -minVal=0 random.f05.CpG.tab

# ERROR: promoters.$STRAND.seq.txt is empty  
grep protein_coding.protein_coding $GENOME/Annotation/Genes/gencode.v37.annotation.bed12 | awk '{OFS="\t"; s=($6=="+")?($2-200):($3-200); if(s<0) s=0; print $1,s,s+400,$4}' | bedtools getfasta -name -tab -fi $GENOME/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo promoters.$STRAND.seq.txt

$pipeline_path/bin/getNormalizedCpGscore.awk promoters.$STRAND.seq.txt | sort -k1,1 > promoters.$STRAND.f05.CpG.txt
#textHistogram -col=2 -real -binSize=0.02 -maxBinCount=50 -minVal=0 promoters.f05.CpG.tab

# ====================================
# TFBS
# ====================================
## TFBS data for 161 transcription factors in 91 cell types was downloaded from ENCODE
# curl -s http://hgdownload.soe.ucsc.edu/goldenPath/hg38/encRegTfbsClustered/encRegTfbsClusteredWithCells.hg38.bed.gz | gunzip | awk '{OFS="\t"; $6=1+gsub(",",",",$6); print}' > count.encRegTfbsClusteredWithCells.hg38.bed
# number of TF ChIP-seq peaks in the region (only if it occurs in at least one cell lines ## only TF peaks supported by >=5 cell lines are counted)
echo "RUNNING ---- TFBS"
awk '$6>=1' $EXTERNAL_FEATURE/TFBS/count.encRegTfbsClusteredWithCells.hg38.bed | intersectBed -a $inputbed -b stdin -c | sort -k4,4 | cut -f4,5 > eRNA.$STRAND.f06.TFBS.txt

# Define TF hotspot (regions with >=5 different TFs bound)
# might be used somewhere else -> dont need 
#sort -k1,1 -k2,2n $EXTERNAL_FEATURE/TFBS/count.encRegTfbsClusteredWithCells.hg38.bed | mergeBed -i stdin -c 4 -o count_distinct | awk '$4>=5' > $EXTERNAL_FEATURE/TFBS/wgEncodeRegTfbsClusteredWithCellsV3.hotspot.bed 

# ====================================
# P300 binding
# ====================================
# if any P300 biding site found in the region
awk '$4=="EP300"' $EXTERNAL_FEATURE/TFBS/count.encRegTfbsClusteredWithCells.hg38.bed | intersectBed -a $inputbed -b stdin -c | sort -k4,4 | cut -f4,5 > eRNA.$STRAND.f07.P300.txt

# ====================================
# CAGE-defined enhancers
# ====================================
echo "RUNNING ---- CAGE"
# annotation with blood expressed enhancers 
intersectBed -a $inputbed -b $EXTERNAL_FEATURE/CAGE/blood_expressed_enhancers_hg38.bed -c | sort -k4,4 | cut -f4,5 > eRNA.$STRAND.f08.CAGEbloodenhancer.txt

# running Fisher test with differentially expressed enhancer tissue facets 
intersectBed -a $inputbed -b $EXTERNAL_FEATURE/CAGE/diff_expressed_enh/all_diff_exp_enhancer_tissues_hg38_sorted.bed -sorted -wb | sort -k4,4  > eRNA.$STRAND.f08.CAGEenhtissue.txt
sort -k 10 eRNA.$STRAND.f08.CAGEenhtissue.txt | bedtools groupby -g 10 -c 10 -o count > eRNA.$STRAND.f08.CAGEenhtissue.counts

# all premissive enhancers -> liftover to hg38 
# curl -s https://fantom.gsc.riken.jp/5/datafiles/latest/extra/Enhancers/human_permissive_enhancers_phase_1_and_2.bed.gz | gzip -d | liftOver stdin hg19ToHg38.over.chain.gz human_permissive_enhancers_phase_1_and_2_hg38.bed unMappedHumanPermissiveEnh # N=65386
intersectBed -a $inputbed -b $EXTERNAL_FEATURE/CAGE/human_permissive_enhancers_phase_1_and_2_hg38.bed -c | sort -k4,4 | cut -f4,5 > eRNA.$STRAND.f08.CAGEenhancer.1n2.txt

# ====================================
# overlap with histone-defined enhancers
# ====================================

# Roadmap Epigenomics enhancers (http://egg2.wustl.edu/roadmap/web_portal/meta.html)
#E062 - Primary mononuclear cells from peripheral blood  
#E034 - Primary T cells from peripheral blood
#E045 - Primary T cells effector/memory enriched from peripheral blood
#E033 - Primary T cells from cord blood
#E044 - Primary T regulatory cells from peripheral blood
#E043 - Primary T helper cells from peripheral blood
#E039 - Primary T helper naive cells from peripheral blood
#E041 - Primary T helper cells PMA-I stimulated
#E042 - Primary T helper 17 cells PMA-I stimulated
#E040 - Primary T helper memory cells from peripheral blood 1
#E037 - Primary T helper memory cells from peripheral blood 2
#E048 - Primary T CD8+ memory cells from peripheral blood
#E038 - Primary T helper naive cells from peripheral blood
#E047 - Primary T CD8+ naive cells from peripheral blood
#E029 - Primary monocytes from peripheral blood
#E031 - Primary B cells from cord blood
#E035 - Primary hematopoietic stem cells
#E051 - Primary hematopoietic stem cells G-CSF-mobilized Male
#E050 - Primary hematopoietic stem cells G-CSF-mobilized Female
#E036 - Primary hematopoietic stem cells short term culture
#E032 - Primary B cells from peripheral blood
#E046 - Primary Natural Killer cells from peripheral blood
#for i in E062 E034 E045 E033 E044 E043 E039 E041 E042 E040 E037 E048 E038 E047 E029 E031 E035 E051 E050 E036 E032 E046; do 
#        echo $i
#        curl -o Segment/${i}_15_coreMarks_hg38lift_segments.bed.gz https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/${i}_15_coreMarks_hg38lift_segments.bed.gz 
#        gzip -d Segment/${i}_15_coreMarks_hg38lift_segments.bed.gz 
#done
#cat /data/bioinformatics/external_data/externalData/Histone/Segment/E*_15_coreMarks_hg38lift_segments.bed | awk '$4~/E6|E7|E12/' > 15_coreMarks_hg38lift_segments.E6E7E12.bed

# E6, E7, and E12 for enhancer (http://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html#core_15state)
echo "RUNNING ---- Histone"
intersectBed -a $inputbed -b $EXTERNAL_FEATURE/Histone/15_coreMarks_hg38lift_segments.E6E7E12.bed -c | sort -k4,4 | cut -f4,5 > eRNA.$STRAND.f09.chromHMM_blood.txt

# ====================================
# overlap with VISTA enhancers
# ====================================
echo "RUNNING ---- VISTA"
intersectBed -a $inputbed -b $EXTERNAL_FEATURE/VISTA/hg38.pos_regions.2022-10-27.bed -wao  | sort -k4,4 | awk '{OFS="\t"; print $4, ($6==-1)?"NA":$8"|"$9}' | groupBy -g 1 -c 2 -o concat > eRNA.$STRAND.f10.VISTA.txt
# only group subset which does not have NA 

# ====================================
# DNase cluster 
# ====================================
# download from ENCODE DNase cluster
# ( DNase cluster in V2: minimal score 500, at least 5 cell lines)
echo "RUNNING ---- DNase"
intersectBed -a $inputbed -b $EXTERNAL_FEATURE/DNase/ENCODE/hg38wgEncodeRegDnaseClustered.filtered.bed -c | sort -k4,4 | cut -f4,5 > eRNA.$STRAND.f12.DNaseENCODE.txt

# Roadmap DNase (http://egg2.wustl.edu/roadmap/web_portal/DNase_reg.html#delieation)
#for i in E062 E034 E045 E033 E044 E043 E039 E041 E042 E040 E037 E048 E038 E047 E029 E031 E035 E051 E050 E036 E032 E046; do echo $i; curl -o regions_enh/regions_enh_${i}.bed.gz https://egg2.wustl.edu/roadmap/data/byDataType/dnase/BED_files_enh/regions_enh_${i}.bed.gz;gzip -d regions_enh/regions_enh_${i}.bed.gz; liftOver regions_enh/regions_enh_${i}.bed hg19ToHg38.over.chain.gz regions_enh/regions_enh_${i}_hg38.bed regions_enh/unMapped_${i}; done
#cat regions_enh/regions_enh_*_hg38.bed | sortBed | mergeBed -i - > regions_enh_merged.hg38.blood.bed 
# same markers -> merge with bedtools

## TODO get the narrowPeak file -> peak information needed for strand peak (heatmap) 
# actually I think I can just run this with a bed file 
intersectBed -a $inputbed -b $EXTERNAL_FEATURE/DNase/ROADMAP/regions_enh_merged.hg38.blood.bed -c | sort -k4,4 | cut -f4,5 > eRNA.$STRAND.f12.DNaseROADMAP.txt

# ====================================
# Conservation - mean phyloP score
# ====================================
#rsync -avz --progress rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw ./
echo "RUNNING ---- Conservation phyloP"
bigWigAverageOverBed $EXTERNAL_FEATURE/Conservation/hg38.phyloP100way.bw $inputbed stdout | sort -k1,1 | cut -f1,6 > eRNA.$STRAND.f13.phyloP.txt

# ====================================
# Conservation - overlap with HCNE or not
# ====================================
# overlap with HCNE
# curl -s http://ancora.genereg.net/downloads/hg38/vs_zebrafish/HCNE_hg38_danRer7_70pc_50col.bed.gz | gunzip > HCNE_hg38_danRer7_70pc_50col.bed 
echo "RUNNING ---- HCNE"
intersectBed -a $inputbed -b $EXTERNAL_FEATURE/Conservation/Ancora/HCNE_hg38_danRer7_70pc_50col.bed -c | sort -k4,4 | cut -f4,5 > eRNA.$STRAND.f15.HCNE.txt

# ====================================
# GWAS - overlap with GWAS snps 
# ====================================
# overlap with GWAS snps 
# NOTE: when converting SNPs to bed file, keep in mind that bed files are on a 0 based coordinate system 

echo "RUNNING ---- GWAS"
intersectBed -a $inputbed -b $EXTERNAL_FEATURE/GWAS_catalog_20220810/GWAS_20220810.v1.02.bed -c -sorted | sort -k4,4 | cut -f4,5 > eRNA.$STRAND.f16.GWAS.txt
intersectBed -a $inputbed -b $EXTERNAL_FEATURE/GWAS_catalog_20220810/GWAS_20220810.v1.02.bed -sorted -wb | sort -k4,4  > eRNA.$STRAND.f16.GWASDisease.txt

# ====================================
# eQTL - overlap with eQTL snps 
# ====================================
# overlap with eQTL snps 
# converting the Caviar file 
#bigBedToBed $EXTERNAL_FEATURE/GTEx_hg38/gtexCaviar.bb stdout | awk '{OFS="\t"; print $1, $7, $8, $14, $19}' | sort -k1,1 -k2,2n | uniq > $EXTERNAL_FEATURE/GTEx_hg38/snp_gtexCaviar_sorted.bed

# converting the Dapg file 
#bigBedToBed $EXTERNAL_FEATURE/GTEx_hg38/gtexDapg.bb stdout | awk '{OFS="\t"; print $1, $7, $8, $14, $19 }' | sort -k1,1 -k2,2n | uniq > $EXTERNAL_FEATURE/GTEx_hg38/snp_gtexDapg_sorted.bed


# signif p-value eQTL file 
#mkdir all_tissues 
#touch ./all_tissues/all_signif_varient_gene_pairs.bed 

#for file in *signif_variant_gene_pairs.txt.gz
#do 
#    tissue="$(echo $file | cut -d"." -f1)"
#    echo $tissue 
#    zcat $file | awk -v tissue=$tissue 'NR>1 {OFS="\t"; split($1, a, "_"); print a[1], a[2] - 1, a[2], $1, tissue }' >> ./all_tissues/all_signif_varient_gene_pairs.bed 
#done 

sort -k1,1 -k2,2n all_signif_varient_gene_pairs.bed > all_signif_varient_gene_pairs_sorted.bed

echo "RUNNING ---- eQTL"
intersectBed -a $inputbed -b $EXTERNAL_FEATURE/GTEx_p_value/GTEx_hg38_UCSC_track/snp_gtexCaviar_sorted.bed -c | sort -k4,4 | cut -f4,5 > eRNA.$STRAND.f18.eSNP.gtexCaviar.txt
intersectBed -a $inputbed -b $EXTERNAL_FEATURE/GTEx_p_value/GTEx_hg38_UCSC_track/snp_gtexDapg_sorted.bed -c | sort -k4,4 | cut -f4,5 > eRNA.$STRAND.f18.eSNP.gtexDapg.txt

intersectBed -a $inputbed -b $EXTERNAL_FEATURE/GTEx_p_value/GTEx_hg38_UCSC_track/snp_gtexCaviar_sorted.bed -sorted -wb | sort -k4,4 > eRNA.$STRAND.f18.eSNP.gtexCaviarDisease.txt
intersectBed -a $inputbed -b $EXTERNAL_FEATURE/GTEx_p_value/GTEx_hg38_UCSC_track/snp_gtexDapg_sorted.bed -sorted -wb | sort -k4,4 > eRNA.$STRAND.f18.eSNP.gtexDapgDisease.txt

intersectBed -a $inputbed -b $EXTERNAL_FEATURE/GTEx_p_value/GTEx_Analysis_v8_eQTL/all_tissues/all_signif_varient_gene_pairs_sorted.bed -sorted -c | sort -k4,4 | cut -f4,5 > eRNA.$STRAND.f18.eSNP.pval.txt
intersectBed -a $inputbed -b $EXTERNAL_FEATURE/GTEx_p_value/GTEx_Analysis_v8_eQTL/all_tissues/all_signif_varient_gene_pairs_sorted.bed -sorted -wb | sort -k4,4 > eRNA.$STRAND.f18.eSNP.pvalDisease.txt

# ====================================
# number of HiTNEs in host genes.  # based on strand information 
# ====================================
# need a new bedtools version to run 
module load bedtools/2.26.0
# enforce strandedness 
inputbedstranded=$pipeline_path/inputs/$STRAND/eRNA_stranded.bed

# if overlap with multiple genes, take the longest one
#intersectBed -a $inputbed -b $GENOME/Annotation/Genes/gencode.v37.annotation.gtf.genes.bed -wao | awk '{OFS="\t"; print $0,$7-$6;}' | sort -k4,4 -k12,12nr | awk '{OFS="\t"; if($4!=id) {print; id=$4;}}' | cut -f1-4,8,12 | sed 's/\t\./\tNA/g' > $TMPFILE
intersectBed -a $inputbedstranded -b $GENOME/Annotation/Genes/gencode.v37.annotation.genes.bed -wao -s | awk '{OFS="\t"; print $0,$9-$8;}' | sort -k4,4 -k14,14nr | awk '{OFS="\t"; if($4!=id) {print; id=$4;}}' | cut -f1-4,10,14 | sed 's/\t\./\tNA/g' > $TMPFILE
# host gene (if any); otherwise NA
cut -f4-5 $TMPFILE > eRNA.$STRAND.f19.Hostgene.txt

grep -vw NA $TMPFILE | cut -f5 | sort | uniq -c | join -1 5 -2 2 -a 1 -e '0' -o '1.1,1.2,1.3,1.4,1.5,1.6,2.1' <(sort -k5,5 $TMPFILE) - | sed 's/ /\t/g' > $TMPFILE.2

cut -f4,7 $TMPFILE.2 | sort -k1,1 > eRNA.$STRAND.f20.nHostgene.txt

# ====================================
# length of host genes
# ====================================
cut -f4,6 $TMPFILE.2 | sort -k1,1 > eRNA.$STRAND.f21.lenHostgene.txt


## move files generated to output directory 
mv *.$STRAND.* $pipeline_path/output/$STRAND/

## merge into a big file: eRNA.characterize.xls
Rscript $pipeline_path/src/eRNA.characterize.merge.R `ls $pipeline_path/output/$STRAND/eRNA.$STRAND.f*.txt`

# ====================================
# Hi-C - overlap with PIR 
# ====================================
echo "RUNNING ---- Hi-C"
# TODO test this because it hasn't been run before 
intersectBed -a $inputbed -b $EXTERNAL_FEATURE/Hi-C/PCHiC_peak_matrix_cutoff5_hg38.bed -c | sort -k4,4 | cut -f4,5 > eRNA.$STRAND.f22.HiC.txt

#reports the intersected region of A and B with A's name 
intersectBed -a $inputbed -b $EXTERNAL_FEATURE/Hi-C/PCHiC_peak_matrix_cutoff5_hg38.bed -wb | sort -k4,4 > eRNA.$STRAND.f22.HiCPromoters.txt

exit;

# ============================
# merging data files for both strands 
# ============================

########## GWAS
cat $pipeline_path/output/minus/eRNA.enrichment/eRNA.minus.f16.GWASDisease.txt $pipeline_path/output/plus/eRNA.enrichment/eRNA.plus.f16.GWASDisease.txt | sort -k 9 | cut -f 5-9 | uniq | bedtools groupby -g 5 -c 5 -o count | sort -k 1 | bedtools groupby -g 1 -c 2 -o sum > $pipeline_path/output/merged/eRNA.f16.GWASDisease.counts
# just for testing 
#cat $pipeline_path/output/minus/eRNA.enrichment/eRNA.minus.f16.GWASDisease.txt $pipeline_path/output/plus/eRNA.enrichment/eRNA.plus.f16.GWASDisease.txt | sort -k 9 | cut -f 5-9 | bedtools groupby -g 5 -c 4 -o concat |less 

########## eQTL 
cat $pipeline_path/output/minus/eRNA.enrichment/eRNA.minus.f18.eSNP.gtexCaviarDisease.txt $pipeline_path/output/plus/eRNA.enrichment/eRNA.plus.f18.eSNP.gtexCaviarDisease.txt | sort -k 9 | cut -f 5-9 | uniq | bedtools groupby -g 5 -c 5 -o count | sort -k 1 | bedtools groupby -g 1 -c 2 -o sum > $pipeline_path/output/merged/eRNA.f18.eSNP.gtexCaviarDisease.txt
cat $pipeline_path/output/minus/eRNA.enrichment/eRNA.minus.f18.eSNP.gtexDapgDisease.txt $pipeline_path/output/plus/eRNA.enrichment/eRNA.plus.f18.eSNP.gtexDapgDisease.txt | sort -k 9 | cut -f 5-9 | uniq | bedtools groupby -g 5 -c 5 -o count | sort -k 1 | bedtools groupby -g 1 -c 2 -o sum > $pipeline_path/output/merged/eRNA.f18.eSNP.gtexDapgDisease.txt
cat $pipeline_path/output/minus/eRNA.enrichment/eRNA.minus.f18.eSNP.pvalDisease.txt $pipeline_path/output/plus/eRNA.enrichment/eRNA.plus.f18.eSNP.pvalDisease.txt | sort -k 9 | cut -f 5-9 | uniq | bedtools groupby -g 5 -c 5 -o count | sort -k 1 | bedtools groupby -g 1 -c 2 -o sum > $pipeline_path/output/merged/eRNA.f18.eSNP.pvalDisease.txt

########## CAGE
cat $pipeline_path/output/minus/eRNA.minus.f08.CAGEenhtissue.txt $pipeline_path/output/plus/eRNA.plus.f08.CAGEenhtissue.txt | sort -k 10 | cut -f 5-9 | uniq | bedtools groupby -g 5 -c 5 -o count | sort -k 1 | bedtools groupby -g 1 -c 2 -o sum > $pipeline_path/output/merged/eRNA.f18.eSNP.gtexCaviarDisease.txt

### TODO 
############################################################################################################################################################

# ====================================
# bidirectional_trans 
# ====================================
#select only the columns you need (name, RPM, pair, RPM, dis)
#cut -f4,7,14,17,21 plus.minus.bed > cut.plus.minus.bed 
#cut -f4,7,14,17,21 minus.plus.bed > cut.minus.plus.bed

#TODO test this 
if [[ $STRAND == "minus" ]]; 
then
Rscript eRNA.bidir.R $EXTERNAL_FEATURE/bidir/cut.plus.minus.bed plus.minus
elif [[ $STRAND == "plus" ]]; 
then
Rscript eRNA.bidir.R $EXTERNAL_FEATURE/bidir/cut.minus.plus.bed minus.plus
fi
############################################################################################################################################################
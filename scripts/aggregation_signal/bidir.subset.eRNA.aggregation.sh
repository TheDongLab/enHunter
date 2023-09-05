#================================================================
# create the bigwig signals bins of smaller size than ROADMAP.subset.eRNA.aggregation.sh
# for bidirectional transcription peak analysis  

## modified from https://github.com/sterding/BRAINcode/blob/7b4c5e816ff1cf9af86041326b71cf3f3e2e4bf6/src/eRNA.aggregation.sh
#================================================================
AMPPD_eRNA=/data/bioinformatics/projects/donglab/AMPPD_eRNA
conda activate /PHShome/rw552/condaenvs/ucsc

# 500 bp upstream and downstream is enough 
intersectBed -a $AMPPD_eRNA/inputs/plus/eRNA_stranded.bed -b regions_enh_merged.blood.narrowPeak -wao | sort -k4,4 -k17,17nr | awk '{if($4!=id) {id=$4;print;}}' | awk '{OFS="\t"; mid=($7==".")?int(($2+$3)/2):($8+$16); print $1, mid-500, mid+500,$4, $5, $6;}'> eRNA.500bp.summit.or.mid.500bp.plus.bed
intersectBed -a $AMPPD_eRNA/inputs/minus/eRNA_stranded.bed -b regions_enh_merged.blood.narrowPeak -wao | sort -k4,4 -k17,17nr | awk '{if($4!=id) {id=$4;print;}}' | awk '{OFS="\t"; mid=($7==".")?int(($2+$3)/2):($8+$16); print $1, mid-500, mid+500,$4, $5, $6;}'> eRNA.500bp.summit.or.mid.500bp.minus.bed

#hg38 
# reduce window size to 20 bp 
bedtools makewindows -b eRNA.500bp.summit.or.mid.500bp.minus.bed -w 20 -i srcwinnum > eRNA.500bp.summit.or.mid.500bp.20bp.windows.minus.bed
bedtools makewindows -b eRNA.500bp.summit.or.mid.500bp.plus.bed -w 20 -i srcwinnum > eRNA.500bp.summit.or.mid.500bp.20bp.windows.plus.bed

## run on RNAseq, CAGE (total ctss), and CAGE PD + HC 

## RNAseq 
bigWigAverageOverBed /data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs/minus/combined.mean.normalized.random.samplesN200.minus.bigwig eRNA.500bp.summit.or.mid.500bp.20bp.windows.minus.bed combined.mean.normalized.random.samplesN200.minus.500bp.tab
bigWigAverageOverBed /data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs/plus/combined.mean.normalized.random.samplesN200.plus.bigwig eRNA.500bp.summit.or.mid.500bp.20bp.windows.plus.bed combined.mean.normalized.random.samplesN200.plus.500bp.tab

## total ctss 
CAGE=/data/bioinformatics/external_data/externalData/CAGE
bigWigAverageOverBed $CAGE/CAGE.FANTOM5.total.rev.bigwig eRNA.500bp.summit.or.mid.500bp.20bp.windows.minus.bed CAGE.FANTOM5.total.rev.minus.500bp.tab
bigWigAverageOverBed $CAGE/CAGE.FANTOM5.total.fwd.bigwig eRNA.500bp.summit.or.mid.500bp.20bp.windows.plus.bed CAGE.FANTOM5.total.fwd.500bp.tab

## CAGE PD + HC 
RIKEN=/data/bioinformatics/external_data/externalData/RIKEN/PD_blood_CTSS
bigWigAverageOverBed $RIKEN/combined.ctss.minus.bw eRNA.500bp.summit.or.mid.500bp.20bp.windows.minus.bed combined.ctss.minus.500bp.tab
bigWigAverageOverBed $RIKEN/combined.ctss.plus.bw eRNA.500bp.summit.or.mid.500bp.20bp.windows.plus.bed combined.ctss.plus.500bp.tab






 
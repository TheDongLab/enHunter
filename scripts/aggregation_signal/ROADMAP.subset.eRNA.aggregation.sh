#================================================================
# create the bigwig signals for making the eRNA aggregation plot
# extract the Histone and DNase blood bigwig files from ROADMAP 

## modified from https://github.com/sterding/BRAINcode/blob/7b4c5e816ff1cf9af86041326b71cf3f3e2e4bf6/src/eRNA.aggregation.sh
#================================================================

## PART 1 DNase

#############  1. download the data 

## ALL SAMPLES ##
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

## SELECTED ##
#E062 - Primary mononuclear cells from peripheral blood (PBMC) <- No DNase file 
#E034 - Primary T cells from peripheral blood (1.1 G)
#E029 - Primary monocytes from peripheral blood (821 M)
#E032 - Primary B cells from peripheral blood (1.2 G)
#E046 - Primary Natural Killer cells from peripheral blood (1.1 G)
cd /data/bioinformatics/external_data/externalData/DNase/ROADMAP

grep DNase macs2signal.list | cut -f2 | grep -E 'E062|E034|E029|E032|E046' | awk '{print "http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/"$1}' | xargs --max-args 1 --max-procs 3 wget -b
mkdir download; mv E* download/

conda activate /PHShome/rw552/condaenvs/ucsc 
#### make sure that ssl in version 1.0.0 
#### conda install -c anaconda openssl=1.0
## bigWigMerge adds together the signal values

bsub -q normal -n 1 -M 4000 -o DNaseMerge.output bigWigMerge download/E*DNase.pval.signal.bigwig merged.DNase.pval.signal.bg 

#############  2. liftover bw files from hg19 to hg38
# bw -> bedGraph -> liftover -> bw
# bigWigToBedGraph input.bw input.bedGraph
# liftOver input.bedGraph hg19ToHg38.over.chain input_hg38.bedgraph unMapped
# fetchChromSizes hg38 > hg38.chrom.sizes
# LC_COLLATE=C sort -k1,1 -k2,2n input_hg38.bedgraph > input_hg38.sorted.bedgraph
# bedGraphToBigWig input_hg38.sorted.bedgraph hg38.chrom.sizes output.bw

module load liftover/1.0 
# liftover from hg19 to hg38
bsub -q normal -o liftOver.output -n 1 -M 2000 liftOver merged.DNase.pval.signal.bg hg19ToHg38.over.chain.gz merged.DNase.pval.signal.hg38.bedgraph unMappedDNaseMerged

#check if file is sorted
# LC_COLLATE=C sort -k1,1 -k2,2n merged.DNase.pval.signal.hg38.bedgraph 
# i think this takes too long 
# find other solutions! https://www.biostars.org/p/66927/#66931 

# conda install -c bioconda bedops --freeze-installed
mkdir tmp
bsub -q normal -o sort.output -n 1 -M 2000 sort-bed --max-mem 2G --tmpdir tmp merged.DNase.pval.signal.hg38.bedgraph > merged.DNase.pval.signal.sorted.hg38.bedgraph
## produces a book-ended bedgraph 

#############  3. move bedgraph to bigwig
# conda install -c bioconda ucsc-bedgraphtobigwig --freeze-installed
# conda activate /PHShome/rw552/condaenvs/ucsc 

# ERROR: Overlapping regions in bedGraph line 2276 of merged.DNase.pval.signal.sorted.hg38.bedgraph #
# must merge overlapping regions in bedgraph... overlapping regions may have resulted from liftover 
# https://www.biostars.org/p/171443/ 

bsub -q normal -o mergebg.output -n 1 -M 1000 bedtools merge -d -1 -c 4 -o sum -i merged.DNase.pval.signal.sorted.hg38.bedgraph > final.merged.DNase.pval.signal.sorted.hg38.bedgraph

bsub -q normal -o bgtobw.output -n 1 -M 4000 bedGraphToBigWig final.merged.DNase.pval.signal.sorted.hg38.bedgraph /PHShome/rw552/Documents/hg38.chrom.sizes merged.DNase.pval.signal.hg38.bigwig

#############  4. find blood DNase enhancer regions and convert to narrowPeak
awk '{OFS="\t"; print $0,"peak"NR }' regions_enh_merged.hg38.blood.bed > named.regions_enh_merged.blood.bed
# split -l : split based on number of lines 
# prefix for new files is tmp.regions_enh_merged.blood.bed
mkdir tmp
split -l 1000 named.regions_enh_merged.blood.bed tmp/tmp.regions_enh_merged.blood.bed

# conda install -c bioconda ucsc-bigwigsummary --freeze-installed

AMPPD_eRNA=/data/bioinformatics/projects/donglab/AMPPD_eRNA
wd=/PHShome/rw552/Documents

# test 
# $AMPPD_eRNA/bin/bed2narrowpeak.sh tmp/tmp.regions_enh_merged.blood.bedab merged.DNase.pval.signal.hg38.bigwig

for i in tmp/tmp.regions_enh_merged.blood.bed*; do echo $i; bsub -q normal -n 1 -M 1000 bash $AMPPD_eRNA/bin/bed2narrowpeak.sh $i merged.DNase.pval.signal.hg38.bigwig; done # for each 1000 length, do bed2narrowpeak
cat tmp/tmp.regions_enh_merged.blood.bed*narrowPeak > regions_enh_merged.blood.narrowPeak
rm -r tmp

############# 5. get the bin signal 
# change to use the narrowPeak based on merged Roadmap blood regions
# finds the "DNase peaks" coordinates 
# if there is no over lap with the Roadmap enhancer region, take the midpoint of the TNE 
# otherwise, take the peak of the overlap region of DNase signal and TNE 

#hg38
# awk '{if($4!=id) {id=$4;print;}}' ### if one TNE interacts with two DNase regions, just pick the first one 
# sort -k4,4 -k17,17nr ### sort by name (col 4) and by numeric reverse order for col 17 (with most amount of overlap)

# awk '{OFS="\t"; mid=($7==".")?int(($2+$3)/2):($8+$16); print $1, mid-1000, mid+1000,$4;}'> eRNA.1kbp.summit.or.mid.1kbp.bed
### $7 is DNase narrowPeak file coordiate (if there is overlap, else it is .)
### $8 start coordinate of DNase narrowPeak file coordinate
### $16 distance from DNase narrowPeak start coordinate to DNase peak 

# check to see if any narrowPeak points could not be called 
### awk '$10 == -1' regions_enh_merged.blood.narrowPeak 
### output was 0 lines

##################################################
intersectBed -a $AMPPD_eRNA/inputs/plus/eRNA_stranded.bed -b regions_enh_merged.blood.narrowPeak -wao | sort -k4,4 -k17,17nr | awk '{if($4!=id) {id=$4;print;}}' | awk '{OFS="\t"; mid=($7==".")?int(($2+$3)/2):($8+$16); print $1, mid-1000, mid+1000,$4, $5, $6;}'> eRNA.1kbp.summit.or.mid.1kbp.plus.bed
intersectBed -a $AMPPD_eRNA/inputs/minus/eRNA_stranded.bed -b regions_enh_merged.blood.narrowPeak -wao | sort -k4,4 -k17,17nr | awk '{if($4!=id) {id=$4;print;}}' | awk '{OFS="\t"; mid=($7==".")?int(($2+$3)/2):($8+$16); print $1, mid-1000, mid+1000,$4, $5, $6;}'> eRNA.1kbp.summit.or.mid.1kbp.minus.bed

#intersectBed -a $AMPPD_eRNA/inputs/eRNA_stranded_sorted.bed -b regions_enh_merged.blood.narrowPeak -wao | sort -k4,4 -k17,17nr | awk '{if($4!=id) {id=$4;print;}}' | awk '{OFS="\t"; mid=($7==".")?int(($2+$3)/2):($8+$16); print $1, mid-1000, mid+1000,$4, $5, $6;}'> eRNA.1kbp.summit.or.mid.1kbp.stranded.bed
### BRAINCODE #intersectBed -a $AMPPD_eRNA/inputs/eRNA_stranded_sorted.bed -b regions_enh_merged.blood.narrowPeak -wao | sort -k4,4 -k15,15nr | awk '{if($4!=id) {id=$4;print;}}' | awk '{OFS="\t"; mid=($5==".")?int(($2+$3)/2):($6+$14); print $1, mid-1000, mid+1000,$4;}'> eRNA.1kbp.summit.or.mid.1kbp.bed

mv eRNA.1kbp.summit.or.mid.1kbp.*.bed $AMPPD_eRNA

#hg19
module load liftover/1.0
liftOver eRNA.1kbp.summit.or.mid.1kbp.minus.bed $wd/hg38ToHg19.over.chain.gz eRNA.1kbp.summit.or.mid.1kbp.hg19.minus.bed unMappedRNA.1kbp.summit.or.mid.minus.1kbp
liftOver eRNA.1kbp.summit.or.mid.1kbp.plus.bed $wd/hg38ToHg19.over.chain.gz eRNA.1kbp.summit.or.mid.1kbp.hg19.plus.bed unMappedRNA.1kbp.summit.or.mid.plus.1kbp
#liftOver eRNA.1kbp.summit.or.mid.1kbp.bed $wd/hg38ToHg19.over.chain.gz eRNA.1kbp.summit.or.mid.1kbp.hg19.bed unMappedRNA.1kbp.summit.or.mid.1kbp

#######################################################################
## PART 2: download H2K4me1, H3K4me3, H3K27ac files from ROADMAP (hg19) 
#######################################################################

## lifting over DNase to hg38 was too difficult... 
## choosing to liftover hg38 blood TNE coordinates to hg19 instead

### not used
#module load liftover/1.0
#liftOver $AMPPD_eRNA/eRNA_stranded_sorted.bed $wd/hg38ToHg19.over.chain.gz $AMPPD_eRNA/eRNA_stranded_hg19.bed $AMPPD_eRNA/unMappedBloodTNEs

for t in E062 E034 E029 E032 E046; do
for i in H3K4me1 H3K4me3 H3K27ac; do
  echo http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/${t}-$i.pval.signal.bigwig >> download.url
done
done

conda activate downloading
conda install -c conda-forge axel
#conda install -c conda-forge parallel

#parallel -a download.url axel -n 5 ## can not create multiple threads?

while read p; do
  bsub -q normal -n 1 -M 1000 axel -n 5 $p
done < <(tail -n+2 download.url)

conda activate /PHShome/rw552/condaenvs/ucsc 

#bigwig merge
for i in H3K4me1 H3K4me3 H3K27ac; do
  bsub -q normal -n 1 -M 4000 bigWigMerge *${i}.pval.signal.bigwig RoadmapBlood${i}Merged.bg
done

hg19=/data/bioinformatics/referenceGenome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta

#sort bedgraph files before converting to bigwig
#this step is needed for newer versions of bedgraph to bigwig 
for i in H3K4me3 H3K27ac H3K4me1; do
  bsub -q normal -n 1 -M 2000 "sort-bed --max-mem 2G --tmpdir /data/bioinformatics/external_data/externalData/Histone/tmp RoadmapBlood${i}Merged.bg > RoadmapBlood${i}Merged.sorted.bg"
done

#bedgraph to bigwig
for i in H3K4me1 H3K4me3 H3K27ac; do
  bsub -q normal -n 1 -M 4000 bedGraphToBigWig RoadmapBlood${i}Merged.sorted.bg $hg19/hg19.chrom.size RoadmapBlood.hg19.${i}Merged.bigwig
done

rm *pval.signal.bigwig *bg

#cd ..
#for i in H3K4me1 H3K4me3 H3K27ac; do ln -fs Roadmap/RoadmapBlood${i}Merged.bigwig RoadmapBlood${i}Merged.bigwig; done

#######################################################################
################### PART 3: split by bins ############################
#######################################################################
# histone marks are in hg19

cd $AMPPD_eRNA
externalData=/data/bioinformatics/external_data/externalData

# split each TNE region into 100 bp windows 
# hg19 
bedtools makewindows -b eRNA.1kbp.summit.or.mid.1kbp.hg19.minus.bed -w 100 -i srcwinnum > eRNA.1kbp.summit.or.mid.1kbp.100bp.windows.minus.hg19.bed
bedtools makewindows -b eRNA.1kbp.summit.or.mid.1kbp.hg19.plus.bed -w 100 -i srcwinnum > eRNA.1kbp.summit.or.mid.1kbp.100bp.windows.plus.hg19.bed

#hg38 
bedtools makewindows -b eRNA.1kbp.summit.or.mid.1kbp.minus.bed -w 100 -i srcwinnum > eRNA.1kbp.summit.or.mid.1kbp.100bp.windows.minus.bed
bedtools makewindows -b eRNA.1kbp.summit.or.mid.1kbp.plus.bed -w 100 -i srcwinnum > eRNA.1kbp.summit.or.mid.1kbp.100bp.windows.plus.bed

#bedtools makewindows -b eRNA.1kbp.summit.or.mid.1kbp.bed -w 100 -i srcwinnum > eRNA.1kbp.summit.or.mid.1kbp.100bp.windows.bed
#bedtools makewindows -b eRNA.1kbp.summit.or.mid.1kbp.hg19.bed -w 100 -i srcwinnum > eRNA.1kbp.summit.or.mid.1kbp.100bp.windows.hg19.bed

# hg19 
for i in $externalData/Histone/*.bigwig;
do
   bsub -q normal -n 1 -M 1000 "bigWigAverageOverBed $i eRNA.1kbp.summit.or.mid.1kbp.100bp.windows.minus.hg19.bed $i.minus.tab" 
   bsub -q normal -n 1 -M 1000 "bigWigAverageOverBed $i eRNA.1kbp.summit.or.mid.1kbp.100bp.windows.plus.hg19.bed $i.plus.tab" 
done

# RNAseq, DNase, TFBS, CAGE, PhyloP are in hg38
# CAGE 
for i in $externalData/CAGE/*.bigwig;
do
    bsub -q normal -n 1 -M 1000 "bigWigAverageOverBed $i eRNA.1kbp.summit.or.mid.1kbp.100bp.windows.minus.bed $i.minus.tab"
    bsub -q normal -n 1 -M 1000 "bigWigAverageOverBed $i eRNA.1kbp.summit.or.mid.1kbp.100bp.windows.plus.bed $i.plus.tab"
done

# DNase
DNase=$externalData/DNase/ROADMAP/merged.DNase.pval.signal.hg38.bigwig
bsub -q normal -n 1 -M 1000 "bigWigAverageOverBed $DNase eRNA.1kbp.summit.or.mid.1kbp.100bp.windows.minus.bed DNase.minus.tab"
bsub -q normal -n 1 -M 1000 "bigWigAverageOverBed $DNase eRNA.1kbp.summit.or.mid.1kbp.100bp.windows.plus.bed DNase.plus.tab"

# TFBS
TFBS=$externalData/TFBS/TFBS.ENCODE.all.count.bigwig
bigWigAverageOverBed $TFBS eRNA.1kbp.summit.or.mid.1kbp.100bp.windows.minus.bed TFBS.minus.tab
bigWigAverageOverBed $TFBS eRNA.1kbp.summit.or.mid.1kbp.100bp.windows.plus.bed TFBS.plus.tab

# PhyloP
PhyloP=$externalData/Conservation/hg38.phyloP100way.bw
bsub -q normal -n 1 -M 1000 "bigWigAverageOverBed $PhyloP eRNA.1kbp.summit.or.mid.1kbp.100bp.windows.minus.bed PhyloP.minus.tab"
bsub -q normal -n 1 -M 1000 "bigWigAverageOverBed $PhyloP eRNA.1kbp.summit.or.mid.1kbp.100bp.windows.plus.bed PhyloP.plus.tab"

# RNAseq 
for i in inputs/*/*.bigwig;
do
    bsub -q normal -n 1 -M 1000 "bigWigAverageOverBed $i eRNA.1kbp.summit.or.mid.1kbp.100bp.windows.minus.bed $i.minus.tab"
    bsub -q normal -n 1 -M 1000 "bigWigAverageOverBed $i eRNA.1kbp.summit.or.mid.1kbp.100bp.windows.plus.bed $i.plus.tab"
done


############### PD vs Healthy Control #######################

# RNA seq
bsub -q normal -n 1 -M 1000 "bigWigAverageOverBed inputs/minus/combined.mean.normalized.control.samplesN200.minus.bigwig eRNA.1kbp.summit.or.mid.1kbp.100bp.windows.minus.bed combined.mean.normalized.control.samplesN200.minus.bigwig.minus.tab"
bsub -q normal -n 1 -M 1000 "bigWigAverageOverBed inputs/plus/combined.mean.normalized.control.samplesN200.plus.bigwig eRNA.1kbp.summit.or.mid.1kbp.100bp.windows.plus.bed combined.mean.normalized.control.samplesN200.plus.bigwig.plus.tab"

bsub -q normal -n 1 -M 1000 "bigWigAverageOverBed inputs/minus/combined.mean.normalized.PD.samplesN200.minus.bigwig eRNA.1kbp.summit.or.mid.1kbp.100bp.windows.minus.bed combined.mean.normalized.PD.samplesN200.minus.bigwig.minus.tab"
bsub -q normal -n 1 -M 1000 "bigWigAverageOverBed inputs/plus/combined.mean.normalized.PD.samplesN200.plus.bigwig eRNA.1kbp.summit.or.mid.1kbp.100bp.windows.plus.bed combined.mean.normalized.PD.samplesN200.plus.bigwig.plus.tab"

# CAGE
RIKEN=/data/bioinformatics/external_data/externalData/RIKEN/PD_blood_CTSS

conda activate /PHShome/rw552/condaenvs/ucsc 

# healthy control
sort -k 1,1 -k2,2n Ct.ctss.bed > Ct.ctss.sorted.bed
bedtools merge -i Ct.ctss.sorted.bed -d -1 -o sum -c 5 -S + > Ct.ctss.plus.bedgraph
bedtools merge -i Ct.ctss.sorted.bed -d -1 -o sum -c 5 -S - > Ct.ctss.minus.bedgraph

bedGraphToBigWig Ct.ctss.plus.bedgraph /PHShome/rw552/Documents/hg38.chrom.sizes Ct.ctss.plus.bw
bedGraphToBigWig Ct.ctss.minus.bedgraph /PHShome/rw552/Documents/hg38.chrom.sizes Ct.ctss.minus.bw

bigWigAverageOverBed $RIKEN/Ct.ctss.minus.bw eRNA.1kbp.summit.or.mid.1kbp.100bp.windows.minus.bed Ct.ctss.minus.tab
bigWigAverageOverBed $RIKEN/Ct.ctss.plus.bw eRNA.1kbp.summit.or.mid.1kbp.100bp.windows.plus.bed Ct.ctss.plus.tab

# PD
bsub -q normal -n 1 -M 1000 "sort -k 1,1 -k2,2n PD.ctss.bed > PD.ctss.sorted.bed"
bedtools merge -i PD.ctss.sorted.bed -d -1 -o sum -c 5 -S + > PD.ctss.plus.bedgraph
bedtools merge -i PD.ctss.sorted.bed -d -1 -o sum -c 5 -S - > PD.ctss.minus.bedgraph

bedGraphToBigWig PD.ctss.plus.bedgraph /PHShome/rw552/Documents/hg38.chrom.sizes PD.ctss.plus.bw
bedGraphToBigWig PD.ctss.minus.bedgraph /PHShome/rw552/Documents/hg38.chrom.sizes PD.ctss.minus.bw

bigWigAverageOverBed $RIKEN/PD.ctss.minus.bw eRNA.1kbp.summit.or.mid.1kbp.100bp.windows.minus.bed PD.ctss.minus.tab
bigWigAverageOverBed $RIKEN/PD.ctss.plus.bw eRNA.1kbp.summit.or.mid.1kbp.100bp.windows.plus.bed PD.ctss.plus.tab

##### combined CAGE signal for PD and Healthy Control 
conda activate /PHShome/rw552/condaenvs/ucsc 

bigWigMerge Ct.ctss.minus.bw PD.ctss.minus.bw combined.ctss.minus.bedGraph
sort -k1,1 -k2,2n combined.ctss.minus.bedGraph > combined.ctss.minus.sorted.bedGraph

bigWigMerge Ct.ctss.plus.bw PD.ctss.plus.bw combined.ctss.plus.bedGraph
sort -k1,1 -k2,2n combined.ctss.plus.bedGraph > combined.ctss.plus.sorted.bedGraph

bedGraphToBigWig combined.ctss.minus.sorted.bedGraph /PHShome/rw552/Documents/hg38.chrom.sizes combined.ctss.minus.bw
bedGraphToBigWig combined.ctss.plus.sorted.bedGraph /PHShome/rw552/Documents/hg38.chrom.sizes combined.ctss.plus.bw

bigWigAverageOverBed $RIKEN/combined.ctss.minus.bw eRNA.1kbp.summit.or.mid.1kbp.100bp.windows.minus.bed combined.ctss.minus.tab
bigWigAverageOverBed $RIKEN/combined.ctss.plus.bw eRNA.1kbp.summit.or.mid.1kbp.100bp.windows.plus.bed combined.ctss.plus.tab

#### CODE DUMP ####

#[ -e $i.eRNA.1kbp.summit.or.mid.1kbp.100bins ] || bsub -q normal -n 1 -M 1000 "bin/toBinRegionsOnBigwig.sh $i eRNA.1kbp.summit.or.mid.1kbp.hg19.bed 100 > $i.eRNA.1kbp.summit.or.mid.1kbp.hg19.100bins"
#[ -e $i.eRNA.1kbp.summit.or.mid.1kbp.100bins ] || bsub -q normal -n 1 -M 1000 "bin/toBinRegionsOnBigwig.sh $i eRNA.1kbp.summit.or.mid.1kbp.bed 100 > $i.eRNA.1kbp.summit.or.mid.1kbp.hg19.100bins"
#[ -e $DNase.eRNA.1kbp.summit.or.mid.1kbp.100bins ] || bsub -q normal -n 1 -M 1000 "bin/toBinRegionsOnBigwig.sh $DNase eRNA.1kbp.summit.or.mid.1kbp.bed 100 > $DNase.eRNA.1kbp.summit.or.mid.1kbp.100bins"
#[ -e $i.eRNA.1kbp.summit.or.mid.1kbp.100bins ] || bsub -q normal -n 1 -M 1000 "toBinRegionsOnBigwig.sh $i eRNA.1kbp.summit.or.mid.1kbp.bed 100 > $i.eRNA.1kbp.summit.or.mid.1kbp.100bins"
#[ -e $TFBS.eRNA.1kbp.summit.or.mid.1kbp.100bins ] || bsub -q normal -n 1 -M 1000 "bin/toBinRegionsOnBigwig.sh $TFBS eRNA.1kbp.summit.or.mid.1kbp.bed 100 > $TFBS.eRNA.1kbp.summit.or.mid.1kbp.100bins"
#[ -e $PhyloP.eRNA.1kbp.summit.or.mid.1kbp.100bins ] || bsub -q normal -n 1 -M 1000 "bin/toBinRegionsOnBigwig.sh $PhyloP eRNA.1kbp.summit.or.mid.1kbp.bed 100 > $PhyloP.eRNA.1kbp.summit.or.mid.1kbp.100bins"

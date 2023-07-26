## script to generate data to draw aggregation plot for HITNE
## Note: update on Oct 21, 2022 (see old archieve in github if needed)

#============================================================
# extract DNase peaks from Roadmap brain samples
#============================================================

# download DNase files for blood samples 
# DNase files were originally downloaded from ROADMAP in braincode, but that is hg19 

# options:
# 1. liftover bw files from hg19 to hg38
# bw -> bedGraph -> liftover -> bw
# bigWigToBedGraph input.bw input.bedGraph
# liftOver input.bedGraph hg19ToHg38.over.chain input_hg38.bedgraph unMapped
# fetchChromSizes hg38 > hg38.chrom.sizes
# LC_COLLATE=C sort -k1,1 -k2,2n input_hg38.bedgraph > input_hg38.sorted.bedgraph
# bedGraphToBigWig input_hg38.sorted.bedgraph hg38.chrom.sizes output.bw

# 2. liftover TNE files from hg38 to hg19 

# 3. download blood DNAse bw files from ENCODE hg38 

#### DOWNLOAD FROM ENCODE (hg38) ####
# 1. http://hgdownload.soe.ucsc.edu/gbdb/hg38/bbi/wgEncodeRegDnase/wgEncodeRegDnaseUwHl60Signal.bw 
# 2. 


#### DOWNLOAD FROM ROADMAP (hg19) ####
## merge the uniform signal of individual samples into the final "DNase bigwig track"
#grep DNase macs2signal.list | cut -f2 | awk '{print "http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/"$1}' | xargs -n 1 -P 8 wget -b
#mkdir download; mv E* download/
#bsub -q big -n 1 -M 10000 bigWigMerge download/E*DNase.pval.signal.bigwig merged.DNase.pval.signal.bg
#bsub -q big -n 1 -M 10000 bedGraphToBigWig merged.DNase.pval.signal.bg $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size merged.DNase.pval.signal.bigwig
# scp to http://panda.partners.org/~xd010/tracks/

## define DNase peaks 
## using the 638304 enhancer-being DNase peaks defined in 22 Roadmap "blood" samples (http://egg2.wustl.edu/roadmap/web_portal/DNase_reg.html#delieation) [SELECTED as final solution]
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
# download regions_enh_*.bed.gz for 22 blood samples (E062 E034 E045 E033 E044 E043 E039 E041 E042 E040 E037 E048 E038 E047 E029 E031 E035 E051 E050 E036 E032 E046)

### RW: does this step only generate narrowPeak files for the merged enhancer regions??? (enhancer region defined by Roadmap)

# downloaded enhancer regions from https://egg2.wustl.edu/roadmap/data/byDataType/dnase/BED_files_enh/ 
# see /data/bioinformatics/external_data/externalData/DNase/ROADMAP

awk '{OFS="\t"; print $0,"peak"NR }' regions_enh_merged.hg38.blood.bed > named.regions_enh_merged.blood.bed
split -l 1000 named.regions_enh_merged.blood.bed tmp.regions_enh_merged.blood.bed # split enhancer regions into 1000 length ??

for i in tmp.regions_enh_merged.blood.bed*; do echo $i; bsub -q short -n 1 -M 1000 bash ~/pipeline/bin/bed2narrowpeak.sh $i merged.DNase.pval.signal.bigwig; done # for each 1000 length, do bed2narrowpeak
cat tmp.regions_enh_merged.blood.bed*narrowPeak > regions_enh_merged.blood.narrowPeak
rm tmp.regions_enh_merged.blood.bed*

#============================================================
# get bin signal
#============================================================
# change to use the narrowPeak based on merged Roadmap blood regions
# finds the "DNase peaks" coordinates 
# if there is no over lap with the Roadmap enhancer region, take the midpoint of the TNE 
# otherwise, take the peak of the overlap region of DNase signal and TNE 
intersectBed -a eRNA.bed -b ../externalData/DNase/regions_enh_merged.blood.narrowPeak -wao | sort -k4,4 -k15,15nr | awk '{if($4!=id) {id=$4;print;}}' | awk '{OFS="\t"; mid=($5==".")?int(($2+$3)/2):($6+$14); print $1, mid-1000, mid+1000,$4;}'> eRNA.1kbp.summit.or.mid.1kbp.bed
intersectBed -a eRNA.bed -b ../externalData/DNase/regions_enh_merged.blood.narrowPeak -wao | sort -k4,4 -k15,15nr | awk '{if($4!=id) {id=$4;print;}}' | awk '{OFS="\t"; mid=($5==".")?int(($2+$3)/2):($6+$14); print $1, mid-2000, mid+2000,$4;}'> eRNA.2kbp.summit.or.mid.2kbp.bed
 
# eRNA.bed = TNE regions defined by algorithm 
# do these TNE regions overlap with narrowPeak bigwig signals defined by Roadmap?

## QUESTION: isn't this a biased view? only takes the DNase signal of ROADMAP enhancer region supported TNEs -> nvm, it also takes the center point of HiTNEs

## below code is similar as bin/bed2narrowPeak.sh + above line
# cat eRNA.bed | while read chr start end name score strand rest
# do
#     N=`expr $end - $start`;
#     s=`bigWigSummary externalData/DNase/DNase.Roadmap.fBrain.pval.bigwig -udcDir=/tmp -type=mean $chr $start $end $N 2>/dev/null | sed 's/n\/a/0/g'`
#     [[ $s == "" ]] && whichmax=`expr $N / 2`
#     whichmax=`echo $s | awk 'BEGIN{max=-9999;}{for(i=1;i<=NF;i++) if($i>max) {which=i;max=$i}}END{print which}'`
#     [[ $name == "" ]] && name="$chr_$start_$end";
#     echo $chr `expr $start + $whichmax - 1000` `expr $start + $whichmax + 1000` $name;
# done | sed 's/ /\t/g' > eRNA.1kbp.summit.or.mid.1kbp.bed

################### make all of the histone bigwig files ####################

## OPTION 1: download H2K4me1, H3K4me3, H3K27ac files from ENCODE 
## cell lines: 
## GM12878 (Epstein-Barr virus transformed lymphoblastoid cell line)
## K562 (immortalized myelogenous leukemia cell line)

# make urls
for t in Gm12878 K562; do
for i in H3K4me1 H3K4me3 H3K27ac; do
  echo http://hgdownload.soe.ucsc.edu/hg38/bbi/wgEncodeReg/wgEncodeRegMarkH3k4me3/wgEncodeBroadHistone${t}${i}StdSig.bigWig >> download.url
done
done

parallel -a download.url axel -n 5

for i in H3K4me1 H3K4me3 H3K27ac; do
  bsub -q big -n 1 -M 10000 bigWigMerge *${i}StdSig.bigWig EncodeBlood${i}Merged.bg
done

#bedgraph to bigwig
# produces the final product used in the agg plot 
for i in H3K4me1 H3K4me3 H3K27ac; do
  bsub -q big -n 1 -M 10000 bedGraphToBigWig EncodeBlood${i}Merged.bg $GENOME/Sequence/WholeGenomeFasta/hg38.chrom.size EncodeBlood${i}Merged.bigwig
done

rm *StdSig.bigWig.bigWig *bg

# split by bins
for i in ../externalData/*/*.bigwig;
do
    [ -e $i.eRNA.1kbp.summit.or.mid.1kbp.100bins ] || bsub -q short -n 1 -M 1000 "toBinRegionsOnBigwig.sh $i eRNA.1kbp.summit.or.mid.1kbp.bed 100 > $i.eRNA.1kbp.summit.or.mid.1kbp.100bins"
done

## OPTION 2: download H2K4me1, H3K4me3, H3K27ac files from ROADMAP (hg19) 

# copied from /data/bioinformatics/external_data/externalData/_copy_from_brincode/externalData/Histone/Roadmap/readme.txt
for t in E067 E068 E069 E070 E071 E072 E073 E074 E081 E082; do
for i in H3K4me1 H3K4me3 H3K27ac; do
  echo http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/${t}-$i.pval.signal.bigwig >> download.url
done
done

parallel -a download.url axel -n 5

#bigwig merge
for i in H3K4me1 H3K4me3 H3K27ac; do
  bsub -q big -n 1 -M 10000 bigWigMerge *${i}.pval.signal.bigwig RoadmapBrain${i}Merged.bg
done

#bedgraph to bigwig
for i in H3K4me1 H3K4me3 H3K27ac; do
  bsub -q big -n 1 -M 10000 bedGraphToBigWig RoadmapBrain${i}Merged.bg $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size RoadmapBrain${i}Merged.bigwig
done

rm *pval.signal.bigwig *bg
cd ..
for i in H3K4me1 H3K4me3 H3K27ac; do ln -fs Roadmap/RoadmapBrain${i}Merged.bigwig RoadmapBrain${i}Merged.bigwig; done



# gets the bigwig signals of the bin regions previously calculated 
#for all the bigwig files (histone, DNase ... etc )
for i in ../externalData/*/*.bigwig;
do
    [ -e $i.eRNA.1kbp.summit.or.mid.1kbp.100bins ] || bsub -q short -n 1 -M 1000 "toBinRegionsOnBigwig.sh $i eRNA.1kbp.summit.or.mid.1kbp.bed 100 > $i.eRNA.1kbp.summit.or.mid.1kbp.100bins"
done

#TODO how does he split the bigwig file by chromosome? I have no idea 
# something like this? awk '{print $0 >> "eRNA.2kbp.summit.or.mid.2kbp.bed."$1}' eRNA.2kbp.summit.or.mid.2kbp.bed

# for externalData/Histone/Roadmap -- too large to run bigWigSummmary as a whole
for i in ../externalData/Histone/Roadmap/*chr*bigwig;
do
    j=${i/*Merged./}
    chr=${j/.bigwig/}
    [ -e $i.eRNA.1kbp.summit.or.mid.1kbp.100bins.$chr ] || bsub -q short -n 1 -M 500 "toBinRegionsOnBigwig.sh $i eRNA.1kbp.summit.or.mid.1kbp.bed.$chr 100 > $i.eRNA.1kbp.summit.or.mid.1kbp.100bins.$chr"
done

# combine the 100bins files for each bin 
for i in H3K4me1 H3K4me3 H3K27ac; do cat ../externalData/Histone/Roadmap/RoadmapBrain${i}Merged.chr*eRNA.1kbp.summit.or.mid.1kbp.100bins.chr* > ../externalData/Histone/$i.Roadmap.brainMerged.bigwig.eRNA.1kbp.summit.or.mid.1kbp.100bins; done
rm ../externalData/Histone/Roadmap/RoadmapBrain*Merged.chr*eRNA.1kbp.summit.or.mid.1kbp.100bins.chr*

# I don't think 2kbp windows are actually being used so this code isn't useful 
# awk '{print $0 >> "eRNA.2kbp.summit.or.mid.2kbp.bed."$1}' eRNA.2kbp.summit.or.mid.2kbp.bed
# for i in ../externalData/Histone/Roadmap/*chr*bigwig;
# do
#     j=${i/*Merged./}
#     chr=${j/.bigwig/}
#     [ -e $i.eRNA.2kbp.summit.or.mid.2kbp.100bins.$chr ] || bsub -q short -n 1 -M 500 "toBinRegionsOnBigwig.sh $i eRNA.2kbp.summit.or.mid.2kbp.bed.$chr 100 > $i.eRNA.2kbp.summit.or.mid.2kbp.100bins.$chr"
# done
# for i in H3K4me1 H3K4me3 H3K27ac; do cat ../externalData/Histone/Roadmap/RoadmapBrain${i}Merged.chr*eRNA.2kbp.summit.or.mid.2kbp.100bins.chr* > ../externalData/Histone/$i.Roadmap.brainMerged2k.bigwig.eRNA.1kbp.summit.or.mid.1kbp.100bins; done
# rm ../externalData/Histone/Roadmap/RoadmapBrain*Merged.chr*eRNA.2kbp.summit.or.mid.2kbp.100bins.chr*

#============================================================
# draw aggregation plot (to run in R console)
#============================================================
Rscript ~/neurogen/pipeline/RNAseq/src/eRNA.aggPlot.R

H3K4me3.Roadmap.brainMerged.bigwig -> Roadmap/RoadmapBrainH3K4me3Merged.bigwig

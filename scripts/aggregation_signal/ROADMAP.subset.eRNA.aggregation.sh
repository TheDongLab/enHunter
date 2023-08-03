#================================================================
# create the bigwig signals for making the eRNA aggregation plot
# extract the Histone and DNase blood bigwig files from ROADMAP 
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
# sort -c merged.DNase.pval.signal.hg38.bedgraph
### TODO i think this takes too long 
# find other solutions! https://www.biostars.org/p/66927/#66931 

bsub -q normal -o sort.output -n 1 -M 2000 LC_COLLATE=C sort -k1,1 -k2,2n merged.DNase.pval.signal.hg38.bedgraph > merged.DNase.pval.signal.sorted.hg38.bedgraph

#############  3. move bedgraph to bigwig
# conda install -c bioconda ucsc-bedgraphtobigwig --freeze-installed
bsub -q normal -o bgtobw.output -n 1 -M 4000 bedGraphToBigWig <(LC_COLLATE=C sort -k1,1 -k2,2n merged.DNase.pval.signal.hg38.bedgraph) /PHShome/rw552/Documents/hg38.chrom.sizes merged.DNase.pval.signal.hg38.bigwig

#############  4. find blood DNase enhancer regions and convert to narrowPeak
awk '{OFS="\t"; print $0,"peak"NR }' regions_enh_merged.hg38.blood.bed > named.regions_enh_merged.blood.bed
split -l 1000 named.regions_enh_merged.blood.bed tmp.regions_enh_merged.blood.bed # split enhancer regions into 1000 length ??

for i in tmp.regions_enh_merged.blood.bed*; do echo $i; bsub -q short -n 1 -M 1000 bash ~/pipeline/bin/bed2narrowpeak.sh $i merged.DNase.pval.signal.bigwig; done # for each 1000 length, do bed2narrowpeak
cat tmp.regions_enh_merged.blood.bed*narrowPeak > regions_enh_merged.blood.narrowPeak
rm tmp.regions_enh_merged.blood.bed*

############# 5. get the bin signal 
# change to use the narrowPeak based on merged Roadmap blood regions
# finds the "DNase peaks" coordinates 
# if there is no over lap with the Roadmap enhancer region, take the midpoint of the TNE 
# otherwise, take the peak of the overlap region of DNase signal and TNE 
intersectBed -a eRNA.bed -b ../externalData/DNase/regions_enh_merged.blood.narrowPeak -wao | sort -k4,4 -k15,15nr | awk '{if($4!=id) {id=$4;print;}}' | awk '{OFS="\t"; mid=($5==".")?int(($2+$3)/2):($6+$14); print $1, mid-1000, mid+1000,$4;}'> eRNA.1kbp.summit.or.mid.1kbp.bed

## PART 2: download H2K4me1, H3K4me3, H3K27ac files from ROADMAP (hg19) 

for t in E062 E034 E029 E032 E046; do
for i in H3K4me1 H3K4me3 H3K27ac; do
  echo http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/${t}-$i.pval.signal.bigwig >> download.url
done
done

parallel -a download.url axel -n 5

#bigwig merge
for i in H3K4me1 H3K4me3 H3K27ac; do
  bsub -q big -n 1 -M 10000 bigWigMerge *${i}.pval.signal.bigwig RoadmapBlood${i}Merged.bg
done

#liftover bedgraph from hg19 to hg38
for i in H3K4me1 H3K4me3 H3K27ac; do
  bsub -q big -n 1 -M 10000 liftOver RoadmapBlood${i}Merged.bg hg19ToHg38.over.chain RoadmapBlood${i}Merged.hg38.bg unMapped
  LC_COLLATE=C sort -k1,1 -k2,2n RoadmapBlood${i}Merged.hg38.bg > RoadmapBlood${i}Merged.sorted.hg38.bedgraph
done 


#bedgraph to bigwig
for i in H3K4me1 H3K4me3 H3K27ac; do
  bsub -q big -n 1 -M 10000 bedGraphToBigWig RoadmapBrain${i}Merged.bg $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size RoadmapBrain${i}Merged.bigwig
done

rm *pval.signal.bigwig *bg
cd ..
for i in H3K4me1 H3K4me3 H3K27ac; do ln -fs Roadmap/RoadmapBrain${i}Merged.bigwig RoadmapBrain${i}Merged.bigwig; done

# split by bins
for i in ../externalData/*/*.bigwig;
do
    [ -e $i.eRNA.1kbp.summit.or.mid.1kbp.100bins ] || bsub -q short -n 1 -M 1000 "toBinRegionsOnBigwig.sh $i  100 > $i.eRNA.1kbp.summit.or.mid.1kbp.100bins"
done


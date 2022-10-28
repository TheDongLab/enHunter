#============================================================
# Script to download bigwig files for various marks
# Usage: bash $pipeline_path/src/get.regulatoryMarks.data.sh
#============================================================

ANNOTATION=$GENOME/Annotation/Genes
externalData=~/projects/PD/results/eRNA/externalData
[ -d $externalData ] || mkdir $externalData
cd $externalData

#========================================
# get bigwig and enhancer for characteristic marks
#========================================

echo "getting CAGE..."
# ----------------------------------
[ -d $externalData/CAGE ] || mkdir $externalData/CAGE
cd $externalData/CAGE

curl -s http://fantom.gsc.riken.jp/5/datahub/hg38/reads/ctssTotalCounts.fwd.bw > CAGE.FANTOM5.total.fwd.bigwig
curl -s http://fantom.gsc.riken.jp/5/datahub/hg38/reads/ctssTotalCounts.rev.bw > CAGE.FANTOM5.total.rev.bigwig

##### enhancers specifically expressed in blood (blood differentially expressed enhancers)
curl -s https://slidebase.binf.ku.dk/human_enhancers/presets/serve/blood > blood_differentially_expressed_enh_hg19.bed
awk '{OFS="\t"; print $1, $2, $3, $4, 0, $6, $7, $8, $9, $10, $11, $12}' blood_differentially_expressed_enh_hg19.bed | liftOver stdin hg19ToHg38.over.chain.gz blood_differentially_expressed_enh_hg38.bed unMappedHumanBlood 

##### for expressed enhancers in all facets
mkdir expressed_enh
cd expressed_enh
curl -s https://slidebase.binf.ku.dk/human_enhancers/presets/serve/facet_expressed_enhancers.tgz | gunzip -c | tar xvf -
rm CL:*
for file in UBERON:*.bed
do
    tissue="$(echo $file | cut -d"_" -f2)"
    val="$(echo $file | cut -d"_" -f3)"
    if [ $val != "expressed" ]
    then
        tissue="$(echo $file | cut -d"_" -f2,3 --output-delimiter='_')"
    fi
    echo $tissue
    liftOver $file ../hg19ToHg38.over.chain.gz stdout '$file + ".unMappedTissueEnh"' |awk -v tissue=$tissue '{OFS="\t"; print $1, $2, $3, $4, $5, tissue, $7, $8, $9, $10, $11, $12 }' $file >> all_enhancer_tissues_hg38.bed
done
#note: all unmapped coordinates can be found in unMappedTissueEnh
sort -k1,1 -k2,2n all_enhancer_tissues_hg38.bed  > all_enhancer_tissues_hg38_sorted.bed
# for counts file
sort -k 6 all_enhancer_tissues_hg38.bed | bedtools groupby -g 6 -c 6 -o count > all_enhancer_tissues_hg38.counts

##### for differentially expressed enhancers in all facets 
mkdir diff_expressed_enh
cd diff_expressed_enh
curl -s https://slidebase.binf.ku.dk/human_enhancers/presets/serve/facet_differentially_expressed_0.05.tgz | gunzip -c | tar xvf -
rm CL:*
for file in UBERON:*.bed
do
    tissue="$(echo $file | cut -d"_" -f2)"
    val="$(echo $file | cut -d"_" -f3)"
    if [ $val != "differentially" ]
    then
        tissue="$(echo $file | cut -d"_" -f2,3 --output-delimiter='_')"
    fi
    echo $tissue
    awk '{OFS="\t"; print $1, $2, $3, $4, 0, $6, $7, $8, $9, $10, $11, $12}' $file | liftOver stdin ../hg19ToHg38.over.chain.gz stdout unMappedTissueEnh |awk -v tissue=$tissue '{OFS="\t"; print $1, $2, $3, $4, $5, tissue, $7, $8, $9, $10, $11, $12 }' $file >> all_diff_exp_enhancer_tissues_hg38.bed
done
#note: all unmapped coordinates can be found in unMappedTissueEnh
sort -k1,1 -k2,2n all_diff_exp_enhancer_tissues_hg38.bed  > all_diff_exp_enhancer_tissues_hg38_sorted.bed
# for counts file
sort -k 6 all_diff_exp_enhancer_tissues_hg38.bed | bedtools groupby -g 6 -c 6 -o count > all_diff_exp_enhancer_tissues_hg38.counts


echo "getting TFBS HOT region"
# ---------------------------------
[ -d $externalData/TFBS ] || mkdir $externalData/TFBS
cd $externalData/TFBS

# clustered TFBS (count all Peaks for 161 transcription factors in 91 cell types)
# only per TF (e.g. if a TF occur in >1 cell types, it's only counted once)
curl -s http://hgdownload.soe.ucsc.edu/goldenPath/hg38/encRegTfbsClustered/encRegTfbsClusteredWithCells.hg38.bed.gz | gunzip > encRegTfbsClusteredWithCells.hg38.bed
sort -k1,1 encRegTfbsClusteredWithCells.hg38.bed | bedItemOverlapCount hg38 -chromSize=$ANNOTATION/ChromInfo.txt stdin | sort -k1,1 -k2,2n | sed 's/ /\t/g' > encRegTfbsClusteredWithCells.hg38.bg
bedGraphToBigWig encRegTfbsClusteredWithCells.hg38.bg $ANNOTATION/ChromInfo.txt TFBS.ENCODE.all.count.bigwig

echo "getting VISTA region"
# ---------------------------------
[ -d $externalData/VISTA ] || mkdir $externalData/VISTA
cd $externalData/VISTA
curl -s "https://enhancer.lbl.gov/cgi-bin/imagedb3.pl?show=1;page_size=100;search.gene=;page=1;search.result=yes;form=ext_search;search.status=Positives;search.org=Human;action=search;search.sequence=1" | sed 's/<[^>]\+>//g' | grep Human | sed 's/:/\t/g;s/-/\t/g;s/>Human|//g;s/ | /\t/g;s/positive /+/g;s/negative /-/g;s/ /_/g' | awk '{OFS="\t"; id=""; if($5=="+") for(i=6;i<=NF;i++) id=$i";"id; else id="NULL"; print $1,$2,$3,$4,NF-5,$5,id;}' > hg19_pos_regions.`date +%F`
cut -f1-4 hg19_pos_regions.`date +%F` | liftOver stdin ../hg19ToHg38.over.chain.gz hg38.pos_regions.`date +%F`.bed unMapped

echo "getting Conservation"
# ---------------------------------
[ -d $externalData/Conservation ] || mkdir $externalData/Conservation
cd $externalData/Conservation

# # phyloP
# # http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/
# rsync -avz --progress rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/phyloP46way/vertebrate ./
# zcat vertebrate/chr*.gz | gzip -c > vertebrate/phyloP46way.wigFix.gz
# rm vertebrate/chr*.gz
# wigToBigWig vertebrate/phyloP46way.wigFix.gz $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size Conservation.phyloP46way.bigwig

echo "getting DNase"
# ---------------------------------
[ -d $externalData/DNase ] || mkdir $externalData/DNase
cd $externalData/DNase

# # ref: http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=377332155_oXR9pqyaZLzzFHxgO4t0YzyQseGN&c=chr1&g=wgEncodeRegDnaseClusteredV2
# Note: V3 differs from V2 as it includes clusters having only 1 cell type contributing to the cluster (previously excluded). For us, V2 is fine.
curl -s http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegDnaseClustered/wgEncodeRegDnaseClusteredV2.bed.gz | gunzip | awk '$5>=500 && $4>=5' > wgEncodeRegDnaseClusteredV2.filtered.bed

# ENCODE Hela
curl -s http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/wgEncodeUwDnaseHelas3RawRep1.bigWig > DNase.ENCODE.Hela.bigwig

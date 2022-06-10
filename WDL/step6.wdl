version 1.0

workflow Step6Task {

   call eRNASigNif
}

# step6: calculate the significance of eRNA
task eRNASigNif {
    input {
        File eRNA5
        File list_bw_file
        Int genome_size
        String SAMPLE_GROUP
        File toExclude
    }

  command <<<
        while read sample filename
    do
        echo " - sample:" $sample;

        format=${filename##*.} # get extension
        
        if [[ $format =~ "bigwig|BIGWIG|bw|BW" ]];
        then
        bedtools shuffle -seed 123 -excl ~{toExclude} -noOverlapping -i ~{eRNA5} -g ~{genome_size} | awk -vOFS="\t" '$4=$1"_"$2"_"$3;' | bigWigAverageOverBed $filename stdin stdout | cut -f1,5 > $filename.~{SAMPLE_GROUP}.rdbg &
        bigWigAverageOverBed $filename ~{eRNA5} stdout | cut -f1,5 | sort -k1,1 -o $filename.~{SAMPLE_GROUP}.eRNA.meanRPM
        elif [[ $format =~ "bam|BAM|cram|CRAM" ]]; 
        then
        bedtools shuffle -seed 123 -excl ~{toExclude} -noOverlapping -i ~{eRNA5} -g ~{genome_size} | awk -vOFS="\t" '$4=$1"_"$2"_"$3;' | bedtools coverage -a stdin -b $filename -d | groupBy -g 4 -c 6 -o mean > $filename.~{SAMPLE_GROUP}.rdbg &
        bedtools coverage -a ~{eRNA5} -b $filename -d | groupBy -g 4 -c 6 -o mean | sort -k1,1 -o $filename.~{SAMPLE_GROUP}.eRNA.meanRPM
        fi
        
    done < ~{list_bw_file}
  
    RScript TNE_caller.consistency.R ~{SAMPLE_GROUP} ~{list_bw_file}
    # output are eRNA.tmp5.meanRPM.xls, eRNA.tmp5.pvalues.xls, and eRNA.tmp5.pvalues.adjusted.xls

    #3. Select TNE with adjusted p <= 0.05: 
    # bonferroni for the major groups
    awk '{OFS="\t"; split($1,a,"_"); if($1~/^chr/) {if($4<=0.05) print a[1],a[2],a[3],$1}}' eRNA.tmp5.pvalues.adjusted.xls | sortBed > eRNA.bonferroni.bed
    # FDR with method <=0.05
    awk '{OFS="\t"; split($1,a,"_"); if($1~/^chr/) {if($5<=0.05) print a[1],a[2],a[3],$1}}' eRNA.tmp5.pvalues.adjusted.xls | sortBed > eRNA.fdr.bed

    # =================
    # use bonferroni for major cell types and FDR for minor cell types.
    ln -fs eRNA.fdr.bed eRNA.bed
    # loci.txt required for fasteQTL
    awk 'BEGIN{OFS="\t"; print "id","chr","s1","s2";}{print $4,$1,$2,$3;}' eRNA.bed > eRNA.loci.txt  
    # meanRPM
    paste eRNA.tmp5.pvalues.adjusted.xls eRNA.tmp5.meanRPM.xls | awk 'NR ==1 || $5<=0.05' | cut -f6- > eRNA.meanRPM.xls

  >>>

  output { 

  }
}
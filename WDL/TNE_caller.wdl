version 1.0

# TODO ADD R scripts as input variable paths?
# TODO remove intermediate eRNA.tmp files? (not sure how)

workflow enHunter {
    input {
        String SAMPLE_GROUP 
        File inputBG
        Int genome_size 
        File toExclude
        File list_bw_file
        File splicing_site 
        Int length_min
    }

    call step0 { 
      input: 
      list_bw_file=list_bw_file,
      SAMPLE_GROUP=SAMPLE_GROUP, 
      genome_size=genome_size, 
      toExclude=toExclude, 
      inputBG=inputBG
    }

    call step1 { 
      input: inputBG=inputBG
    }

    call step2 { 
      input: 
      eRNA1=step1.eRNA1,
      Dsig=step0.Dsig
    }

    call step3 { 
      input:
      eRNA2=step2.eRNA2, 
      toExclude=toExclude
    }

    call step4 { 
      input: 
      eRNA3=step3.eRNA3,
      length_min=length_min
    }

    call step5 { 
      input:
      eRNA4=step4.eRNA4, 
      splicing_site=splicing_site
    }

    call step6 { 
      input: 
      eRNA5=step5.eRNA5, 
      list_bw_file=list_bw_file,
      genome_size=genome_size, 
      SAMPLE_GROUP=SAMPLE_GROUP, 
      toExclude=toExclude
    }

    output { 
      File bed = step6.bed
      File loci = step6.loci
      File meanRPM = step6.meanRPM
    }
}

# step0: measure transcriptional noise in background genomic regions
task step0 {
  input {
    File list_bw_file
    String SAMPLE_GROUP
    Int genome_size
    File toExclude
    File inputBG
  }
    
    command { 
      bedtools random -seed 3 -g ~{genome_size} -l 1 -n 1000000 | sortBed | intersectBed -a - -b ~{toExclude} -v -sorted | intersectBed -a ~{inputBG} -b - -sorted -u | cut -f4 > transcriptional.noise.rpm.txt
      RScript TNE_caller.fit.Tx.noise.R transcriptional.noise.rpm.txt 
      tail -n1 transcriptional.noise.rpm.pvalues.txt
    }

    output { 
      #stdout could contain more than Dsig (check)
      Int Dsig = read_int(stdout())
    }
}

#step1: any regions with value > baseLevel (i.e. average sequencing depth)
task step1 {
  input {
    File inputBG 
  }

  command <<<
     basalLevel=$(tail -n1 ~{inputBG} | cut -f2 -d'=' | cut -f1)
     awk -vmin=basalLevel '{OFS="\t"; if($4>=min) print $1,$2,$3,".",$4}' ~{inputBG} | mergeBed -c 5 -o max > eRNA.tmp1
  >>>
  
  output {
     File eRNA1 = "eRNA.tmp1"
  }
}

# step2: summit RPM >=Dsig (density with p<0.05)
task step2 {
  input {
    File eRNA1
    Int Dsig 
  }

  command <<<
    awk -vD=~{Dsig} '{OFS="\t"; if($4>=D) print $1,$2,$3,".",$4}' ~{eRNA1} | mergeBed -d 100 -c 5 -o max > eRNA.tmp2
  >>>

  output { 
    File eRNA2 = "eRNA.tmp2"
  }
}

# step3: located in non-generic regions (e.g. 500bp away from any annotated exons)
task step3 {
  input {
    File eRNA2
    File toExclude
  }

  command { 
    intersectBed -a ~{eRNA2} -b  ~{toExclude} -v > eRNA.tmp3
  }

  output { 
    File eRNA3 = "eRNA.tmp3" 
  }
}

# step4: length > $length_min
task step4 {
  input {
    File eRNA3
    Int length_min
  }

  command <<< 
    awk -v len_min=~{length_min} '{OFS="\t"; if(($3-$2)>len_min) print $1,$2,$3,$1"_"$2"_"$3}' ~{eRNA3} > eRNA.tmp4
  >>>

  output { 
    File eRNA4 = "eRNA.tmp4" 
  }
}

# step5: don't contain any splicing sites (donor or acceptor from trinity/cufflinks de novo assembly)
task step5 {
  input {
    File eRNA4
    File splicing_site
  }

  command { 
    intersectBed -a ~{eRNA4} -b ~{splicing_site} -v > eRNA.tmp5
  }

  output { 
    File eRNA5 = "eRNA.tmp5"  
  }
}

# step6: calculate the significance of eRNA
task step6 {
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
    # meanRPM (RPM = reads per million)
    paste eRNA.tmp5.pvalues.adjusted.xls eRNA.tmp5.meanRPM.xls | awk 'NR ==1 || $5<=0.05' | cut -f6- > eRNA.meanRPM.xls

  >>>

  output { 
    #final outputs: 
    # eRNA.bed, eRNA.loci.txt, eRNA.meanRPM.xls
    File bed = "eRNA.bed"
    File loci = "eRNA.loci.txt"
    File meanRPM = "eRNA.meanRPM.xls"
  }

}

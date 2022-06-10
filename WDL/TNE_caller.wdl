version 1.0

#TODO ADD R scripts as input variable paths?

import "step6.wdl" as Step6Task

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
      #TODO make inputBG a conditional variable 
      init_inputBG=inputBG
    }

    call step1 { 
      input: inputBG=step0.inputBG
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

    call Step6Task.eRNASigNif { 
      input: 
      eRNA5=step5.eRNA5

    }
}

# step0: measure transcriptional noise in background genomic regions
task step0 {
  input {
    File list_bw_file
    String SAMPLE_GROUP
    Int genome_size
    File toExclude
    File init_inputBG
  }
    
    command { 
      #TODO merged signal might already be implemented 
      #also this might run into issues bc you are calling a .sh script 
      inputBG=$(TNE_caller.combine_bigwig.v2.sh ~{list_bw_file} ~{SAMPLE_GROUP})

      bedtools random -seed 3 -g ~{genome_size} -l 1 -n 1000000 | sortBed | intersectBed -a - -b ~{toExclude} -v -sorted | intersectBed -a ~{inputBG} -b - -sorted -u | cut -f4 > transcriptional.noise.rpm.txt
      RScript TNE_caller.fit.Tx.noise.R transcriptional.noise.rpm.txt 

      tail -n1 transcriptional.noise.rpm.pvalues.txt
    }

    output { 
      File inputBG="trimmedmean.uniq.normalized.${SAMPLE_GROUP}.bedGraph"

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

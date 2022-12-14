version 1.0

# https://hub.docker.com/repository/docker/rwang429/enhunter-filter 

workflow enHunter {

    input {
        Int n_CPU
        Float GB_memory
        Int GB_disk

        Float GB_memory_step6
        Int GB_disk_step6

        File R_fit         # TNE_caller.fit.Tx.noise.R
        File R_consistency # TNE_caller.consistency.R

        File genome_size 
        String SAMPLE_GROUP 
        File inputBG
        File toExclude
        File list_bw_file
        File splicing_site 
        Int length_min
    }

    call step0 { 
      input: 
        n_CPU = n_CPU,
        GB_memory = GB_memory,
        GB_disk = GB_disk, 

        R_fit = R_fit,
        genome_size=genome_size, 
        toExclude=toExclude, 
        inputBG=inputBG
    }

    call step1 { 
      input: 
        n_CPU = n_CPU,
        GB_memory = GB_memory,
        GB_disk = GB_disk, 

        inputBG=inputBG
    }

    call step2 { 
      input: 
        n_CPU = n_CPU,
        GB_memory = GB_memory,
        GB_disk = GB_disk, 

        eRNA1=step1.eRNA1,
        Dsig=step0.Dsig
    }

    call step3 { 
      input:
        n_CPU = n_CPU,
        GB_memory = GB_memory,
        GB_disk = GB_disk, 

        eRNA2=step2.eRNA2, 
        toExclude=toExclude
    }

    call step4 { 
      input: 
        n_CPU = n_CPU,
        GB_memory = GB_memory,
        GB_disk = GB_disk, 

        eRNA3=step3.eRNA3,
        length_min=length_min
    }

    call step5 { 
      input:
        n_CPU = n_CPU,
        GB_memory = GB_memory,
        GB_disk = GB_disk, 
        eRNA4=step4.eRNA4, 
        splicing_site=splicing_site
    }

 Array[File] input_bigwig = read_lines(list_bw_file)

    scatter (file in input_bigwig) {

      call sampleNameParser {
        input: 
        samplePath=file,
        n_CPU = n_CPU,
        GB_memory = 1.5,
        GB_disk = 10
      }

      call step6A { 
          input: 
            inputFile=file, 
            sampleName=sampleNameParser.sampleName,
            genome_size=genome_size, 
            SAMPLE_GROUP=SAMPLE_GROUP, 
            toExclude=toExclude, 
            eRNA5=step5.eRNA5, 
            n_CPU = n_CPU,
            GB_memory = GB_memory,
            GB_disk = GB_disk
        }
    }


    call step6 { 
      input: 
        SampleNames=sampleNameParser.sampleName,
        RDBG=step6A.rdbg,
        meanRPM=step6A.meanRPM,
        n_CPU = n_CPU,
        GB_memory = GB_memory_step6,
        GB_disk = GB_disk_step6,
        R_consistency = R_consistency, 
        SAMPLE_GROUP=SAMPLE_GROUP
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
    Int n_CPU
    Float GB_memory
    Int GB_disk
    File R_fit 
    File genome_size
    File toExclude
    File inputBG
  }
    
    command { 
      bedtools random -seed 3 -g ~{genome_size} -l 1 -n 1000000 | sortBed | intersectBed -a - -b ~{toExclude} -v -sorted | intersectBed -a ~{inputBG} -b - -sorted -u | cut -f4 > transcriptional.noise.rpm.txt
      Rscript ~{R_fit} transcriptional.noise.rpm.txt 
      tail -n1 transcriptional.noise.rpm.pvalues.txt > Dsig.txt
    }

    output { 
      Float Dsig = read_float("Dsig.txt")
    }

    runtime {
        docker: "rwang429/enhunter-filter:latest"
        memory: "${GB_memory}GB"
        disks: "local-disk ${GB_disk} SSD"
        cpu: "${n_CPU}"
        preemptible: 1
        zones: "us-central1-a us-central1-b"
    }
}

#step1: any regions with value > baseLevel (i.e. average sequencing depth)
task step1 {
  input {
    Int n_CPU
    Float GB_memory
    Int GB_disk

    File inputBG 
  }

  command <<<
     basalLevel=$(tail -n1 ~{inputBG} | cut -f2 -d'=' | cut -f1)
     awk -vmin=$basalLevel '{OFS="\t"; if($4>=min) print $1,$2,$3,".",$4}' ~{inputBG} | mergeBed -c 5 -o max > eRNA.tmp1
  >>>
  
  output {
     File eRNA1 = "eRNA.tmp1"
  }

  runtime {
    docker: "rwang429/enhunter-filter:latest"
    memory: "${GB_memory}GB"
    disks: "local-disk ${GB_disk} SSD"
    cpu: "${n_CPU}"
    preemptible: 1
    zones: "us-central1-a us-central1-b"
  }
}

# step2: summit RPM >=Dsig (density with p<0.05)
task step2 {
  input {
    Int n_CPU
    Float GB_memory
    Int GB_disk

    File eRNA1
    Float Dsig 
  }

  command <<<
    awk -vD=~{Dsig} '{OFS="\t"; if($4>=D) print $1,$2,$3,".",$4}' ~{eRNA1} | mergeBed -d 100 -c 5 -o max > eRNA.tmp2
  >>>

  output { 
    File eRNA2 = "eRNA.tmp2"
  }

  runtime {
    docker: "rwang429/enhunter-filter:latest"
    memory: "${GB_memory}GB"
    disks: "local-disk ${GB_disk} SSD"
    cpu: "${n_CPU}"
    preemptible: 1
    zones: "us-central1-a us-central1-b"
  }
}

# step3: located in non-generic regions (e.g. 500bp away from any annotated exons)
task step3 {
  input {
    Int n_CPU
    Float GB_memory
    Int GB_disk


    File eRNA2
    File toExclude
  }

  command { 
    intersectBed -a ~{eRNA2} -b  ~{toExclude} -v > eRNA.tmp3
  }

  output { 
    File eRNA3 = "eRNA.tmp3" 
  }

  runtime {
    docker: "rwang429/enhunter-filter:latest"
    memory: "${GB_memory}GB"
    disks: "local-disk ${GB_disk} SSD"
    cpu: "${n_CPU}"
    preemptible: 1
    zones: "us-central1-a us-central1-b"
  }
}

# step4: length > $length_min
task step4 {
  input {
    Int n_CPU
    Float GB_memory
    Int GB_disk


    File eRNA3
    Int length_min
  }

  command <<< 
    awk -v len_min=~{length_min} '{OFS="\t"; if(($3-$2)>len_min) print $1,$2,$3,$1"_"$2"_"$3}' ~{eRNA3} > eRNA.tmp4
  >>>

  output { 
    File eRNA4 = "eRNA.tmp4" 
  }

  runtime {
    docker: "rwang429/enhunter-filter:latest"
    memory: "${GB_memory}GB"
    disks: "local-disk ${GB_disk} SSD"
    cpu: "${n_CPU}"
    preemptible: 1
    zones: "us-central1-a us-central1-b"
  }
}

# step5: don't contain any splicing sites (donor or acceptor from trinity/cufflinks de novo assembly)
task step5 {
  input {
    Int n_CPU
    Float GB_memory
    Int GB_disk


    File eRNA4
    File splicing_site
  }

  command { 
    intersectBed -a ~{eRNA4} -b ~{splicing_site} -v > eRNA.tmp5
  }

  output { 
    File eRNA5 = "eRNA.tmp5"  
  }


  runtime {
    docker: "rwang429/enhunter-filter:latest"
    memory: "${GB_memory}GB"
    disks: "local-disk ${GB_disk} SSD"
    cpu: "${n_CPU}"
    preemptible: 1
    zones: "us-central1-a us-central1-b"
  }
}

task sampleNameParser {
  input {
    String samplePath
    Int n_CPU
    Float GB_memory
    Int GB_disk
  }
 
  command <<<
      sample=$(basename ~{samplePath})
      sampleName=$(echo $sample | cut -d '.' -f 1)

      echo $sampleName > "sampleName.txt"
  >>>

  output {
    String sampleName = read_string("sampleName.txt")
  }

  runtime {
    docker: "rwang429/enhunter-filter:latest"
    memory: "~{GB_memory}GB"
    disks: "local-disk ~{GB_disk} SSD"
    cpu: "~{n_CPU}"
    preemptible: 1
    zones: "us-central1-a us-central1-b"
  }
}

task step6A {
  input {
    Int n_CPU
    Float GB_memory
    Int GB_disk

    File inputFile
    File toExclude 
    File eRNA5
    File genome_size
    String SAMPLE_GROUP
    String sampleName 
  }

  command <<<
    bedtools shuffle -seed 123 -excl ~{toExclude} -noOverlapping -i ~{eRNA5} -g ~{genome_size} | awk -vOFS="\t" '$4=$1"_"$2"_"$3;' | bigWigAverageOverBed ~{inputFile} stdin stdout | cut -f1,5 > ~{sampleName}.~{SAMPLE_GROUP}.rdbg
    bigWigAverageOverBed ~{inputFile} ~{eRNA5} stdout | cut -f1,5 | sort -k1,1 -o ~{sampleName}.~{SAMPLE_GROUP}.eRNA.meanRPM
  >>>

  output {
    File rdbg = "~{sampleName}.~{SAMPLE_GROUP}.rdbg"
    File meanRPM = "~{sampleName}.~{SAMPLE_GROUP}.eRNA.meanRPM"
  }

  runtime {
    docker: "rwang429/enhunter-filter:latest"
    memory: "~{GB_memory}GB"
    disks: "local-disk ~{GB_disk} SSD"
    cpu: "~{n_CPU}"
    preemptible: 1
    zones: "us-central1-a us-central1-b"
  }
}


# step6: calculate the significance of eRNA
task step6 {
    input {
      Array[String] SampleNames

      Array[File] RDBG
      Array[File] meanRPM
      Int n_CPU
      Float GB_memory
      Int GB_disk

      File R_consistency
      String SAMPLE_GROUP

    }

  File SampleIDs = write_lines(SampleNames)

  command <<<
    echo ~{sep=' ' meanRPM} | sed 's/ /\n/g' > sampleRPMs.txt
    echo ~{sep=' ' RDBG} | sed 's/ /\n/g'> sampleRDBGs.txt

    Rscript ~{R_consistency} ~{SAMPLE_GROUP} ~{SampleIDs} sampleRDBGs.txt sampleRPMs.txt

    awk '{OFS="\t"; split($1,a,"_"); if($1~/^chr/) {if($4<=0.05) print a[1],a[2],a[3],$1}}' eRNA.tmp5.pvalues.adjusted.xls | sortBed > eRNA.bonferroni.bed
    awk '{OFS="\t"; split($1,a,"_"); if($1~/^chr/) {if($5<=0.05) print a[1],a[2],a[3],$1}}' eRNA.tmp5.pvalues.adjusted.xls | sortBed > eRNA.fdr.bed

    ln -fs eRNA.fdr.bed eRNA.bed
    awk 'BEGIN{OFS="\t"; print "id","chr","s1","s2";}{print $4,$1,$2,$3;}' eRNA.bed > eRNA.loci.txt  
    paste eRNA.tmp5.pvalues.adjusted.xls eRNA.tmp5.meanRPM.xls | awk 'NR ==1 || $5<=0.05' | cut -f6- > eRNA.meanRPM.xls

  >>>

  output { 
    #final outputs: 
    # eRNA.bed, eRNA.loci.txt, eRNA.meanRPM.xls
    File bed = "eRNA.bed"
    File loci = "eRNA.loci.txt"
    File meanRPM = "eRNA.meanRPM.xls"
  }

  runtime {
    docker: "rwang429/enhunter-filter:latest"
    memory: "${GB_memory}GB"
    disks: "local-disk ${GB_disk} SSD"
    cpu: "${n_CPU}"
    preemptible: 1
    zones: "us-central1-a us-central1-b"
  }

}
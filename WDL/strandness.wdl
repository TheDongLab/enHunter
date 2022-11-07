version 1.0

## An investigation into the antisense mapping still prevalent in deepTools bamCoverage
## 1. localize AMP-PD bam file, subset to chr1:25553420-25555318 region. convert bam to sam 
## 2. grep for XS:A:+ and XS:A:- 

workflow strandedness_workflow {
  input {
        Int n_CPU
        Float GB_memory
        Int GB_disk

        File list_bam_files
    }
    
    Array[File] input_bams = read_lines(list_bam_files)
    
     scatter (file in input_bams) {
      call sampleNameParser {
        input: 
        inputFile=file,
        n_CPU = 1,
        GB_memory = 1,
        GB_disk = 5
      }

      call bam2sam {
        input: 
        sampleName=sampleNameParser.sampleName,
        inputFile=file,
        n_CPU = n_CPU,
        GB_memory = GB_memory,
        GB_disk = GB_disk
      }
    }
    
    call strand {
      input: 
        sam_files=bam2sam.sam, 
        n_CPU = n_CPU,
        GB_memory = GB_memory,
        GB_disk = GB_disk
    }
    
    output {
      File minus=strand.plus 
      File plus=strand.minus
    }
  
}

task sampleNameParser {
  input {
    File inputFile 
    Int n_CPU
    Float GB_memory
    Int GB_disk
  }
  String samplePath = inputFile

  command <<<
      sample=$(basename ~{samplePath})
      sampleName=$(echo $sample | cut -d '.' -f 1)

      echo $sampleName > "sampleName.txt"
  >>>

  output {
    String sampleName = read_string("sampleName.txt")
  }

  runtime {
    docker: "ubuntu:latest"
    memory: "~{GB_memory}GB"
    disks: "local-disk ~{GB_disk} SSD"
    cpu: "~{n_CPU}"
    preemptible: 1
    zones: "us-central1-a us-central1-b"
  }
}

task bam2sam {
  input {
        String sampleName 
        Int n_CPU
        Float GB_memory
        Int GB_disk
        
        File inputFile 
  }
    String samplePath = inputFile
      
  command <<<
      samtools view -o ~{sampleName}.sam -h ~{inputFile} chr1:25553450-25555318
  >>>
  
  output {
    File sam = "~{sampleName}.sam"
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

task strand {
  
  input {
    Array[File] sam_files 
    
    Int n_CPU
    Float GB_memory
    Int GB_disk
  }
  
  command {
    for x in ~{sep=' ' sam_files} 
    do 
  
      plus=$(grep -c XS:A:+ $x)
      minus=$(grep -c XS:A:- $x) 
      
      echo $x $plus '\n' >> plus_spliced_exons.txt
      echo $x $minus '\n' >> minus_spliced_exons.txt
      
    done 
  }
  
  output {
    File plus = "plus_spliced_exons.txt"
    File minus = "minus_spliced_exons.txt"
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
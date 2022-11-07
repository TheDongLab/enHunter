version 1.0

## An investigation into the antisense mapping still prevalent in deepTools bamCoverage
## 1. localize AMP-PD bam file, subset to chr1:25553420-25555318 region. convert bam to sam 
## 2. grep for XS:A:+ and XS:A:- 

workflow strandedness_workflow {
  input {
        Int n_CPU
        Float GB_memory
        Int GB_disk
        
        File BAM 
        File BAI
        String SampleID
        
    }

      call bam2sam {
        input: 
        BAM=BAM, 
        BAI=BAI, 
        SampleID=SampleID,
        n_CPU = n_CPU,
        GB_memory = GB_memory,
        GB_disk = GB_disk
      }
    
    output {
    File bam = bam2sam.bam
    File index = bam2sam.index
    File antisense = bam2sam.antisense 
    }
  
}

task bam2sam {
  input {
        Int n_CPU
        Float GB_memory
        Int GB_disk
        
        File BAM
        File BAI
        String SampleID
  }
      
  command <<<
      samtools view -o subset~{SampleID}.bam -b ~{BAM} chr1:25553450-25555318
      samtools index subset~{SampleID}.bam > subset~{SampleID}.bam.bai 
      samtools view -h -o subset~{SampleID}.sam subset~{SampleID}.bam
      grep XS:A:- subset~{SampleID}.sam > antisense.txt
  >>>
      
  output {
    File bam = "subset~{SampleID}.bam"
    File index = "subset~{SampleID}.bam.bai"
    File antisense = "antisense.txt"
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

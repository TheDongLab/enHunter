version 1.0

## WDL script to get the count of alignments in each input BAM file using featureCounts 

workflow getTranscriptReadCount {
      
      input {
        String sample_name
        File bams_list
        File TNEs
        Int n_CPU
        Float GB_memory
        Int GB_disk
      }
      
      Array[File] input_bams = read_lines(bams_list)
       
      call ReadCount {
        input: 
          sample_name = sample_name, 
          input_bams = input_bams, 
          TNE = TNEs, 
          n_CPU = n_CPU,
          GB_memory = GB_memory,
          GB_disk = GB_disk
      }
      
      output {
        File countTable = ReadCount.read_counts
        File summary = ReadCount.summary
      }

}

task ReadCount {
  input {
    String sample_name
    Array[File] input_bams
    File TNE
    Int n_CPU
    Float GB_memory
    Int GB_disk
  }
    
  command <<< 
    featureCounts -a ~{saf} -F SAF -f -s 1 -p --verbose -o ~{sample_name} ~{sep=' ' input_bams}
  >>>
  
  output {
    File read_counts = "~{sample_name}"
    File summary = "~{sample_name}.summary"
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
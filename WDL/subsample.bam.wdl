workflow subsample_bams { 
    input {
        File all_bam_files
    }

    Array[File] inputBams = read_lines(all_bam_files)

    scatter(oneFile in inputFiles) {
        call subsample_each_bam {
            input: 
                cur_bam = oneFile, 
                n_CPU = n_CPU,
                GB_memory = GB_memory,
                GB_disk = GB_disk
        }
    }

}

task subsample_each_bam {
    input {
        File cur_bam
    }

    command <<<
          for p in `seq 0.1 0.1 0.9`;
  do
         echo “== $p”
         samtools view -b -s $p -o $p.bam ~{cur_bam}
  done
    >>>

    output {
        Array[File] output_bam = glob("*.bam")
    }


    runtime {
        docker: "bwhbioinformaticshub/bioinformatics-basic:latest"
        memory: "${GB_memory}GB"
        disks: "local-disk ${GB_disk} SSD"
        cpu: "${n_CPU}"
        preemptible: 1
        zones: "us-central1-a us-central1-b"
    }
}
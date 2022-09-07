version 1.0 

workflow calculate_depth {

    input {
        Int n_CPU
        Float GB_memory
        Int GB_disk

        File bam_file_list
    }

    Array[File] bam_files = read_lines(bam_file_list)

    scatter(bam in bam_files) {
        call findDepth {
            input: 
            inputFile=bam,
            n_CPU = n_CPU,
            GB_memory = GB_memory,
            GB_disk = GB_disk

        }
    }

    call combineDepths {
    input:
        allDepths = findDepth.depth, 
        n_CPU = n_CPU,
        GB_memory = GB_memory,
        GB_disk = GB_disk
    }


}

task findDepth {
    input {
        File inputFile
        Int n_CPU
        Float GB_memory
        Int GB_disk
    }

    command <<<
    samtools view -c -F 2308 ~{inputFile} > numreads.txt
    >>>

    output {
        String depth = read_string("numreads.txt")
    }
    runtime {
        memory: "${GB_memory}GB"
        disks: "local-disk ${GB_disk} SSD"
        cpu: "${n_CPU}"
        docker: "bwhbioinformaticshub/bioinformatics-basic:latest"
    }
}

task combineDepths { 
    input {
        Array[String] allDepths
        Int n_CPU
        Float GB_memory
        Int GB_disk
    }

    command <<<

    >>>

    runtime {
        memory: "${GB_memory}GB"
        disks: "local-disk ${GB_disk} SSD"
        cpu: "${n_CPU}"
        docker: "bwhbioinformaticshub/bioinformatics-basic:latest"
    }

    output {
        File numReads = write_lines(allDepths)
    }
}

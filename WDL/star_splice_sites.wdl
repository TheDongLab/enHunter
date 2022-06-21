version 1.0

workflow Scatter {

    input {
        File bigwig_files_list
        Int n_CPU
        Float GB_memory
        Int GB_disk
    }

    Array[File] inputFiles = read_lines(bigwig_files_list)

    scatter(oneFile in inputFiles) {
        call md5 { 
            input: 
            inputFile=oneFile, 
            n_CPU = n_CPU,
            GB_memory = GB_memory,
            GB_disk = GB_disk
        }
    }

    call md5ofmd5s { 
        input: 
        inputFiles=md5.value, 
        n_CPU = n_CPU,
        GB_memory = GB_memory,
        GB_disk = GB_disk 
    }

}

task md5 {
    input {
        File inputFile
        Int n_CPU
        Float GB_memory
        Int GB_disk
    }

    command  <<<
        awk '{
            OFS="\t"; 
            if($7>10) { 
            print $1,$2-10,$2+10; 
            print $1,$3-10,$3+10;
            }
        }' ~{inputFile} | sortBed | uniq > "md5sum.tab"

    >>>

    output {
        File value = "md5sum.tab"
    }

    runtime {
        memory: "${GB_memory}GB"
        disks: "local-disk ${GB_disk} SSD"
        cpu: "${n_CPU}"
        docker: "bwhbioinformaticshub/bioinformatics-basic:latest"
    }
}

task md5ofmd5s {

    input {
        Array[File] inputFiles
        Int n_CPU
        Float GB_memory
        Int GB_disk
    }


    command <<<
        cat ~{sep=" " inputFiles} | sort | uniq -c | awk '{OFS="\t"; if($1>5) print $2,$3,$4}' > starjunctions.merged.splicingsites.flanking20nt.bed
    >>>

    output {
        File value = "starjunctions.merged.splicingsites.flanking20nt.bed"
    }

    runtime {
        memory: "${GB_memory}GB"
        disks: "local-disk ${GB_disk} SSD"
        cpu: "${n_CPU}"
        docker: "bwhbioinformaticshub/bioinformatics-basic:latest"
    }
}
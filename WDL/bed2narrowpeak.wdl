version 1.0

workflow bed2narrowpeak_workflow {

    input {
        File minus_bw_file
        File minus_bed_file
        String minus_file_name

        File plus_bw_file
        File plus_bed_file
        String plus_file_name

        Int n_CPU
        Float GB_memory
        Int GB_disk
    }

    call convert_to_narrowpeak as minusPeak {
        input: 
        n_CPU = n_CPU, 
        strand = "-", 
        GB_memory = GB_memory, 
        GB_disk = GB_disk, 
        input_bw = minus_bw_file, 
        input_bed = minus_bed_file, 
        name = minus_file_name
    }

    call convert_to_narrowpeak as plusPeak {
        input: 
        n_CPU = n_CPU, 
        strand = "+", 
        GB_memory = GB_memory, 
        GB_disk = GB_disk, 
        input_bw = plus_bw_file, 
        input_bed = plus_bed_file, 
        name = plus_file_name
    }

    call closest {
        input:
        plus_peaks = plusPeak.narrowPeak, 
        minus_peaks = minusPeak.narrowPeak, 
        n_CPU = n_CPU,
        GB_memory = GB_memory, 
        GB_disk = GB_disk
    }

    output {
        File plus_peaks = plusPeak.narrowPeak
        File minus_peaks = minusPeak.narrowPeak
        File plus_minus_bed = closest.plus_minus
        File minus_plus_bed = closest.minus_plus
    }

}

task convert_to_narrowpeak {

    input {
        String strand
        File input_bw
        File input_bed
        Int n_CPU
        Float GB_memory
        Int GB_disk
        String name 
    }

    command <<< 

    cat ~{input_bed} | while read chr start end name rest
    do
        l=$(expr $end - $start)
        bigWigSummary ~{input_bw} $chr $start $end $l | awk -vstrand=~{strand} -vchr=$chr -vstart=$start -vend=$end -vname=$name '{OFS="\t"; for(i=1;i<=NF;i++) if($i!="n/a" && $i>max) {imax=i;max=$i}}END{print chr, start+imax-1, start+imax, name, 0, strand, max, -1, -1, 0}'
    done > ~{name}.narrowPeak

    >>>

    output {
        File narrowPeak = "~{name}.narrowPeak"
    }

    runtime {
        docker: "rwang429/narrowpeak:latest"
        memory: "${GB_memory}GB"
        disks: "local-disk ${GB_disk} SSD"
        cpu: "${n_CPU}"
        preemptible: 1
        zones: "us-central1-a us-central1-b"

    }
}

task closest {

    input {
        File plus_peaks
        File minus_peaks
        Int n_CPU
        Float GB_memory
        Int GB_disk
    }

    command <<< 

    bedtools closest -a ~{minus_peaks} -b ~{plus_peaks} -t all -D b > minus.plus.bed 
    bedtools closest -a ~{plus_peaks} -b ~{minus_peaks} -t all -D a > plus.minus.bed
    >>>

    output {
        File plus_minus = "plus.minus.bed"
        File minus_plus = "minus.plus.bed"
    }

    runtime {
        docker: "rwang429/narrowpeak:latest"
        memory: "${GB_memory}GB"
        disks: "local-disk ${GB_disk} SSD"
        cpu: "${n_CPU}"
        preemptible: 1
        zones: "us-central1-a us-central1-b"
    }


}
## rna-bam2bigwig.wdl
##
## This workflow converted alignment file (.bam) to the normalized genomic coverage file (.bigwig).
##
## Inputs:
## - Per-sample:
##   - Alignment BAM file (.bam)
##   - Sample name
##   - Total mapped reads
##
## - Reference:
##   - genome chromosome size (a tab file with chr and size for each chromosome)
##
## Outputs:
## - Per-sample
##   - <sample-id>.plus.normalized.bw
##   - <sample-id>.minus.normalized.bw

workflow bam2bigwig_workflow {
    String sample_name
    File bam_file
    Int total_mapped_reads
    File chr_sizes
    Int n_CPU
    Float GB_memory
    Int GB_disk
    

    call bam2bigwig {
        input:
            bam_file = bam_file,
            sample_name = sample_name,
            total_mapped_reads = total_mapped_reads,
            chr_sizes = chr_sizes,
            n_CPU = n_CPU,
            GB_memory = GB_memory,
            GB_disk = GB_disk
    }
    output {
        File coverage_bigwig_plus = bam2bigwig.coverage_bigwig_plus
        File coverage_bigwig_minus = bam2bigwig.coverage_bigwig_minus
    }
}

task bam2bigwig {

    String sample_name
    File bam_file
    Int total_mapped_reads
    File chr_sizes
    Int n_CPU
    Float GB_memory
    Int GB_disk
	
    command {
        set -o errexit
        set -o nounset
        set -o pipefail
        set -o xtrace

        RPMscale=$(bc <<< "scale=6;1000000/${total_mapped_reads}") 

        echo "bam-->bw, by strand..."
        bedtools genomecov -ibam ${bam_file} -bg -scale $RPMscale -split -strand + | LC_COLLATE=C sort -k1,1 -k2,2n > ${sample_name}.plus.normalized.bedGraph && \
        bedGraphToBigWig ${sample_name}.plus.normalized.bedGraph ${chr_sizes} ${sample_name}.plus.normalized.bw && \
        rm ${sample_name}.plus.normalized.bedGraph

        bedtools genomecov -ibam ${bam_file} -bg -scale $RPMscale -split -strand - | LC_COLLATE=C sort -k1,1 -k2,2n > ${sample_name}.minus.normalized.bedGraph && \
        bedGraphToBigWig ${sample_name}.minus.normalized.bedGraph ${chr_sizes} ${sample_name}.minus.normalized.bw && \
        rm ${sample_name}.minus.normalized.bedGraph

    }

    output {
        File coverage_bigwig_plus = "${sample_name}.plus.normalized.bw"
        File coverage_bigwig_minus = "${sample_name}.minus.normalized.bw"
    }

    runtime {
        docker: "bwhbioinformaticshub/bioinformatics-basic:latest"
        memory: "${GB_memory}GB"
        disks: "local-disk ${GB_disk} SSD"
        cpu: "${n_CPU}"
        preemptible: 1
        zones: "us-central1-a us-central1-b"
    }

    meta {
        author: "Xianjun Dong, PhD"
        email: "xdong@rics.bwh.harvard.edu"
    }
}
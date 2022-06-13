## combine_bigwig.wdl
##
## WDL script to combine bigwig for multiple samples from the same group
##
## Inputs:
## - a file including a list of bigwig files
## - group label
## - genome chromosome size (a tab file with chr and size for each chromosome)
##
## Outputs:
## - combined bigwig: combined.mean.normalized.${group_label}.bigwig

task combine_biwig {

    String group_label
    File chr_sizes
    Int n_CPU
    Float GB_memory
    Int GB_disk
	Array[File] bigwig_files
	Int N = length(bigwig_files)
    
    command <<<
        set -o errexit
        set -o nounset
        set -o pipefail
        set -o xtrace

        echo "["`date`"] ## using ucsc-bigWigMerge to add signal values of multiple bigWigs together into a single output bedGraph"
        
        ## once Cromwell has localized the files, you'll then need to create the file list that you pass to bigWigMerge in the WDL task itself
        #echo "["`date`"] create file list"
        #echo ${sep=' ' bigwig_files} | tr " " "\n" > bigwig_file_list.txt
        
        echo "["`date`"] get the nuclear genome size"    
        GENOME_SIZE=$(grep -v "_" ${chr_sizes} | grep -v chrM | awk '{s+=$2}END{print s}') 
        
        echo "["`date`"] merge bigwig files togethr"    
		# bigWigMerge -inList bigwig_file_list.txt combined.sum.normalized.${group_label}.bedGraph.unsorted
        bigWigMerge ${sep=' ' bigwig_files} combined.sum.normalized.${group_label}.bedGraph.unsorted
        
        echo "["`date`"] sort"    
        LC_ALL=C sort -S 1G -k1,1 -k2,2n combined.sum.normalized.${group_label}.bedGraph.unsorted > combined.sum.normalized.${group_label}.bedGraph
          
        echo "["`date`"] average the signal" 
        awk -v N=${N} -v TOTAL=$GENOME_SIZE '{OFS="\t";$4=$4/N; S+=$4*($3-$2); if(id!=$4 || e!=$2 || chr!=$1) {if(chr!="") print chr,s,e,id; chr=$1;s=$2;e=$3;id=$4;} else {e=$3;}}END{print chr,s,e,id; bc=S/TOTAL;print "#basalCoverage="bc, "#"S"/"TOTAL;}' combined.sum.normalized.${group_label}.bedGraph \
        > combined.mean.normalized.${group_label}.bedGraph
        
        echo "["`date`"] bedgraph --> bigwig"
        bedGraphToBigWig combined.mean.normalized.${group_label}.bedGraph ${chr_sizes}  combined.mean.normalized.${group_label}.bigwig

    >>>

    output {
        File combined_biwig = "combined.mean.normalized.${group_label}.bigwig"
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


workflow combine_biwig_workflow {
    String group_label
    File bigwig_files_list
    File chr_sizes
    Int n_CPU
    Float GB_memory
    Int GB_disk
	Array[File] bigwig_files = read_lines(bigwig_files_list)

    call combine_biwig {
        input:
            bigwig_files = bigwig_files,
            chr_sizes = chr_sizes,
            group_label = group_label,
            n_CPU = n_CPU,
            GB_memory = GB_memory,
            GB_disk = GB_disk
    }
    output {
        File combined_bigwig = combine_biwig.combined_biwig
    }
}
module load gsutil/default
module load ucsc/default

bigWigAverageOverBed combined.mean.normalized.random.samplesN200.minus.bigwig /data/bioinformatics/referenceGenome/Homo_sapiens/UCSC/hg38/Annotation/Genes/gencode.v37.annotation.gtf.genes.bed minus.tab  
bigWigAverageOverBed combined.mean.normalized.random.samplesN200.plus.bigwig /data/bioinformatics/referenceGenome/Homo_sapiens/UCSC/hg38/Annotation/Genes/gencode.v37.annotation.gtf.genes.bed plus.tab
 
bigWigAverageOverBed combined.mean.normalized.random.samplesN200.plus.bigwig /data/bioinformatics/referenceGenome/Homo_sapiens/UCSC/hg38/Annotation/Genes/gencode.v37.annotation.bed12 transcripts_plus.tab
bigWigAverageOverBed combined.mean.normalized.random.samplesN200.minus.bigwig /data/bioinformatics/referenceGenome/Homo_sapiens/UCSC/hg38/Annotation/Genes/gencode.v37.annotation.bed12 transcripts_minus.tab

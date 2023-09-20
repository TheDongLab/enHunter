#!/bin/bash

conda activate /PHShome/rw552/condaenvs/ucsc 
inputs=/data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs/eRNA_exp/chr16
GENOME=/data/bioinformatics/referenceGenome/Homo_sapiens/UCSC/hg38

bedtools getfasta -name -fi $GENOME/Sequence/WholeGenomeFasta/genome.fa -bed $inputs/eRNA.LITAF.loci.bed -fo eRNA.LITAF.loci.fasta

## run fimo on website 
fimo --oc . --verbosity 1 --bgfile --nrdb-- --thresh 1.0E-4 complete-factorbook-catalog.meme eRNA.LITAF.loci.fasta

## find the genomic locations of snps in LITAF locus : eRNA.LITAF.snps.bed 
bedtools intersect -a eRNA.LITAF.snps.sorted.bed -b fimo.motifs.bed -wao

join fimo.id.sorted.tsv eRNA.LIAF.snps.and.motifs.sorted.txt

join eRNA.LIAF.snps.and.motifs.uniq.txt fimo.id.sorted.tsv | sed "s/ /\t/g" > fimo.and.LITAF.snps.txt 
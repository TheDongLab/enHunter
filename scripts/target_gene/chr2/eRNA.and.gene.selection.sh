#!/bin/bash

# generating gene and eRNA TPM/RPM table

# CXCR2
cd /data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs/gene_exp
cat <(head -n 1 gene_expr_matrix_tpm_row_genes.txt) <(grep ENSG00000180871 gene_expr_matrix_tpm_row_genes.txt) > CXCR2_gene_tpm_table.txt

### file: /data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs/eRNA.merged.readCounts.v2.xls
RPM=/data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs/eRNA.merged.readCounts.v2.xls

### stored in the /data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs/eRNA_exp directory 

# chr2_218155750_218156200 (minus)
cat <(head -n 1 $RPM) <(grep chr2_218155750_218156200_minus $RPM) > chr2_218155750_218156200_minus_RPM_table.txt

# chr2_218156440_218156680 (plus)
cat <(head -n 1 $RPM) <(grep chr2_218156440_218156680_plus $RPM) > chr2_218156440_218156680_plus_RPM_table.txt

### find neighboring gene TPM expression 
# CXCR2_loci.bed
# chr2  218155750 218156680 chr2_218155750_218156680_eRNA
./neighboring_genes.sh CXCR2_loci.bed 1000000 /data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs/gene_exp/gene_expr_matrix_tpm_row_genes.txt CXCR2_locus_expr




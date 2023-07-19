#!/bin/bash

# generating gene and eRNA TPM/RPM table

# CXCR6 
cd /data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs/gene_exp
cat <(head -n 1 gene_expr_matrix_tpm_row_genes.txt) <(grep ENSG00000189067 gene_expr_matrix_tpm_row_genes.txt) > CXCR6_gene_exp_table.txt

### file: /data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs/eRNA.merged.readCounts.v2.xls
RPM=/data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs/eRNA.merged.readCounts.v2.xls

### stored in the /data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs/eRNA_exp directory 

# chr3_46089640_46089850 (minus)
cat <(head -n 1 $RPM) <(grep chr3_46089640_46089850_minus $RPM) > chr3_46089640_46089850_minus_exp_table.txt

# chr3_46090120_46090230 (plus)
cat <(head -n 1 $RPM) <(grep chr3_46090120_46090230_plus $RPM) > chr3_46090120_46090230_plus_exp_table.txt
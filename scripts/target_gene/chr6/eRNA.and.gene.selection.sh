#!/bin/bash

# generating gene and eRNA TPM/RPM table

# HLA-E
### stored in the /data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs/eRNA_exp directory 
RPM_plus=/data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs/plus/eRNA.meanRPM.xls
# chr6_30497480_30497870 (plus)
cat <(head -n 1 $RPM_plus) <(grep chr6_30497480_30497870 $RPM_plus) > chr6_30497480_30497870_plus_RPM_table.txt

### find neighboring gene TPM expression 
# HLA_E_loci.bed
# chr6  30497480 30497870 chr6_30497480_30497870_plus_eRNA
./neighboring_genes.sh HLA_E_loci.bed 1000000 /data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs/gene_exp/gene_expr_matrix_tpm_row_genes.txt HLA_E_locus_expr




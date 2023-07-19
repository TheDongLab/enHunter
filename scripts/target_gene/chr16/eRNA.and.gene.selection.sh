#!/bin/bash

# generating gene and eRNA TPM/RPM table

# LITAF 
cd /data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs
cat <(head -n 1 gene_expr_matrix_tpm_row_genes.txt) <(grep ENSG00000189067 gene_expr_matrix_tpm_row_genes.txt) > LITAF_gene_exp_table.txt

# CLEC16A
cd /data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs
cat <(head -n 1 gene_expr_matrix_tpm_row_genes.txt) <(grep ENSG00000038532 gene_expr_matrix_tpm_row_genes.txt) > CLEC16A_gene_exp_table.txt


# stored in the /data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs/eRNA_exp directory

# chr16_11611980_11612400 (eRNA 1 minus)
cd /data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs/minus
cat <(head -n 1 eRNA.meanRPM.xls) <(grep chr16_11611980_11612400 eRNA.meanRPM.xls) > chr16_11611980_11612400_exp_table.txt

# chr16_11612780_11613560 (eRNA 1 plus)
cd /data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs/plus
cat <(head -n 1 eRNA.meanRPM.xls) <(grep chr16_11612780_11613560 eRNA.meanRPM.xls) > chr16_11612780_11613560_exp_table.txt

# chr16_11613470_11613780 (eRNA 2 minus)
cd /data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs/minus
cat <(head -n 1 eRNA.meanRPM.xls) <(grep chr16_11613470_11613780 eRNA.meanRPM.xls) > chr16_11613470_11613780_exp_table.txt

# chr16_11613950_11614560 (eRNA 2 plus)
cd /data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs/plus
cat <(head -n 1 eRNA.meanRPM.xls) <(grep chr16_11613950_11614560 eRNA.meanRPM.xls) > chr16_11613950_11614560_exp_table.txt

# chr16_11639850_11640300 (eRNA 3 minus)
cd /data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs/minus
cat <(head -n 1 eRNA.meanRPM.xls) <(grep chr16_11639850_11640300 eRNA.meanRPM.xls) > chr16_11639850_11640300_exp_table.txt

# chr16_11640470_11641120 (eRNA 3 plus)
cd /data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs/plus
cat <(head -n 1 eRNA.meanRPM.xls) <(grep chr16_11640470_11641120 eRNA.meanRPM.xls) > chr16_11640470_11641120_exp_table.txt


# keeps track of the genes tested with the eRNAs  

scripts="/Users/rw552/Documents/amp-pd/enHunter/scripts/pair-exp"

## generate paired RPM information 
cd /Users/rw552/Documents/amp-pd/enHunter/input_files/eRNAs
Rscript $scripts/make_pairs.R chr16_11611980_11612400.txt chr16_11612780_11613560.txt eRNA1.txt

Rscript $scripts/make_pairs.R chr16_11613470_11613780.txt chr16_11613950_11614560.txt eRNA2.txt

Rscript $scripts/make_pairs.R chr16_11639850_11640300.txt chr16_11640470_11641120.txt eRNA3.txt

# run the analysis for each gene 
input_dir="/Users/rw552/Documents/amp-pd/enHunter/input_files"

Rscript $scripts/target_gene_exp_cor.R $input_dir/genes/SNN_gene_exp.txt ENSG00000184602.5 $input_dir/eRNAs/eRNA1.txt chr16_11611980_11612400_chr16_11612780_11613560 $input_dir/eRNAs/eRNA2.txt chr16_11613470_11613780_chr16_11613950_11614560 $input_dir/eRNAs/eRNA3.txt chr16_11639850_11640300_chr16_11640470_11641120
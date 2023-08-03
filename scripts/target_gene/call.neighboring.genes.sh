#!/bin/bash

dir=./input_files/eRNAs/tables

# see neighboringgenes.sh to generate the *locus_expr.txt files


### LITAF loci
Rscript ./scripts/target_gene/neighboring.genes.cor.R $dir/neighboring_LITAF_locus_expr.txt "eRNA1_LITAF" $dir/chr16_11611980_11612400_exp_table.txt $dir/chr16_11612780_11613560_exp_table.txt
Rscript ./scripts/target_gene/neighboring.genes.cor.R $dir/neighboring_LITAF_locus_expr.txt "eRNA2_LITAF" $dir/chr16_11613470_11613780_exp_table.txt $dir/chr16_11613950_11614560_exp_table.txt
Rscript ./scripts/target_gene/neighboring.genes.cor.R $dir/neighboring_LITAF_locus_expr.txt "eRNA3_LITAF" $dir/chr16_11639850_11640300_exp_table.txt $dir/chr16_11640470_11641120_exp_table.txt

### HLA-E loci 
Rscript ./scripts/target_gene/neighboring.genes.cor.R $dir/HLA_E_locus_expr.txt "eRNA_HLA_E" $dir/chr6_30497480_30497870_plus_RPM_table.txt
# chr6:29,497,480-31,497,870
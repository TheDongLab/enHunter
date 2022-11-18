#!/bin/bash
# keeps track of all the graphs and analysis to generate 
wd=/Users/rosanwang/Documents/college/dong_lab/code/enHunter

Rscript eRNA.GWAS.enrichment.R $wd/input_files/characterization/feature.enrichment/counts/minus/eRNA.minus.f16.GWASDisease.counts $wd/input_files/characterization/feature.enrichment/counts/GWAS_20220810.v1.02.counts.v2 ./minus/GWAS_minus_signif.png
Rscript eRNA.GWAS.enrichment.R $wd/input_files/characterization/feature.enrichment/counts/plus/eRNA.plus.f16.GWASDisease.counts $wd/input_files/characterization/feature.enrichment/counts/GWAS_20220810.v1.02.counts.v2 ./plus/GWAS_plus_signif.png

Rscript fisher_test.R $wd/input_files/characterization/feature.enrichment/counts/minus/eRNA.minus.f08.CAGEenhtissue.counts $wd/input_files/characterization/feature.enrichment/counts/all_enhancer_tissues_hg38_CAGE.counts ./minus/CAGE_minus.png
Rscript fisher_test.R $wd/input_files/characterization/feature.enrichment/counts/plus/eRNA.plus.f08.CAGEenhtissue.counts $wd/input_files/characterization/feature.enrichment/counts/all_enhancer_tissues_hg38_CAGE.counts ./plus/CAGE_plus.png 

Rscript fisher_test.R $wd/input_files/characterization/feature.enrichment/counts/minus/eRNA.minus.f18.eSNP.gtexCaviarDisease.counts $wd/input_files/characterization/feature.enrichment/counts/snp_gtexCaviar.counts ./minus/gtexCaviar_minus.png
Rscript fisher_test.R $wd/input_files/characterization/feature.enrichment/counts/plus/eRNA.plus.f18.eSNP.gtexCaviarDisease.counts $wd/input_files/characterization/feature.enrichment/counts/snp_gtexCaviar.counts ./plus/gtexCaviar_plus.png

Rscript fisher_test.R $wd/input_files/characterization/feature.enrichment/counts/minus/eRNA.minus.f18.eSNP.gtexDapgDisease.counts $wd/input_files/characterization/feature.enrichment/counts/snp_gtexDapg.counts ./minus/gtexDapg_minus.png
Rscript fisher_test.R $wd/input_files/characterization/feature.enrichment/counts/plus/eRNA.plus.f18.eSNP.gtexDapgDisease.counts $wd/input_files/characterization/feature.enrichment/counts/snp_gtexDapg.counts ./plus/gtexDapg_plus.png

Rscript fisher_test.R $wd/input_files/characterization/feature.enrichment/counts/minus/eRNA.minus.f18.eSNP.Disease.counts $wd/input_files/characterization/feature.enrichment/counts/all_signif_varient_gene_pairs_pval_eQTL.counts ./minus/pval_minus.png
Rscript fisher_test.R $wd/input_files/characterization/feature.enrichment/counts/plus/eRNA.plus.f18.eSNP.Disease.counts $wd/input_files/characterization/feature.enrichment/counts/all_signif_varient_gene_pairs_pval_eQTL.counts ./plus/pval_plus.png

#### ONLY CLASS 1 PAIRS 
wd=/Users/rosanwang/Documents/college/dong_lab/code/enHunter
Rscript eRNA.GWAS.enrichment.R $wd/input_files/characterization/feature.enrichment/counts/minus/class1/eRNA.minus.class1.f16.GWASDisease.counts.v2 $wd/input_files/characterization/feature.enrichment/counts/GWAS_20220810.v1.02.counts.v2 GWAS_class1_minus_signif.png

Rscript fisher_test.R  $wd/input_files/characterization/feature.enrichment/counts/minus/class1/eRNA.minus.class1.f08.CAGEenhtissue.counts $wd/input_files/characterization/feature.enrichment/counts/all_enhancer_tissues_hg38_CAGE.counts ./minus/class1/CAGE_class1_minus.png

Rscript fisher_test.R  $wd/input_files/characterization/feature.enrichment/counts/minus/class1/eRNA.minus.class1.f18.eSNP.gtexCaviarDisease.counts $wd/input_files/characterization/feature.enrichment/counts/snp_gtexCaviar.counts ./minus/class1/Caviar_class1_minus.png
Rscript fisher_test.R  $wd/input_files/characterization/feature.enrichment/counts/minus/class1/eRNA.minus.class1.f18.eSNP.gtexDapgDisease.counts $wd/input_files/characterization/feature.enrichment/counts/snp_gtexDapg.counts ./minus/class1/Dapg_class1_minus.png
Rscript fisher_test.R  $wd/input_files/characterization/feature.enrichment/counts/minus/class1/eRNA.minus.class1.f18.eSNP.pvalDisease.counts $wd/input_files/characterization/feature.enrichment/counts/all_signif_varient_gene_pairs_pval_eQTL.counts ./minus/class1/pval_class1_minus.png


Rscript eRNA.GWAS.enrichment.R $wd/input_files/characterization/feature.enrichment/counts/plus/class1/eRNA.plus.class1.f16.GWASDisease.counts.v2 $wd/input_files/characterization/feature.enrichment/counts/GWAS_20220810.v1.02.counts.v2 ./plus/class1/GWAS_class1_plus_signif.png

Rscript fisher_test.R  $wd/input_files/characterization/feature.enrichment/counts/plus/class1/eRNA.plus.class1.f18.eSNP.gtexCaviarDisease.counts $wd/input_files/characterization/feature.enrichment/counts/snp_gtexCaviar.counts ./plus/class1/Caviar_class1_plus.png
Rscript fisher_test.R  $wd/input_files/characterization/feature.enrichment/counts/plus/class1/eRNA.plus.class1.f18.eSNP.gtexDapgDisease.counts $wd/input_files/characterization/feature.enrichment/counts/snp_gtexDapg.counts ./plus/class1/Dapg_class1_plus.png
Rscript fisher_test.R  $wd/input_files/characterization/feature.enrichment/counts/plus/class1/eRNA.plus.class1.f08.CAGEenhtissue.counts $wd/input_files/characterization/feature.enrichment/counts/all_enhancer_tissues_hg38_CAGE.counts ./plus/class1/CAGE_class1_plus.png
Rscript fisher_test.R  $wd/input_files/characterization/feature.enrichment/counts/plus/class1/eRNA.plus.class1.f18.eSNP.pvalDisease.counts $wd/input_files/characterization/feature.enrichment/counts/all_signif_varient_gene_pairs_pval_eQTL.counts ./plus/class1/pval_class1_plus.png


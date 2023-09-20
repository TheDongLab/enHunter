### finding the significant PD GWAS snps from overlapping blood and dopamine TNEs 


### SNPs: /data/bioinformatics/external_data/externalData/signif.Nallsetal_META_ANALYSIS_10K23_beta_se_wo23AndMe1_pdgene_sharing.tbl 
### convert tbl to bed file 

conda activate /PHShome/rw552/condaenvs/ucsc
AMPPD_eRNA=/data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs
externalData=/data/bioinformatics/external_data/externalData

# plus strand 
intersectBed -a $AMPPD_eRNA/plus/eRNA.bed -b $externalData/signif.Nallsetal_META_ANALYSIS_10K23_beta_se_wo23AndMe1_pdgene_sharing.tbl -wao | less

# minus strand
intersectBed 
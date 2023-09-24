# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# 
# Be sure to use an up to date version of R and Matrix eQTL.

# source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
library(MatrixEQTL)

## Location of the package with the data files.
base.dir = "/data/bioinformatics/projects/donglab/AMPPD_eRNA/inputs/eQTL"

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = paste0(base.dir, "/snps/eRNA.SNP_Matrix.revised.txt");
snps_location_file_name = paste0(base.dir, "/snps/snps.coordinates");

# Gene expression file name
expression_file_name = paste0(base.dir, "/genes/eRNA.gene.Matrix.reordered.txt");
gene_location_file_name = paste0(base.dir, "/genes/gene.eRNA.coordinates");

# Covariates file name
# Set to character() for no covariates
covariates_file_name = paste0(base.dir, "/covariate.reordered.txt");


# Output file name
output_file_name_cis = "cis.eQTL";
output_file_name_tra = "trans.eQTL";

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 2e-2;
pvOutputThreshold_tra = 1e-2;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 1e6;

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

# for some reason it keeps reading snpid as a column 
snpcols <- snps$columnNames
colnames(snps) <- snpcols[ !snpcols == "snpid"]

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}


## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt,
  output_file_name      = "", #running only cis-eQTLs
  pvOutputThreshold     = 0, #running only cis-eQTLs
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE, 
  output_file_name.cis  = output_file_name_cis,
  pvOutputThreshold.cis = 1, # record all pairs: ntests = neqtls
  snpspos = snpspos, 
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name_cis);
unlink(output_file_name_tra);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
write.table(me$cis$eqtls, file = "cis.eQTLs.xls")


## Make the histogram of local and distant p-values
pdf(paste("cis.eQTL", "pdf", sep="."))
plot(me)
dev.off()







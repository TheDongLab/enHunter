### making a sample manhattan plot ###

# Load the library
library(qqman)

# Make the Manhattan plot on the gwasResults dataset
manhattan(gwasResults, chr="CHR", bp="BP", snp="SNP", p="P" )

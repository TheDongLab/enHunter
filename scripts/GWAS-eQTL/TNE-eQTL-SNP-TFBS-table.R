# creating merge table for ALL eQTL and GWAS snp information 

# usage: Rscript TNE-eQTL-SNP-TFBS-table.R eRNA.plus.f18.eSNP.gtexDapgDisease.intersect eRNA.minus.f18.eSNP.gtexDapgDisease.intersect eRNA.plus.f18.eSNP.gtexCaviarDisease.intersect eRNA.minus.f18.eSNP.gtexCaviarDisease.intersect eRNA.mplus.f16.GWASDisease.intersect eRNA.minus.f16.GWASDisease.intersect 

args<-commandArgs(trailingOnly=TRUE)

#plus strand Dapg SNPs
system(paste0("awk 'OFS=\"", "\t\"", " {print $4, \"+\", $8, $5\"_\"$6\"_\"$7, $9, \"Yes\"}' ", args[1], "> eRNA.f18.plus.eSNP.gtexDapg.txt"))

# minus strand Dapg SNPs
system(paste0("awk 'OFS=\"\t\" {print $4, \"-\", $8, $5\"_\"$6\"_\"$7, $9, \"Yes\"}' ", args[2], "> eRNA.f18.minus.eSNP.gtexDapg.txt"))

system(paste0("cat eRNA.f18.plus.eSNP.gtexDapg.txt eRNA.f18.minus.eSNP.gtexDapg.txt > eRNA.f18.eSNP.gtexDapg.txt"))

#plus strand Caviar SNPs
system(paste0("awk 'OFS=\"\t\" {print $4, \"+\", $8, $5\"_\"$6\"_\"$7, $9, \"Yes\"}' ", args[3], "> eRNA.f18.plus.eSNP.gtexCaviar.txt") )

# minus strand Caviar SNPs
system(paste0("awk 'OFS=\"\t\" {print $4, \"-\", $8, $5\"_\"$6\"_\"$7, $9, \"Yes\"}' ", args[4], "> eRNA.f18.minus.eSNP.gtexCaviar.txt") )

system(paste0("cat eRNA.f18.plus.eSNP.gtexCaviar.txt eRNA.f18.minus.eSNP.gtexCaviar.txt > eRNA.f18.eSNP.gtexCaviar.txt") )

Dapg <- read.table("eRNA.f18.eSNP.gtexDapg.txt")
colnames(Dapg) <- c("TNE", "strand", "rsid", "eQTL_pos", "tissue", "dapg")

Caviar <- read.table("eRNA.f18.eSNP.gtexCaviar.txt")
colnames(Caviar) <- c("TNE", "strand", "rsid", "eQTL_pos", "tissue", "caviar")

all_eQTL <- merge(Dapg, Caviar, by=c("TNE", "strand", "rsid" ,"eQTL_pos", "tissue"), all=TRUE)

write.table(all_eQTL, "eRNA.eQTL.snps", sep="\t", quote = F, col.names = T, row.names = F)

## combining the GWAS information 
system(paste0("awk 'OFS=\"\t\"; FS=\"\t\" {print $4, \"+\", $8, $5\"_\"$6\"_\"$7, $9}' ", args[5], "> GWAS_plus.txt"))
system(paste0("awk 'OFS=\"\t\"; FS=\"\t\" {print $4, \"-\", $8, $5\"_\"$6\"_\"$7, $9}' ", args[6], "> GWAS_minus.txt"))

system("cat GWAS_minus.txt GWAS_plus.txt > GWAS.txt")

GWAS <- read.table("GWAS.txt")
colnames(GWAS) <- c("TNE", "strand", "rsid", "eQTL_pos", "disease")

eQTL_GWAS <- merge(all_eQTL, GWAS, sep="\t", all=T)

write.table(eQTL_GWAS, "eRNA.eQTL.snps.GWAS", sep="\t", quote = F, col.names = T, row.names = F)



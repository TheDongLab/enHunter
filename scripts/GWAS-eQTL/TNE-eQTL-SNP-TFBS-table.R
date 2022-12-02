# creating merge table for ALL eQTL and GWAS snp information 

# usage: Rscript TNE-eQTL-SNP-TFBS-table.R caviar.table dapg.table
# eRNA.plus.f18.eSNP.gtexDapgDisease.intersect eRNA.minus.f18.eSNP.gtexDapgDisease.intersect eRNA.plus.f18.eSNP.gtexCaviarDisease.intersect eRNA.minus.f18.eSNP.gtexCaviarDisease.intersect eRNA.plus.f16.GWASDisease.intersect eRNA.minus.f16.GWASDisease.intersect 

args<-commandArgs(trailingOnly=TRUE)

#plus strand Dapg SNPs
system(paste0("awk 'OFS=\"", "\t\"", " {print $4, \"+\", $8, $5\"_\"$6\"_\"$7, $9, \"Yes\"}' eRNA.plus.f18.eSNP.gtexDapgDisease.intersect > eRNA.f18.plus.eSNP.gtexDapg.txt"))

# minus strand Dapg SNPs
system(paste0("awk 'OFS=\"\t\" {print $4, \"-\", $8, $5\"_\"$6\"_\"$7, $9, \"Yes\"}' eRNA.minus.f18.eSNP.gtexDapgDisease.intersect > eRNA.f18.minus.eSNP.gtexDapg.txt"))

system(paste0("cat eRNA.f18.plus.eSNP.gtexDapg.txt eRNA.f18.minus.eSNP.gtexDapg.txt > eRNA.f18.eSNP.gtexDapg.txt"))

#plus strand Caviar SNPs
system(paste0("awk 'OFS=\"\t\" {print $4, \"+\", $8, $5\"_\"$6\"_\"$7, $9, \"Yes\"}' eRNA.plus.f18.eSNP.gtexCaviarDisease.intersect > eRNA.f18.plus.eSNP.gtexCaviar.txt") )

# minus strand Caviar SNPs
system(paste0("awk 'OFS=\"\t\" {print $4, \"-\", $8, $5\"_\"$6\"_\"$7, $9, \"Yes\"}' eRNA.minus.f18.eSNP.gtexCaviarDisease.intersect > eRNA.f18.minus.eSNP.gtexCaviar.txt") )

system(paste0("cat eRNA.f18.plus.eSNP.gtexCaviar.txt eRNA.f18.minus.eSNP.gtexCaviar.txt > eRNA.f18.eSNP.gtexCaviar.txt") )

# read in the Caviar and Dapg table 
# V19: tissue 
# V1: chromosome 
# V7: eQTL start 
# V8: eQTL end 
# V14: rsid 
# V16: gene name 
caviar <- read.table(args[1])
caviar <- caviar[, c(1,7,8,14,16,19)]
colnames(caviar) <- c("chr", "start", "end", "rsid", "gene", "tissue")
caviar$eQTL_pos <- paste(caviar$chr, caviar$start, caviar$end, sep="_")
#caviar[,-c("chr", "start", "end")]
caviar[,-c(1:3)]

dapg <- read.table(args[2])
dapg <- dapg[, c(1,7,8,14,16,19)]
colnames(dapg) <- c("chr", "start", "end", "rsid", "gene", "tissue")
dapg$eQTL_pos <- paste(dapg$chr, dapg$start, dapg$end, sep="_")
#dapg <- dapg[, -c("chr", "start", "end")]
dapg <- dapg[, -c(1:3)]

TNE_Dapg <- read.table("eRNA.f18.eSNP.gtexDapg.txt")
colnames(TNE_Dapg) <- c("TNE", "strand", "rsid", "eQTL_pos", "tissue", "dapg")

TNE_Dapg <- merge(TNE_Dapg, dapg, by=c("eQTL_pos", "tissue", "rsid"), all.x = TRUE)

TNE_Caviar <- read.table("eRNA.f18.eSNP.gtexCaviar.txt")
colnames(TNE_Caviar) <- c("TNE", "strand", "rsid", "eQTL_pos", "tissue", "caviar")

TNE_Dapg <- merge(TNE_Caviar , caviar, by=c("eQTL_pos", "tissue", "rsid"), all.x = TRUE)

all_eQTL <- merge(TNE_Dapg, TNE_Caviar, by=c("TNE", "strand", "rsid" ,"eQTL_pos", "tissue"), all=TRUE)

write.table(all_eQTL, "eRNA.eQTL.snps", sep="\t", quote = F, col.names = T, row.names = F)

## combining the GWAS information 
system(paste0("awk 'OFS=\"\t\"; FS=\"\t\" {print $4, \"+\", $8, $5\"_\"$6\"_\"$7, $9}' eRNA.plus.f16.GWASDisease.intersect > GWAS_plus.txt"))
system(paste0("awk 'OFS=\"\t\"; FS=\"\t\" {print $4, \"-\", $8, $5\"_\"$6\"_\"$7, $9}'  eRNA.minus.f16.GWASDisease.intersect > GWAS_minus.txt"))

system("cat GWAS_minus.txt GWAS_plus.txt > GWAS.txt")

GWAS <- read.table("GWAS.txt")
colnames(GWAS) <- c("TNE", "strand", "rsid", "eQTL_pos", "disease")

eQTL_GWAS <- merge(all_eQTL, GWAS, sep="\t", all=T)

write.table(eQTL_GWAS, "eRNA.eQTL.snps.GWAS", sep="\t", quote = F, col.names = T, row.names = F)



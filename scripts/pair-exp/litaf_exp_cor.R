library(ggplot2)
library(dplyr)
library(GGally)

make_vector <- function(file, name) {
  vector <- scan(file = file, sep = "\t", nlines = 1, na.strings = name) %>% log(base=10) 
  vector[is.infinite(vector)] <- NA
  #return(vector[!is.na(vector) & is.finite(vector)])
  return(vector)
}

LITAF <- make_vector("./input_files/LITAF/LITAF_gene_exp.tsv", "ENSG00000189067.12")[1:8258]

# eRNA 1 minus 
chr16_11611980_11612400 <- make_vector("./input_files/LITAF/chr16_11611980_11612400.txt", 
                                       "chr16_11611980_11612400") 
# eRNA 1 plus 
chr16_11612780_11613560 <- make_vector("./input_files/LITAF/chr16_11612780_11613560.txt", 
                                       "chr16_11612780_11613560") 
# eRNA2 plus 
chr16_11613950_11614560 <- make_vector("./input_files/LITAF/chr16_11613950_11614560.txt", 
                                       "chr16_11613950_11614560") 
# eRNA2 minus 
chr16_11613470_11613780 <- make_vector("./input_files/LITAF/chr16_11613470_11613780.txt", 
                                       "chr16_11613470_11613780") 
# eRNA3 plus 
chr16_11640470_11641120 <- make_vector("./input_files/LITAF/chr16_11640470_11641120.txt", 
                                       "chr16_11640470_11641120") 
# eRNA3 minus 
chr16_11639850_11640300 <- make_vector("./input_files/LITAF/chr16_11639850_11640300.txt", 
                                       "chr16_11639850_11640300") 


max_ln <- max(c(length(chr16_11611980_11612400), length(chr16_11612780_11613560), 
                length(chr16_11613950_11614560), length(chr16_11640470_11641120), 
                length(chr16_11639850_11640300), length(chr16_11613470_11613780)))

uni_len <- function(vector, ln) {
  return(c(vector,rep(NA, ln - length(vector))))
}

# each eRNA run seperately 
df <- data.frame(LITAF = uni_len(LITAF, max_ln),
                 eRNA_1_minus = uni_len(chr16_11611980_11612400, max_ln), 
                 eRNA_1_plus = uni_len(chr16_11612780_11613560, max_ln), 
                 eRNA_2_minus = uni_len(chr16_11613470_11613780, max_ln),
                 eRNA_2_plus = uni_len(chr16_11613950_11614560, max_ln), 
                 eRNA_3_plus = uni_len(chr16_11640470_11641120, max_ln), 
                 eRNA_3_minus = uni_len(chr16_11639850_11640300, max_ln))


eRNAs <- ggpairs(df)
eRNAs

# running eRNAs as pairs 


########################## TESTING ##########################################
# why does this keep giving me NaN
correlation <- cor(chr16_11612780_11613560, chr16_11639850_11640300, method = "pearson")

# contains infinite values 
# contians -Inf
any(!is.finite(chr16_11612780_11613560))
test <- scan(file = "./input_files/LITAF/chr16_11612780_11613560.txt", sep = "\t",
             nlines = 1, na.strings = "chr16_11612780_11613560")
test <- test[!is.na(test) & test != 0] 

index <- which(!is.finite(chr16_11612780_11613560)) 

test[1686] 
# input value of 0 results in -Inf
########################TEST END############################################

read_vector <- function(file, name) {
  vector <- scan(file = file, sep = "\t", nlines = 1, na.strings = name)[-1] 
  return(vector)
}

# eRNA 1
eRNA_1 <- (read_vector("./input_files/LITAF/chr16_11611980_11612400.txt", "chr16_11611980_11612400") + 
  read_vector("./input_files/LITAF/chr16_11612780_11613560.txt", "chr16_11612780_11613560")) %>% log(base=10) 

# eRNA2 
eRNA_2 <- (read_vector("./input_files/LITAF/chr16_11613950_11614560.txt", "chr16_11613950_11614560") + 
             read_vector("./input_files/LITAF/chr16_11613470_11613780.txt", "chr16_11613470_11613780")) %>% log(base=10) 

# eRNA3  
eRNA_3 <- (read_vector("./input_files/LITAF/chr16_11639850_11640300.txt", "chr16_11639850_11640300") + 
             read_vector("./input_files/LITAF/chr16_11640470_11641120.txt", "chr16_11640470_11641120")) %>% log(base=10) 

df_eRNA <- data.frame(LITAF, eRNA_1, eRNA_2, eRNA_3)

eRNA_pairs <- ggpairs(df_eRNA)
eRNA_pairs

## testing transcripts instead of LITAF gene 
ENST00000571688.5 <- read.table("./input_files/LITAF/notebooks_ENST00000571688.5.tsv")[["TPM"]] %>% log(base=10) 

df_transcript <- data.frame(ENST00000571688.5[1:8258], eRNA_1, eRNA_2, eRNA_3)
eRNA_transcript <- ggpairs(df_transcript)
eRNA_transcript

cor(df_transcript, method = "pearson")
##### ALL Transcripts #####
transcripts <- read.table("./input_files/LITAF/LITAF_transcripts.tsv")
split_by_transcript <- split(transcripts, f = transcripts$Name)  

for (df in split_by_transcript) {
  
  TPM <- df[["TPM"]] %>% log(base=10) 
  TPM <- TPM[!is.na(TPM) & is.finite(TPM)]
  transcript <- df[["Name"]][[1]]
  
  
  max_ln <- max(c(length(eRNA_1), length(eRNA_2), 
                  length(eRNA_3), length(TPM)))
  
  full_df <- data.frame(uni_len(TPM, max_ln), uni_len(eRNA_1, max_ln), 
                        uni_len(eRNA_2, max_ln), uni_len(eRNA_3, max_ln) )
  
  eRNA_transcript <- ggpairs(full_df)
  ggsave(paste0(transcript, ".pdf"), plot = eRNA_transcript, width = 10, height = 10)
  
}

##### CLEC16A Gene ########
ENSG00000038532.15 <- read.table("./input_files/CLEC16A/clec16a_genes.tsv")[["TPM"]] %>% log(base=10) 
ENSG00000038532.15 <- ENSG00000038532.15[!is.na(ENSG00000038532.15) & is.finite(ENSG00000038532.15)]

max_ln_clec16a <- max(c(length(eRNA_1), length(eRNA_2), length(eRNA_3), length(ENSG00000038532.15)))

clec16a <- data.frame(clec16a = uni_len(ENSG00000038532.15, max_ln_clec16a), eRNA_1 = uni_len(eRNA_1, max_ln_clec16a), 
           eRNA_2 = uni_len(eRNA_2, max_ln_clec16a), eRNA_3 = uni_len(eRNA_3, max_ln_clec16a) )

ggpairs(clec16a)

###### testing a new method ... #######

read_vector <- function(file, name) {
  vector <- scan(file = file, sep = "\t", nlines = 1, na.strings = name)
}

# eRNA 1
eRNA_1 <- (read_vector("./input_files/LITAF/chr16_11611980_11612400.txt", "chr16_11611980_11612400") + 
             read_vector("./input_files/LITAF/chr16_11612780_11613560.txt", "chr16_11612780_11613560")) %>% log(base=10) 

# eRNA2 
eRNA_2 <- (read_vector("./input_files/LITAF/chr16_11613950_11614560.txt", "chr16_11613950_11614560") + 
             read_vector("./input_files/LITAF/chr16_11613470_11613780.txt", "chr16_11613470_11613780")) %>% log(base=10) 

# eRNA3  
eRNA_3 <- (read_vector("./input_files/LITAF/chr16_11639850_11640300.txt", "chr16_11639850_11640300") + 
             read_vector("./input_files/LITAF/chr16_11640470_11641120.txt", "chr16_11640470_11641120")) %>% log(base=10)

LITAF <- read_vector("./input_files/LITAF/LITAF_gene_exp.tsv", "ENSG00000189067.12")[1:length(eRNA_1)]

df <- data.frame(eRNA1 = eRNA_1, 
                 eRNA2 = eRNA_2, 
                 eRNA3 = eRNA_3,
                 LITAF = LITAF)






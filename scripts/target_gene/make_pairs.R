# =============================
## given the RPM values of two eRNAs (minus and plus), creates the summation of the two 
## usage: Rscript make_pairs.R eRNA_minus eRNA_plus output
### NOT USED 
# =============================
args<-commandArgs(trailingOnly=TRUE)

read_vector <- function(file) {
  vector <- scan(file = file, what = character(), sep = "\t", nlines = 1)
  return(vector)
}

eRNA_minus <- read_vector(args[1])
eRNA_plus <- read_vector(args[2])

#eRNA_minus <- read_vector("./input_files/eRNAs/chr16_11611980_11612400.txt")
#eRNA_plus <- read_vector("./input_files/eRNAs/chr16_11612780_11613560.txt")

eRNA_minus_name <- eRNA_minus[1]
eRNA_plus_name <- eRNA_plus[1]

eRNA_minus <- as.numeric(eRNA_minus[-1])
eRNA_plus <- as.numeric(eRNA_plus[-1])


pair <- (eRNA_minus + eRNA_plus) 
pair <- append(pair, paste(eRNA_minus_name, eRNA_plus_name, sep = "_"), after = 0)

output <- args[3]
#output <- "test.txt"

write(pair, output, sep = "\t", ncolumns=length(pair))

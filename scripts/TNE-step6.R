sample <- "gs://fc-secure-4300ea7f-8e77-4f20-a20b-38417286eaaf/60e7508a-7741-447f-a871-ada0e1daff47/bam2bigwig_workflow/00512029-da84-4de5-8185-3738cd6f5ce1/call-bam2bigwig/PP-3502-SVM24T1.minus.normalized.bw"

file <- strsplit(sample, split="/")
index <- length(file[[1]])

file <- file[[1]][index]

sampleName <- strsplit(file, split="[.]")[[1]][1]

library(tidyverse)


test <- c()
for (x in 1:10) {
  test <- c(test, x)
}


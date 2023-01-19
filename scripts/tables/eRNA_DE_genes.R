# script to massage FINAL.xls file into desired output 

library(data.table)
library(dplyr)

all_data <- fread("scripts/tables/FINAL.xls")

DE_TNE_only <- all_data %>% filter(DE_TNE != FALSE)
copy <- DE_TNE_only

long_DE_TNE <- apply(DE_TNE_only, 1, function(x) {
  
  # genes and ensgs should ALWAYS be the same length 
  genes <- unlist(strsplit(x[["dapg_genes"]], ","))
  ensgs <- unlist(strsplit(x[["dapg_ensgs"]], ","))
  
  if (length(genes) > 1) {
    DE_genes <- list()
    for(i in 1:length(genes)) {

      # mutating 
      x[["dapg_genes"]] <- genes[i]
      x[["dapg_ensgs"]] <- ensgs[i]
      DE_genes[[i]] <- data.frame(as.list(x)) %>% as.data.table()
    }
    temp <- rbindlist(DE_genes)
    print(nrow(temp))
    return(temp)
    
   } else {
    return(data.frame(as.list(x)) %>% as.data.table())
   }
})

dapg_transformed <- rbindlist(long_DE_TNE)

dapg_genes <- dapg_transformed %>% select(TNE, strand, dapg_genes, dapg_ensgs) %>% 
  mutate(dapg_ensgs_no_version = sapply(dapg_ensgs, function(x){
    unlist(strsplit(x, split = ".", fixed = T)[[1]][1])
    }))

no_dups <- unique(dapg_genes)
  
DE_genes <- fread("scripts/tables/DEresult.all.xls") %>% 
  select(V1, padj) %>% mutate(dapg_ensgs_no_version = sapply(V1, function(x){
  unlist(strsplit(x, split = ".", fixed = T)[[1]][1])
}))

test <- merge(DE_genes, no_dups, by = "dapg_ensgs_no_version") %>% select(-c(dapg_ensgs, V1)) 

DE_dapg <- test %>% mutate(strand_text = lapply(strand, function(x) {
  if(x == "-") {"minus"} else {"plus"} 
  })) %>% 
  mutate(TNE = paste(TNE, strand_text, sep = "_")) %>% select(-c(strand, strand_text)) %>% 
  rename("dapg_padj" = "padj")

DE_eRNA <- fread("scripts/tables/DE_eRNA_padj.xls") %>% 
  rename("eRNA_padj" = "padj", "TNE" = "V1")

DE_final_dapg <- merge(DE_eRNA, DE_dapg, by = "TNE")

# TODO make the table with caviar 

##### making the DE genes table into a function 



DE_gene_eRNA_table <- function(df, col) {
  
  long_DE_TNE <- apply(df, 1, function(x) {
    
    # genes and ensgs should ALWAYS be the same length 
    genes <- unlist(strsplit(x[[paste0(col, "_genes")]], ","))
    ensgs <- unlist(strsplit(x[[paste0(col, "_ensgs")]], ","))
    
    if (length(genes) > 1) {
      DE_genes <- list()
      for(i in 1:length(genes)) {
        
        # mutating 
        x[[paste0(col, "_genes")]] <- genes[i]
        x[[paste0(col, "_ensgs")]] <- ensgs[i]
        DE_genes[[i]] <- data.frame(as.list(x)) %>% as.data.table()
      }
      temp <- rbindlist(DE_genes)
      print(nrow(temp))
      return(temp)
      
    } else {
      return(data.frame(as.list(x)) %>% as.data.table())
    }
  })
  
  long_merged <- rbindlist(long_DE_TNE)
  
  vars <- c("TNE", "strand", paste0(col, "_genes"), paste0(col, "_ensgs"))
  
  ## TODO find a way to select the columns with dapg ugh 
  dcol_genes <- long_merged %>% select(all_of(vars)) %>% 
    mutate(dapg_ensgs_no_version = sapply(dapg_ensgs, function(x){
      unlist(strsplit(x, split = ".", fixed = T)[[1]][1])
    })) %>% unique()
  
  DE_genes <- fread("scripts/tables/DEresult.all.xls") %>% 
    select(V1, padj) %>% mutate(dapg_ensgs_no_version = sapply(V1, function(x){
      unlist(strsplit(x, split = ".", fixed = T)[[1]][1])
    }))
  
  # TODO change to column instead of dapg
  test <- merge(DE_genes, col_genes, by = paste0(col,"_ensgs_no_version") ) %>% 
    select(-c(dapg_ensgs, V1)) 
  
  DE_col <- test %>% mutate(strand_text = lapply(strand, function(x) {
    if(x == "-") {"minus"} else {"plus"} 
  })) %>% mutate(TNE = paste(TNE, strand_text, sep = "_")) %>% 
    select(-c(strand, strand_text)) %>% 
    rename("dapg_padj" = "padj")
  
  DE_eRNA <- fread("scripts/tables/DE_eRNA_padj.xls") %>% 
    rename("eRNA_padj" = "padj", "TNE" = "V1")
  
  DE_final<- merge(DE_eRNA, DE_col, by = "TNE")
  
}






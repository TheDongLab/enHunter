# script to massage snp_table.xls file into desired output 

library(data.table)
library(dplyr)

all_data <- fread("input_files/snp_table.xls")

DE_TNE_only <- all_data %>% filter(DE_TNE != FALSE)

# long_DE_TNE <- apply(DE_TNE_only, 1, function(x) {
#   
#   # genes and ensgs should ALWAYS be the same length 
#   genes <- unlist(strsplit(x[["dapg_genes"]], ","))
#   ensgs <- unlist(strsplit(x[["dapg_ensgs"]], ","))
#   
#   if (length(genes) > 1) {
#     DE_genes <- list()
#     for(i in 1:length(genes)) {
# 
#       # mutating 
#       x[["dapg_genes"]] <- genes[i]
#       x[["dapg_ensgs"]] <- ensgs[i]
#       DE_genes[[i]] <- data.frame(as.list(x)) %>% as.data.table()
#     }
#     temp <- rbindlist(DE_genes)
#     print(nrow(temp))
#     return(temp)
#     
#    } else {
#     return(data.frame(as.list(x)) %>% as.data.table())
#    }
# })
# 
# dapg_transformed <- rbindlist(long_DE_TNE)
# 
# # select TNE and genes 
# dapg_genes <- dapg_transformed %>% select(TNE, strand, dapg_genes, dapg_ensgs) %>% 
#   mutate(dapg_ensgs_no_version = sapply(dapg_ensgs, function(x){
#     unlist(strsplit(x, split = ".", fixed = T)[[1]][1])
#     }))
# 
# # unique TNE and gene combinations only 
# no_dups <- unique(dapg_genes)
#   
# DE_genes <- fread("input_files/DE/genes/DEresult.all.xls") %>% 
#   select(V1, padj) %>% mutate(dapg_ensgs_no_version = sapply(V1, function(x){
#   unlist(strsplit(x, split = ".", fixed = T)[[1]][1])
# }))
# 
# test <- merge(DE_genes, no_dups, by = "dapg_ensgs_no_version") %>% select(-c(dapg_ensgs, V1)) 
# 
# DE_dapg <- test %>% mutate(strand_text = lapply(strand, function(x) {
#   if(x == "-") {"minus"} else {"plus"} 
#   })) %>% 
#   mutate(TNE = paste(TNE, strand_text, sep = "_")) %>% select(-c(strand, strand_text)) %>% 
#   rename("dapg_padj" = "padj")
# 
# # table with p values for eRNA
# DE_eRNA <- fread("input_files/DE/eRNA/DE_eRNA_padj.xls") %>% 
#   rename("eRNA_padj" = "padj", "TNE" = "V1")
# 
# DE_final_dapg <- merge(DE_eRNA, DE_dapg, by = "TNE")

# TODO make the table with caviar 

##### making the DE genes table into a function 

# TODO drop the _genes column 

DE_gene_eRNA_table <- function(df, col) {
  
  long_DE_TNE <- apply(df, 1, function(x) {
    
    # genes and ensgs should ALWAYS be the same length 
    #genes <- unlist(strsplit(x[[paste0(col, "_genes")]], ","))
    ensgs <- unlist(strsplit(x[[paste0(col, "_ensgs")]], ","))
    
    if (length(ensgs) > 1) {
      DE_genes <- list()
      for(i in 1:length(ensgs)) {
        
        # mutating 
        #x[[paste0(col, "_genes")]] <- genes[i]
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
  
  #paste0(col, "_genes"),
  vars <- c("TNE", "strand", paste0(col, "_ensgs"))
  
  col_genes <- long_merged %>% select(all_of(vars)) %>% 
    mutate(
      !!sym(paste0(col, "_ensgs_no_version")) := sapply(!!sym(paste0(col, "_ensgs")), function(x){
      unlist(strsplit(x, split = ".", fixed = T)[[1]][1])
    })) %>% unique()
  
  # DE information for all genes (run by Ruifeng)
  DE_genes <- fread("input_files/DE/genes/DEresult.all.xls") %>% 
    select(V1, padj) %>% mutate(!!sym(paste0(col, "_ensgs_no_version")) := sapply(V1, function(x){
      unlist(strsplit(x, split = ".", fixed = T)[[1]][1])
    }))
  
  # merge drops rows in col_genes (matched by ensgs)
  # which are not in DE_genes
  test <- merge(DE_genes, col_genes, by = paste0(col,"_ensgs_no_version") ) %>% 
    select(-c(!!sym(paste0(col, "_ensgs")), V1)) 
  
  DE_col <- test %>% mutate(strand_text = lapply(strand, function(x) {
    if(x == "-") {"minus"} else {"plus"} 
  })) %>% mutate(TNE = paste(TNE, strand_text, sep = "_")) %>% 
    select(-c(strand, strand_text)) %>% 
    rename(!!paste0(col, "_padj") := "padj")
  
  DE_eRNA <- fread("input_files/DE/eRNA/DE_eRNA_padj.xls") %>% 
    rename("eRNA_padj" = "padj", "TNE" = "V1")
  
  DE_final<- merge(DE_eRNA, DE_col, by = "TNE")
  
}

#--------- notes ---------------
# rows with a col_padj (ex. dapg_padj, caviar_padj, prom_padj) of NA 
# are in the DE_genes file but were deemed to have a padj value of NA 

# genes in the col_genes which are not in DE_genes are dropped 
#-------------------------------

test_dapg <- DE_gene_eRNA_table(DE_TNE_only, "dapg") %>% rename(gene = dapg_ensgs_no_version, gene_padj = dapg_padj) 
test_dapg$source <- "dapg"

test_caviar <- DE_gene_eRNA_table(DE_TNE_only, "caviar") %>% rename(gene = caviar_ensgs_no_version, gene_padj = caviar_padj )
test_caviar$source <- "caviar"

# for the Promoter column, the gene and ensgs column don't always have the same length 
# (technically this can happen for all dapg and caviar too)
# potentially just drop the gene column
test_prom <- DE_gene_eRNA_table(DE_TNE_only, "Prom") %>% rename(gene = Prom_ensgs_no_version, gene_padj = Prom_padj)
test_prom$source <- "PCHiC"

# merge the dapg, caviar, prom tables 
# then convert to long format 
#test_dapg_caviar_prom <- merge(merge(test_dapg, test_caviar, by=c("TNE", "eRNA_padj"), all = T, allow.cartesian=TRUE), 
#                               test_prom, by=c("TNE", "eRNA_padj"), all = T, allow.cartesian=TRUE)

# merging in long format instead 
long_eRNA_gene <- rbindlist(list(test_dapg, test_caviar, test_prom))

fwrite(long_eRNA_gene, file = "eRNA_gene_signif.csv")


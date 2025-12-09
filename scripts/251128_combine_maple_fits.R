# This R script takes an input directory of maple fit files and produces a combines them into matrices of effect sizes, standard errors
# and pleiotropy effects.
# Author: Derek Lamb
# Modified: November 26th, 2025

# get arguments from the command line, input and output directories
args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

start_time <- Sys.time()

# setup
setwd("~/Documents/network_discovery_project/")
library(magrittr)

# load file list and gene list 
## the gene list is output by 251121_split_bigbrain.rmd
input_files <- list.files(input_dir, full.names = TRUE)
gene_list <- data.table::fread("data/251128_all_feature_list.tsv", header = FALSE) %>% 
  dplyr::pull(1)

# initialize matrices
tce_matrix <- se_matrix <- correlated_effect_matrix <- matrix(NA, 
                    nrow = length(gene_list), ncol = length(gene_list),
                    dimnames = list(gene_list, gene_list))

# loop over files in directory
for (maple_file in input_files) {
  working_file <- readRDS(maple_file)
  
  # loop over genes in working file
  for (gene_pair in names(working_file)) {
    working_pair <- gene_pair %>% 
      stringr::str_split("_") %>% 
      unlist()
    
    tce_matrix[working_pair[1], working_pair[2]] <- working_file[[gene_pair]]$causal_effect
    se_matrix[working_pair[1], working_pair[2]] <- working_file[[gene_pair]]$cause.se
    correlated_effect_matrix[working_pair[1], working_pair[2]] <- working_file[[gene_pair]]$correlated_pleiotropy_effect
  }
}

# save output
saveRDS(tce_matrix, paste0(output_dir, "combined_maple_tce_matrix.rds"))
saveRDS(se_matrix, paste0(output_dir, "combined_maple_se_matrix.rds"))
saveRDS(correlated_effect_matrix, paste0(output_dir, "combined_maple_correlated_pleiotropy_matrix.rds"))

# finish
end_time <- Sys.time()
print(paste0(
  "Finished processing ",
  length(input_files),
  " MAPLE fit files."
))
print(end_time - start_time)
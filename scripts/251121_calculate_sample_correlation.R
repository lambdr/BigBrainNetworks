# This script computes sample relatedness matrices for gene pairs in the bigbrain data
# Author: Derek Lamb
# Modified: November 21th, 2025

## This script was written assuming that there is a directory of gene pairs, separated into parts (the first argument), 
## along with a directory of summary stats, split by feature, and a single file containing all LD scores. To slightly speed up 
## computation, I removed some redundant calculations in the MAPLE::est_SS() function, included in the `est_ss_code_file` argument.
## An alternative is to simply use MAPLE::est_SS() and skip sourcing this file.

# get arguments from command line - there should be 6 arguments
args <- commandArgs(trailingOnly = TRUE)

part <- args[1] # four digit part indicating which file this script should loop over
gene_pair_dir <- args[2] # directory containing the gene-pair files
sumstats_dir <- args[3] # directory containing summary stats from mmQTL
ldscore_file <- args[4] # path to LD score file
est_ss_code_file <- args[5] # path to 251121_modified_est_SS.R
output_file <- args[6] # location and basename name of output file


# load and source code
library(magrittr)
source(est_ss_code_file)

# load data
df_ldscore <- data.table::fread(ldscore_file)
df_gene_pair <- readRDS(paste0(gene_pair_dir, "/gene_pairs_part_", part, ".rds"))

# initialize lists
omega_list <- list()


# loop
for (i in 1:nrow(df_gene_pair)) {
  # Start timing
  start_time <- Sys.time()
  
  working_pair <- df_gene_pair %>% 
    dplyr::slice(i) %>% 
    dplyr::pull(gene_pair)
  
  print(paste0("Beginning analysis of ", working_pair))
  
  omega_list[[working_pair]] <- est_SS(
    working_pair, 
    df_ldscore = df_ldscore, 
    sumstats_dir = sumstats_dir)
  
  end_time <- Sys.time()
  print(paste0("Finished with ", working_pair))
  print(end_time-start_time)
}


# save data
saveRDS(omega_list, paste0(output_file, "_part_", part, ".rds"))
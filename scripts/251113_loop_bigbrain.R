# This script loops over split parts of the summary stats, with only QTLs extracted and performs read_bigbrain

# Get arguments
args <- commandArgs(trailingOnly = TRUE)

input_dir = args[1]
output_dir = args[2]
start = args[3]
end = args[4]

# Setup
library(magrittr)
setwd("/Users/lambdr/Documents/network_discovery_project/")
source("~/Documents/BigBrainNetworks/scripts/process_bigbrain_functions.R")

# Generate sequence to run over
run_seq <- sprintf('%0.3d', as.numeric(start):as.numeric(end))

# Run over sequence
for (i in run_seq) {
  input_file <- paste0(input_dir, "/part_", i, "_with_header.tsv")
  output_file <- paste0(output_dir, "/part_", i, "_processed_data.rds")
  
  df_matched <- read_bigbrain(
    filename = input_file,
    beta_col = "fixed_beta",
    p_col = "Fixed_P",
    id_col = "variant_id",
    se_col = "fixed_sd",
    p_thresh = 5e-8,
    delim = " " # update for space delim
  )
  
  saveRDS(df_matched, output_file)
  
}
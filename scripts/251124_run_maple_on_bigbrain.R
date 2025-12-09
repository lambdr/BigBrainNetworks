# This script runs MAPLE on the split bigbrain data.
# Author: Derek Lamb
# Modified: November 26th, 2025

## This script was written assuming that there is a directory of RDS files split as in Step 3, separated into parts (the first argument), 
## along with a directory of matched sample-correlation matrices (omega), and the LD matrix. There is also a file of QC'd instruments
## output by the RMD file in Step 3.


# get arguments from command line - there should be 6 arguments
args <- commandArgs(trailingOnly = TRUE)

part <- args[1] # four digit part indicating which file this script should loop over
omega_dir <- args[2] # path to sample structure matrices
split_bigbrain_dir <- args[3] # directory containing summary stats from mmQTL
ld_file <- args[4] # path to LD file
instrument_file <- args[5] # path to instrument list 
output_dir <- args[6] # location to write output

# load for piping reasons
library(magrittr)

# define variables
input_file = paste0(
  split_bigbrain_dir,
  "/bigbrain_part_",
  part,
  ".rds")

# load data
inst_list <- data.table::fread(instrument_file, header=F) %>% 
  dplyr::pull(1)

ld_mat <- readRDS(ld_file)

print(paste0("Finished loading LD matrix."))

omega_list <- readRDS(paste0(omega_dir, 
                             "/gene_correlation_results_part_", 
                             part, ".rds"))


print(paste0("Loading ", input_file))
df_matched <- readRDS(input_file)

# initialize lists
maple_list <- list()

# Loop time
for (i in 1:nrow(df_matched)) {
  # start timer
  start_time <- Sys.time()
  
  # load and process data
  working_data <- df_matched %>% 
    dplyr::slice(i)
  
  working_pair <- working_data$exp_out
  instruments <-  working_data$inst[[1]]
  
  # QC instrument set 
  valid_indic <- (instruments %in% inst_list)
  
  if (sum(valid_indic) < 2) {
    print(paste0("Only ", sum(valid_indic), " valid instruments for ", working_pair, 
                 ". Moving to next gene-pair."))
    next
  }
  
  print(paste0("Beginning analysis for ", working_pair, " exposure-outcome pair with ", 
               sum(valid_indic), " instruments."))
  
  # update with valid instruments
  instruments <- instruments[valid_indic]
  Zscore_1 <- working_data$beta_exp[[1]]/working_data$se_exp[[1]]
  Zscore_1 <- Zscore_1[valid_indic]
  Zscore_2 <- working_data$beta_out[[1]]/working_data$se_out[[1]]
  Zscore_2 <- Zscore_2[valid_indic]
  
  # get working ld and sumstats
  working_sigma <- ld_mat[instruments, instruments]
  
  # get omega matrix and reverse if needed
  paras <- omega_list[[working_data$gene_pair]]
  if (working_data$flipped_pair) paras <- paras[2:1, 2:1]
  
  # run MAPLE
  maple_list[[working_pair]] <- MAPLE::MAPLE(
    Zscore_1 = Zscore_1,
    Zscore_2 = Zscore_2,
    Sigma1in = working_sigma,
    Sigma2in = working_sigma,
    samplen1 = 4656,
    samplen2 = 4656,
    Gibbsnumber = 1000,
    burninproportion = 0.2, 
    pi_beta_shape = 0.5,
    pi_beta_scale = 4.5,
    pi_c_shape = 0.5,
    pi_c_scale = 9.5,
    pi_1_shape = 0.5,
    pi_1_scale = 1.5,
    pi_0_shape = 0.05,
    pi_0_scale = 9.95,
    t1 = paras[1,1],
    t2 = paras[2,2],
    t12 = paras[1,2])
  
  # End time and write status to log
  end_time <- Sys.time()
  print(paste0("Finished with ", working_pair))
  print(end_time-start_time)
}

# save output
saveRDS(maple_list, paste0(output_dir, "bigbrain_maple_fits_part_", part, ".rds"))
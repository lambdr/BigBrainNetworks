# This is an example shell script for Step 2 of the MAPLE analysis workflow

# Read bigbrain
## This assumes you have split the results of Step 1B into chunks and would parallelize running the following:
Rscript code/251113_loop_bigbrain.R data/251113_bigbrain_instruments_split data/251113_processed_qtl_data 0 24


# Join bigbrain
## You would run the following with the input file below being a text file of paths to the outputs from the Read Bigbrain step above
Rscript ~/Documents/BigBrainNetworks/scripts/join_bigbrain.R \
  --inputs ../network_discovery_project/data/list_of_processed_rds_fils.txt \
  --inst ../network_discovery_project/processed_data/251115_full_data_instruments.tsv.gz \
  --result ../network_discovery_project/251115_full_data_results.rds
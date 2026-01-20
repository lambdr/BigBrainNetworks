# This script combines the fits from MAPLE into some summary matrices. It takes three arguments:
## First, the input directory of MAPLE fits
## Second, the output directory location to save the matrices
## Third, a list of all features to extract, this is the output of Step 3

Rscript 251128_combine_maple_fits.R \
	../processed_data/maple_fits/ \
	../results/ \
	../data/251128_all_feature_list.tsv 
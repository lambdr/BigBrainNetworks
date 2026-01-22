#!/bin/bash

# directory for generated LSF scripts
script_dir="260122_run_maple_scripts"
mkdir -p "$script_dir"

# directory for log files
log_dir="/home/lambdr/logs/260122_complete_maple_run/"
mkdir -p "$log_dir"

# define inputs
omega_dir="~/data/260121_gene_correlation_results/"
split_bigbrain_dir="~/data/251124_split_rds_files/"
ld_file="~/data/251126_qc_ld_matrix_for_20711_snps.rds"
instrument_file="~/data/251126_qc_20711_instruments.txt"
output_file="~/results/"

for i in $(seq -w 0001 1000); do
    job_script="${script_dir}/run_maple_${i}.lsf"

    cat <<EOF > "$job_script"
#!/bin/bash
#BSUB -J run_maple_${i}
#BSUB -o ${log_dir}/run_maple_${i}.out
#BSUB -e ${log_dir}/run_maple_${i}.err
#BSUB -q cceb_normal
#BSUB -n 1
#BSUB -M 12G

module load R/4.4

R CMD BATCH \
  --no-save \
  --no-restore \
  "--args ${i} ${omega_dir} ${split_bigbrain_dir} ${ld_file} ${instrument_file} ${output_file}" \
  ~/code/251124_run_maple_on_bigbrain.R \
  ${log_dir}/run_maple_${i}.Rout
EOF
done

echo "Generated 1000 LSF job scripts in $script_dir/"

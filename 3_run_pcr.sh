#!/bin/bash

# Parameters to control the jobs
start_yr=1995
end_yr=2020
main_path="../output/"
tmp_file_imp="$main_path/imp_tmp.csv"
tmp_file_bc="$main_path/bc_tmp.csv"

# Parameters for principal component regressions
signals_file="$main_path/signals_best125_1985.txt" # file with list of signals to use
n_pcs=75 # number of PCs in maximal regression
n_years=10 # number of years for principal components/predictive regs
quantile_prob=0.1 # quantile to form long/short portfolios (0.1 means deciles)

# The sample start year is however many years before start_yr (to make PCs)
sample_start_yr=$(($start_yr-$n_years))

# Paths
params_path="$main_path/impute_ests/" 
output_path="$main_path/pcr_returns/"

# Submit the data prepping job
out=$(Rscript --grid_submit=batch \
    --grid_hold=6726172 \
    --grid_ncpus=12 \
    --grid_mem=200G \
    --grid_email="jmccoy26@gsb.columbia.edu" \
    3a_prep_data_em.R \
        --impute_vec=$signals_file \
        --sample_start_year=$sample_start_yr \
        --sample_end_year=$end_yr \
        --params_path=$params_path \
        --tmp_file_bc=$tmp_file_bc \
        --tmp_file_imp=$tmp_file_imp)

# Extract the job ID for monitoring
jobid=$( echo $out | grep -o -E '[0-9]+' )

# Run the PCR with the new data
# These take up so much memory that I need to do one year at a time
# But they each run relatively quickly
iter=0
for _yr in ` $start_yr $end_yr `; do

    # Simple mean imputations
    Rscript --grid_submit=batch \
        --grid_hold=$jobid \
        --grid_ncpus=10 \
        --grid_mem=250G \
        --grid_SGE_TASK_ID=1-4 \
        --grid_email="jmccoy26@gsb.columbia.edu" \
        3b_pcr.R \
            --signals_keep=$signals_file \
            --data_file=$tmp_file_bc \
            --prefix="pcr_mn_" \
            --out_path=$output_path \
            --iter_year=$_yr \
            --n_yrs=$n_years \
            --quantile_prob=$quantile_prob \
            --n_pcs=$n_pcs

    # EM imputations
    Rscript --grid_submit=batch \
        --grid_hold=$jobid \
        --grid_ncpus=10 \
        --grid_mem=250G \
        --grid_SGE_TASK_ID=1-4 \
        --grid_email="jmccoy26@gsb.columbia.edu" \
        3b_pcr.R \
            --signals_keep=$signals_file \
            --data_file=$tmp_file_imp \
            --prefix="pcr_em_" \
            --out_path=$output_path \
            --iter_year=$_yr \
            --n_yrs=$n_years \
            --quantile_prob=$quantile_prob \
            --n_pcs=$n_pcs

    # Increment for next year
    ((iter+=1))
done


#!/bin/bash

# Parameters to control the jobs
start_yr=1995
end_yr=2020
main_path="../output/"

# Parameters for principal component regressions
signals_file="${main_path}signals_best125_1985.txt" # file with list of signals to use
n_pcs=75 # number of PCs in maximal regression
n_years=10 # number of years for principal components/predictive regs
quantile_prob=0.1 # quantile to form long/short portfolios (0.1 means deciles)

# The sample start year is however many years before start_yr (to make PCs)
sample_start_yr=$(($start_yr-$n_years))

# Paths
params_path="${main_path}impute_ests/" 

# Big datasets
tmp_file_mn="${main_path}/bcsignals/bcsignals_none.csv"
tmp_file_em="${main_path}/bcsignals/bcsignals_em.csv"
tmp_file_ac="${main_path}/bcsignals/bcsignals_availcase.csv"

# Submit the data prepping job
out_mn=$(Rscript --grid_submit=batch \
    --grid_ncpus=12 \
    --grid_mem=200G \
    --grid_email="jmccoy26@gsb.columbia.edu" \
    3a_prep_big_data.R \
        --impute_type="none" \
        --impute_vec=$signals_file \
        --sample_start_year=$sample_start_yr \
        --sample_end_year=$end_yr \
        --params_path=$params_path \
        --bcsignals_filename=$tmp_file_mn)

out_em=$(Rscript --grid_submit=batch \
    --grid_ncpus=12 \
    --grid_mem=200G \
    --grid_email="jmccoy26@gsb.columbia.edu" \
    3a_prep_big_data.R \
        --impute_type="em" \
        --impute_vec=$signals_file \
        --sample_start_year=$sample_start_yr \
        --sample_end_year=$end_yr \
        --params_path=$params_path \
        --bcsignals_filename=$tmp_file_em)

out_ac=$(Rscript --grid_submit=batch \
    --grid_ncpus=12 \
    --grid_mem=200G \
    --grid_email="jmccoy26@gsb.columbia.edu" \
    3a_prep_big_data.R \
        --impute_type="availcase" \
        --impute_vec=$signals_file \
        --sample_start_year=$sample_start_yr \
        --sample_end_year=$end_yr \
        --params_path=$params_path \
        --bcsignals_filename=$tmp_file_ac)

# Extract the job ID for monitoring
jobid_mn=$( echo $out_mn | grep -o -E '[0-9]+' )
jobid_em=$( echo $out_em | grep -o -E '[0-9]+' )
jobid_ac=$( echo $out_ac | grep -o -E '[0-9]+' )

# Run the PCR with the new data
# These take up so much memory that I need to do one year at a time
# But they each run relatively quickly
for _yr in $(eval echo "{$start_yr..$end_yr}"); do
    for forecast in "pca" "spca1" "spca2"; do 
        
        if [ $forecast = "pca" ]; then
            scaled_pca="FALSE"
            scaled_pca_weight="ew"
        elif [ $forecast = "spca1" ]; then
            scaled_pca="TRUE"
            scaled_pca_weight="ew"
        elif [ $forecast = "spca1" ]; then
            scaled_pca="TRUE"
            scaled_pca_weight="vw"
        fi
        
        # Simple mean imputations
        Rscript --grid_submit=batch \
            --grid_hold=$jobid_mn \
            --grid_ncpus=10 \
            --grid_mem=250G \
            --grid_SGE_TASK_ID=1-12 \
            --grid_email="jmccoy26@gsb.columbia.edu" \
            3b_pcr.R \
                --signals_keep=$signals_file \
                --data_file=$tmp_file_mn \
                --out_path="$main_path/pcr_returns/${forecast}_mn/" \
                --scaled_pca=$scaled_pca \
                --scaled_pca_weight=$scaled_pca_weight \
                --iter_year=$_yr \
                --n_yrs=$n_years \
                --quantile_prob=$quantile_prob \
                --n_pcs=$n_pcs
        
        # EM imputations
        Rscript --grid_submit=batch \
            --grid_hold=$jobid_em \
            --grid_ncpus=10 \
            --grid_mem=250G \
            --grid_SGE_TASK_ID=1-12 \
            --grid_email="jmccoy26@gsb.columbia.edu" \
            3b_pcr.R \
                --signals_keep=$signals_file \
                --data_file=$tmp_file_em \
                --out_path="$main_path/pcr_returns/${forecast}_em/" \
                --scaled_pca=$scaled_pca \
                --scaled_pca_weight=$scaled_pca_weight \
                --iter_year=$_yr \
                --n_yrs=$n_years \
                --quantile_prob=$quantile_prob \
                --n_pcs=$n_pcs
        
        # available case
        Rscript --grid_submit=batch \
            --grid_hold=$jobid_ac \
            --grid_ncpus=10 \
            --grid_mem=250G \
            --grid_SGE_TASK_ID=1-12 \
            --grid_email="jmccoy26@gsb.columbia.edu" \
            3b_pcr.R \
                --signals_keep=$signals_file \
                --data_file=$tmp_file_ac \
                --out_path="$main_path/pcr_returns/${forecast}_ac/" \
                --scaled_pca=$scaled_pca \
                --scaled_pca_weight=$scaled_pca_weight \
                --iter_year=$_yr \
                --n_yrs=$n_years \
                --quantile_prob=$quantile_prob \
                --n_pcs=$n_pcs
    done
done


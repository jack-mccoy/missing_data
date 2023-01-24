#!/bin/bash

# User email for error and completion notifications
user_email="$USER@gsb.columbia.edu"

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
scratch_path="/scratch/jpm2223/bcsignals"
tmp_file_mn="${scratch_path}/bcsignals_none.csv"
tmp_file_em="${scratch_path}/bcsignals_em.csv"
tmp_file_ac="${scratch_path}/bcsignals_availcase.csv"
tmp_file_ind="${scratch_path}/bcsignals_ind.csv"
tmp_file_indsize="${scratch_path}/bcsignals_indsize.csv"
tmp_file_lastval="${scratch_path}/bcsignals_lastval.csv"

# Submit the data prepping job
out_mn=$(Rscript --grid_submit=batch \
    --grid_ncpus=12 \
    --grid_mem=200G \
    --grid_email=$user_email \
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
    --grid_email=$user_email \
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
    --grid_email=$user_email \
    3a_prep_big_data.R \
        --impute_type="availcase" \
        --impute_vec=$signals_file \
        --sample_start_year=$sample_start_yr \
        --sample_end_year=$end_yr \
        --params_path=$params_path \
        --bcsignals_filename=$tmp_file_ac)

out_ind=$(Rscript --grid_submit=batch \
    --grid_ncpus=12 \
    --grid_mem=200G \
    --grid_email=$user_email \
    3a_prep_big_data.R \
        --impute_type="ind" \
        --impute_vec=$signals_file \
        --sample_start_year=$sample_start_yr \
        --sample_end_year=$end_yr \
        --params_path=$params_path \
        --bcsignals_filename=$tmp_file_ind)

out_indsize=$(Rscript --grid_submit=batch \
    --grid_ncpus=12 \
    --grid_mem=200G \
    --grid_email=$user_email \
    3a_prep_big_data.R \
        --impute_type="indsize" \
        --impute_vec=$signals_file \
        --sample_start_year=$sample_start_yr \
        --sample_end_year=$end_yr \
        --params_path=$params_path \
        --bcsignals_filename=$tmp_file_indsize)

out_lastval=$(Rscript --grid_submit=batch \
    --grid_ncpus=12 \
    --grid_mem=200G \
    --grid_email=$user_email \
    3a_prep_big_data.R \
        --impute_type="lastval" \
        --impute_vec=$signals_file \
        --sample_start_year=$sample_start_yr \
        --sample_end_year=$end_yr \
        --params_path=$params_path \
        --bcsignals_filename=$tmp_file_lastval)

# Extract the job ID for monitoring
jobid_mn=$( echo $out_mn | grep -o -E '[0-9]+' )
jobid_em=$( echo $out_em | grep -o -E '[0-9]+' )
jobid_ac=$( echo $out_ac | grep -o -E '[0-9]+' )
jobid_ind=$( echo $out_ind | grep -o -E '[0-9]+' )
jobid_indsize=$( echo $out_indsize | grep -o -E '[0-9]+' )
jobid_lastval=$( echo $out_lastval | grep -o -E '[0-9]+' )

# Run the PCR with the new data
# These take up so much memory that I need to do one year at a time
# But they each run relatively quickly
for _yr in $(eval echo "{$start_yr..$end_yr}"); do
    for imp in "mn" "em" "ac" "ind" "indsize" "lastval"; do
        prep_jobid=jobid_${imp}
        data_file=tmp_file_${imp}
        for forecast in "pca" "spca1" "spca2"; do 

            out_path="$main_path/pcr_returns/${forecast}_${imp}/"
            mkdir -p $out_path # -p says make only if doesn't exist
            
            if [ $forecast = "pca" ]; then
                scaled_pca="FALSE"
                scaled_pca_weight="ew"
            elif [ $forecast = "spca1" ]; then
                scaled_pca="TRUE"
                scaled_pca_weight="ew"
            elif [ $forecast = "spca2" ]; then
                scaled_pca="TRUE"
                scaled_pca_weight="vw"
            fi
            
            # Run the PCR for the given imputation and forecast type
            Rscript --grid_submit=batch \
                --grid_hold=${!prep_jobid} \
                --grid_ncpus=12 \
                --grid_mem=250G \
                --grid_SGE_TASK_ID=1-12 \
                --grid_email=$user_email \
                3b_pcr.R \
                    --signals_keep=$signals_file \
                    --data_file=${!data_file} \
                    --out_path=$out_path \
                    --scaled_pca=$scaled_pca \
                    --scaled_pca_weight=$scaled_pca_weight \
                    --iter_year=$_yr \
                    --n_yrs=$n_years \
                    --quantile_prob=$quantile_prob \
                    --n_pcs=$n_pcs
        done
    done
done


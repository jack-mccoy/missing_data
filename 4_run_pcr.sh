#!/bin/bash

# User email for error and completion notifications
user_email="$USER@gsb.columbia.edu"

# Parameters to control the jobs
start_yr=1995
end_yr=2021

# Parameters for principal component regressions
n_pcs=160 # number of PCs in maximal regression
skip_n=10 # skipper in sequence
n_years=10 # number of years for principal components/predictive regs
quantile_prob=0.1 # quantile to form long/short portfolios (0.1 means deciles)

# The sample start year is however many years before start_yr (to make PCs)
sample_start_yr=$(($start_yr-$n_years))

# Run the PCR with the new data
# These take up so much memory that I need to do one year at a time
# But they each run relatively quickly
for _yr in $(eval echo "{$start_yr..$end_yr}"); do
    for imp in "none" "em" "bllp6"; do
        for forecast in "pca" "spca1"; do 
            
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
                --grid_ncpus=12 \
                --grid_mem=250G \
                --grid_SGE_TASK_ID=1-12 \
                --grid_email=$user_email \
                4a_pcr.R \
                    --signal_file="bcsignals_${imp}.csv" \
                    --scaled_pca=$scaled_pca \
                    --scaled_pca_weight=$scaled_pca_weight \
                    --iter_year=$_yr \
                    --n_yrs=$n_years \
                    --quantile_prob=$quantile_prob \
                    --n_pcs=$n_pcs \
                    --skip_n=$skip_n
        done
    done
done


#!/bin/bash

# Parameters to control the jobs
start_yr=1990
end_yr=2018
tmp_file_imp="imp_tmp.csv"
tmp_file_bc="bc_tmp.csv"

# Parameters for principal component regressions
signals_file="signals.txt" # file with list of signals to use
n_pcs=75 # number of PCs in maximal regression
n_years=10 # number of years for principal components/predictive regs
quantile_prob=0.1 # quantile to form long/short portfolios (0.1 means deciles)

# Paths
params_path="/user/jpm2223/Documents/missing_data/output/impute_ests/" 
output_path="/user/jpm2223/Documents/missing_data/output/pcr_returns/"

# Note the number of months for when we submit PCR array job
n_mons=$((( $end_yr - $start_yr + 1 ) * 12 )) 

## Submit the data prepping job
#out=$(Rscript --grid_submit=batch \
#    --grid_ncpus=20 \
#    --grid_mem=200G \
#    --grid_email="jmccoy26@gsb.columbia.edu" \
#    3a_prep_data_em.R \
#        --params_path=$params_path \
#        --tmp_file_bc=$tmp_file_bc \
#        --tmp_file_imp=$tmp_file_imp)
#
## Extract the job ID for monitoring
#jobid=$( echo $out | grep -o -E '[0-9]+' )
#
## Monitor the data prepping job to wait until the data is ready, the run PCR
#prep_count=$( qstat | grep -i $jobid | wc -l )
#while [ $prep_count -gt 0 ]; do
#    sleep 1m # check back every minute
#    prep_count=$( qstat | grep -i $jobid | wc -l )
#done
#
## Just to keep track in log
#echo "Data is ready"

# Run the PCR with the new data
# These take up so much memory that I need to do one year at a time
# But they each run relatively quickly
iter=0
for _yr in `seq $start_yr $end_yr`; do
    Rscript --grid_submit=batch \
        --grid_ncpus=12 \
        --grid_mem=1000G \
        --grid_SGE_TASK_ID=1-12 \
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

    # Wait for year to finish
    #jobid=$( echo $out | grep -o -E '[0-9]+' | head -1 )
    #echo $jobid
    #job_count=$( qstat | grep -i "$jobid" | wc -l )
    #while [ $job_count -gt 0 ]; do
    #    sleep 1m # check back every minute
    #    job_count=$( qstat | grep -i "$jobid" | wc -l ) 
    #done
    
    # Increment for next year
    ((iter+=1))
done


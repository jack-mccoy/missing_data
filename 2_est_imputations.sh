#!/bin/bash

# They are all the same so write it here
user_email="$USER@gsb.columbia.edu"
signals_list="../data/signals_best158_1985.txt"

if [ $1 = "em" ]; then
    Rscript --grid_submit=batch \
        --grid_SGE_TASK_ID=1985-2020 \
        --grid_mem=100G \
        --grid_ncpus=12 \
        --grid_email=$user_email \
        2a_ar1_em_est.R \
            --impute_vec=$signals_list \
            --em_type="$2" \
            --maxiter=10000 \
            --out_path="/scratch/jpm2223/" \
            --ar1_sample_length=5 \
            --tol=1e-4 \
            -f
fi

if [ $1 = "ppca" ]; then
    Rscript --grid_submit=batch \
        --grid_mem=800G \
        --grid_ncpus=12 \
        --grid_email=$user_email\
        2c_ppca_est.R \
            --start_yr=1985 \
            --end_yr=2020 \
            --impute_vec=$signals_list \
            --maxiter=100000 \
            --out_path="/scratch/jpm2223/bcsignals/" \
            --n_pcs=$2
fi

if [ $1 = "bllp" ]; then
    Rscript --grid_submit=batch \
        --grid_mem=100G \
        --grid_ncpus=12 \
        --grid_email=$user_email \
        2d_bllp_imp.R \
            --num_PCs=6 \
            --out_path="/scratch/jpm2223/" \
            --out_name="bcsignals_bllp.csv" \
            --impute_vec=$signals_list
fi

if [ $1 = "adhoc" ]; then
    Rscript --grid_submit=batch \
        --grid_mem=500G \
        --grid_ncpus=12 \
        --grid_email=$user_email \
        2e_adhoc_imps.R \
            --impute_type="$2" \
            --impute_vec=$signals_list \
            --data_path="/scratch/jpm2223/bcsignals/" \
            --sample_start_year=1985 \
            --sample_end_year=2020 \
            --cores_frac=1
fi



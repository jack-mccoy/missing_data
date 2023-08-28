#!/bin/bash

if [ $1 = "em" ]; then
    Rscript --grid_submit=batch \
        --grid_SGE_TASK_ID=1985-2020 \
        --grid_mem=50G \
        --grid_ncpus=12 \
        --grid_email="jmccoy26@gsb.columbia.edu" \
        2a_ar1_em_est.R \
            --impute_vec="../output/signals_best125_1985.txt" \
            --em_type="$2" \
            --maxiter=10000 \
            --out_path="/scratch/jpm2223/" \
            --ar1_sample_length=5 \
            --tol=1e-4 \
            -f
fi

if [ $1 = "ppca" ]; then
    Rscript --grid_submit=batch \
        --grid_mem=500G \
        --grid_ncpus=12 \
        --grid_email="$USER@gsb.columbia.edu" \
        2c_ppca_est.R \
            --start_yr=1985 \
            --end_yr=2020 \
            --impute_vec="../data/signals_best125_1985.txt" \
            --maxiter=100000 \
            --out_path="/scratch/jpm2223/bcsignals/" \
            --n_pcs=$2
fi

if [ $1 = "bllp" ]; then
    Rscript --grid_submit=batch \
        --grid_mem=100G \
        --grid_ncpus=12 \
        --grid_email="$USER@gsb.columbia.edu" \
        2d_bllp_imp.R \
            --num_PCs=6 \
            --out_path="/scratch/jpm2223/" \
            --out_name="bcsignals_bllp.csv" \
            --impute_vec="../output/signals_best125_1985.txt"
fi


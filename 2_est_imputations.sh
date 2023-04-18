#!/bin/bash

if [ $1 = "em" ]; then
    Rscript --grid_submit=batch \
        --grid_SGE_TASK_ID=1987-2020 \
        --grid_mem=40G \
        --grid_ncpus=12 \
        --grid_email="jmccoy26@gsb.columbia.edu" \
        2a_em_est.R \
            --impute_vec="../output/signals_best125_1985.txt" \
            --maxiter=10000 \
            --out_path="../output/impute_ests/" \
            --tol=1e-4 \
            --boxcox \
            -f
fi

if [ $1 = "ppca" ]; then
    Rscript --grid_submit=batch \
        --grid_mem=500G \
        --grid_ncpus=12 \
        --grid_email="$USER@gsb.columbia.edu" \
        2b_ppca_est.R \
            --start_yr=1985 \
            --end_yr=2020 \
            --impute_vec="../output/signals_best125_1985.txt" \
            --maxiter=100000 \
            --out_path="/scratch/jpm2223/bcsignals/" \
            --n_pcs=$2
fi


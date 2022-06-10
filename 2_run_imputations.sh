#!/bin/bash

Rscript --grid_submit=batch \
    --grid_SGE_TASK_ID=1980-2018 \
    --grid_mem=1000G \
    --grid_ncpus=39 \
    --grid_email="jmccoy26@gsb.columbia.edu" \
    2a_em_est.R \
        --impute_vec="signals.txt" \
        --maxiter=10000 \
        --out_path="/user/jpm2223/Documents/missing_data/output/impute_ests/" \
        --tol=1e-4 \
        -f # force the convergence
        # Note that --boxcox IS NOT SPECIFIED right now. Need to specify explicitly now

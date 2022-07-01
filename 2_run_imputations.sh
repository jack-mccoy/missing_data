#!/bin/bash

Rscript --grid_submit=batch \
    --grid_SGE_TASK_ID=1980-2020 \
    --grid_mem=40G \
    --grid_ncpus=2 \
    --grid_email="jmccoy26@gsb.columbia.edu" \
    2a_em_est.R \
        --impute_vec="signals.txt" \
        --maxiter=10000 \
        --out_path="../output/impute_ests/" \
        --tol=1e-4 \
        --boxcox \
        -f # force the convergence

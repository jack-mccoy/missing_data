#!/bin/bash

Rscript --grid_submit=batch \
    --grid_SGE_TASK_ID=1985-2020 \
    --grid_mem=40G \
    --grid_ncpus=12 \
    --grid_email="jmccoy26@gsb.columbia.edu" \
    2a_em_est.R \
        --impute_vec="../output_best157/signals_best157.txt" \
        --maxiter=10000 \
        --out_path="../output_best157/impute_ests/" \
        --tol=1e-4 \
        --boxcox \
        -f # force the convergence

#!/bin/bash

Rscript --grid_submit=batch \
    --grid_SGE_TASK_ID=1985-1986 \
    --grid_mem=40G \
    --grid_ncpus=12 \
    --grid_email="jmccoy26@gsb.columbia.edu" \
    2a_em_est.R \
        --impute_vec="../output_best125_1985_fixedmns/signals_best125_1985.txt" \
        --maxiter=10000 \
        --out_path="../output_best125_1985_fixedmns/impute_ests/" \
        --tol=1e-4 \
        --boxcox \
        -f # force the convergence

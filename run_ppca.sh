#/bin/bash

Rscript --grid_submit=batch \
    --grid_SGE_TASK_ID=1985-2020 \
    --grid_mem=200G \
    --grid_ncpus=12 \
    --grid_email="$USER@gsb.columbia.edu" \
    ppca_test.R \
        --impute_vec="../output/signals_best125_1985.txt" \
        --maxiter=10000 \
        --out_path="/scratch/jpm2223/bcsignals/" \
        --n_pcs=40



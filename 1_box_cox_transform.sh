#!/bin/bash

Rscript \
    --grid_submit=batch \
    --grid_mem=500G \
    --grid_email="$USER@gsb.columbia.edu" \
    --grid_ncpus=10 \
    1a_box_cox_transform.R \
        --cores_frac=1 \
        --impute_vec="../data/signals_best125_1985.txt" \
        --data_path="../data/" \
        --out_path="/scratch/jpm2223/bcsignals/"

#!/bin/bash

Rscript \
    --grid_submit=batch \
    --grid_mem=500G \
    --grid_email="$USER@gsb.columbia.edu" \
    --grid_ncpus=10 \
    1a_box_cox_transform.R \
        --cores_frac=1

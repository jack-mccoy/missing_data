#!/bin/bash

# Submit job
Rscript --grid_submit=batch \
    --grid_mem=20G \
    1a_select_signals.R \
        --sample_start_year=1980 \
        --sample_end_year=2020 \
        --selection_year=1980 \
        --selection_month=6 \
        -n 100

#!/bin/bash

# Submit job
Rscript --grid_submit=batch \
    --grid_email="jmccoy26@gsb.columbia.edu" \
    --grid_mem=20G \
    1a_select_signals.R

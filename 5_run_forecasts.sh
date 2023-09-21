#/bin/bash

# Array of model names
mods=("keras1" "keras2" "keras3" "keras4" "keras5" "lightgbm" "lm" "pcr" "spcr") 
# Array of imputation types
imps=("none" "em" "emar1" "bllp6" "ppca6" "ppca10" "ppca40" "ppca55" "ind" "indsize" "lastval")
firmsets=("micro" "small" "big" "all")

# Start and end dates
yearm_begin="1995-06"
yearm_end="2020-06"

# Paths
output_folder="auto"

# Loop through each model
for mod in ${mods[@]}; do
    for imp in ${imps[@]}; do
        for firmset in ${firmsets[@]}; do
            Rscript \
                --grid_submit=batch \
                --grid_email="$USER@gsb.columbia.edu" \
                --grid_mem=50G \
                5a_one_forecast.R \
                    --model=$mod \
                    --signal_file="bcsignals_${imp}.csv" \
                    --firmset=$firmset \
                    --output_folder=$output_folder \
                    --yearm_begin=$yearm_begin \
                    --yearm_end=$yearm_end
        done
    done
done


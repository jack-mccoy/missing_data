# Repo for missing data imputations

## Setting up file paths 
User should create a file `FILEPATHS.R` that contains a `list` named `FILEPATHS`, with entries `(data_path, out_path, signal_list)`, i.e.

```r
FILEPATHS <- list(
    data_path = "../data/",
    out_path = "../output/",
    signal_list = "../data/signals_best125_1985.txt"
)
```

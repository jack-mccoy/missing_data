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

Note that the scripts will automatically generate data and output folders.
It will also download the Chen-Zimmerman and other datasets, and output fully
imputed datasets. The user should make sure the specified `data_path` and
`out_path` folders are large enough to store files that will take up 
several dozen GB.

## Running the code
The code generally requires a HPC set up to run the imputations and forecasting
algorithms. The various `.sh` files are wrappers that are tailored to the CBS
grid, but hopefully provide examples.

Most of the scripts accept trailing arguments using the `optparse` 
package in `R`.




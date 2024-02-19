---
contributors:
  - Andrew Chen
  - Jack McCoy
title: Replication Documentation for 
    "Missing Values Handling For Machine Learning Portfolios"
author: Andrew Chen, Jack McCoy
geometry: margin=0.5in
---

## Overview

The code in this replication package 

* Downloads return data from CRSP via WRDS and signal data from the Chen-Zimmerman 
    Open Source Asset Pricing website;
* Provides figures and tables describing missingness in the dataset (Section 2
    in the paper);
* Applies a Box-Cox transformation to the signals to produce `bcsignals_none.csv`
    (Section 3.1);
* Imputes the transformed data to construct `bcsignals_IMP.csv`, 
    where `IMP` is one of the imputation methods explored in the paper;
    (Sections 3.2 and 5.1);
* For various imputed datasets, produces the principal component (PC) regression
    long-short portfolios for different numbers of PCs (Section 4 in the paper); 
* Constructs long-short portfolios for hypertuned machine learning methods
    (Section 5 in the paper); and
* Produces figures and tables for the paper.

## Data Availability and Provenance Statements

The raw data comes from three main sources:

* The Chen-Zimmerman Open Source Asset Pricing dataset, the construction of which
    is described in Chen and Zimmerman (2022), 
* CRSP/Compustat as hosted on WRDS, and
* Ken French's data repository.
    - Market equity decile breakpoints are downloaded directly from his website
    - Fama and French (2015) and Carhart (1997) monthly factor portfolio data
        are downloaded via WRDS

The Chen-Zimmerman data is used for the file `signed_predictors_dl_wide.csv`,
which has 212 monthly firm-level signal replications replicated in Chen and
Zimmerman (2022). We also use the file `SignalDoc.csv` to filter our signal
set to continuous predictors. Both datasets are available publicly at the 
[Open Source Asset Pricing website](https://www.openassetpricing.com/). 

The CRSP/Compustat data is available to any user with the necessary WRDS 
account. It is downloaded in `0a_download_data.R`. We use this for monthly
firm-level return and market cap data, as well as several stock signals that are
immediately constructed from the CRSP data, hence why they are not already
available in the Chen-Zimmerman dataset.

For users that are not able to access CRSP, we provide the pseudo data file
`crsp_data_dummy.csv`. This matches the structure of the downloaded and cleaned
dataset in the paper `crsp_data.csv`, but it *will not* replicate the actual
results. It only exists to allow users without access to CRSP run and test
the paper's code.

### Statement about Rights

- [X] I certify that the author(s) of the manuscript have legitimate access to and permission to use the data used in this manuscript. 
- [X] I certify that the author(s) of the manuscript have documented permission to redistribute/publish the data contained within this replication package. Appropriate permission are documented in the [LICENSE.txt](LICENSE.txt) file.

### Summary of Availability

- [ ] All data **are** publicly available.
- [X] Some data **cannot be made** publicly available.
- [ ] **No data can be made** publicly available.

### Public use data sourced from elsewhere and provided

Data on firm-level monthly stock signals are available via the 
[Open Source Asset Pricing website](https://www.openassetpricing.com/).
We use the file `signed_predictors_dl_wide.csv` in the August 2023 release. 
We also use the file `SignalDoc.csv` to filter the raw signals data.
Data can be downloaded from 
https://drive.google.com/drive/folders/1EP6oEabyZRamveGNyzYU0u6qJ-N43Qfq.
Select "Firm Level Characteristics" then "Full Sets" for firm-level predictors.
A copy of the data is provided as part of this archive.
The data are in the public domain.

Datafile:  `signed_predictors_dl_wide.csv`; `SignalDoc.csv`

### Proprietary data

Monthly firm-level stock return data from CRSP/Compustat are downloaded via
WRDS. Since these require a subscription, we are not able to provide them in
the replication packet. For those users with a WRDS login and CRSP access,
the CRSP file we use can be downloaded by running `0a_download_data.R`. It is helpful
to have a `~/.pgpass` set up to access WRDS's Postgres database, but there is
also an interactive username and password entry option for users who do not
have this set up. Note that the script also downloads and stores the signals data
from the Open Source Asset Pricing website.

For those users without a WRDS subscription, we provide the pseudo data file
`crsp_data_dummy.csv`. This *will not* replicate the paper's results, but will
allow you to run the code. You will either need to replace the filename
`crsp_data.csv` wherever it appears in the code with `crsp_data_dummy.csv`.

### Example for preliminary code during the editorial process

Code for data cleaning and analysis is provided as part of the replication package. 
It is available at https://github.com/jack-mccoy/missing_data for review.
It has also been uploaded to the Journal of Financial Economics data repository
for the accepted manuscript.

## Dataset list

All of the below datasets are stored in `$data_path/raw/`, where `$data_path`
is the path in to the list entry of `FILEPATHS.R`.

| Data file | Source | Notes    |Provided |
|-----------|--------|----------|---------|
| `signed_predictors_dl_wide.csv` | Chen and Zimmerman (2022) | Signal acronyms described in `SignalDoc.csv` | Yes |
| `signed_predictors_dl_wide_filtered.csv` | Chen and Zimmerman (2022) | Same as unfiltered file, but with only common shares on major exchanges. Provided for users who are not able to run `0a_download_data.R`, which also filters data, because of lack of WRDS login. | Yes |
| `SignalDoc.csv` | Chen and Zimmerman (2022) | Describes signals contained in `signed_predictors_dl_wide.csv` | Yes |
| `ME_Breakpoints_20_50_clean.csv` | Kenneth French | Contains 20th and 50th percentile monthly NYSE breakpoints, cleaned from Ken French's data repository. | Yes |
| `ff5_factors.csv` | Kenneth French via WRDS WRDS | Contains monthly returns for Fama and French (2015) five-factor portfolios, plus Carhart (1997) momentum factor. | Yes |
| `crsp_data_dummy.csv` | CRSP/Compustat via WRDS | Contains pseudo data for permno-month level returns, plus derived signals which cannot be shared in `signed_predictors_dl_wide.csv` | Yes |
| `crsp_data.csv` | CRSP/Compustat via WRDS | Contains actual permno-month level return data, and derived signals which cannot be shared in `signed_predictors_dl_wide.csv` | No |

The variable names of the signals are provided in `SignalDoc.csv`. Additional
variables are described in the table below. When the same name appears in
a `bcsignals_IMP.csv` file, the units of the variable are then Box-Cox transformed
plus scaled to unit variance, so they differ from the original documentation.

| Variable | Description |
|----------|-------------|
|`permno`| CRSP share-level PERMNO identifier|
|`yyyymm`| `zoo::yearmon`-formatted year-month of observation. Note that this is stored as a float. |
|`hsiccd`| Firm SIC code from CRSP.|
|`ret`| Monthly percent return out of 100 from CRSP. |
|`bh1m`| The one-month ahead return. That is, if `yyyymm` is June 2020, `bh1m` will be the return for July 2020.|
|`me`| Market equity from CRSP |
|`mktrf`| Fama-French monthly market return less RF rate.|
|`smb`| Fama-French monthly size factor portfolio return.|
|`hml`| Fama-French monthly value factor portfolio return.|
|`rmw`| Fama-French monthly profitability factor portfolio return.|
|`cma`| Fama-French monthly investment factor portfolio return.|
|`umd`| Carhart monthly momentum factor portfolio return.|
|`p20`| Monthly NYSE 20th percentile of market cap breakpoint from Ken French's website.|
|`p50`| Monthly NYSE 50th percentile of market cap breakpoint from Ken French's website.|

## Computational requirements

### Software Requirements

- [X] The replication package contains one or more programs to install all 
    dependencies and set up the necessary directory structure. You can 
    run this with `Rscript 00_install_packages.R`.

- The Fortran shared library `mvn_emf.so` was compiled with GNU Fortran 8.3.0
    from the source file `mvn_emf.f90`. Linux users can re-compile the 
    file with `compile_mvn_emf.sh`
- R code was run on R version 4.1.2. Packages used include:
    - `car` (3.1-0)
    - `data.table` (1.14.2)
    - `doParallel` (1.0.17)
    - `dplyr` (1.0.10)
    - `extrafont` (0.18)
    - `foreach` (1.5.2)
    - `getPass` (0.2-2)
    - `ggplot2` (3.3.6)
    - `googledrive` (2.0.0)
    - `gridExtra` (2.3)
    - `janitor` (2.2.0)
    - `kableExtra` (1.3.4)
    - `keras` (2.11.0)
    - `latex2exp` (0.9.5)
    - `lightgbm` (3.3.5)
    - `lubridate` (1.8.0)
    - `optparse` (1.7.3)
    - `pcaMethods` (1.79.1)
    - `pracma` (2.4.2)
    - `ranger` (0.14.1)
    - `RPostgres` (1.4.4)
    - `stringr` (1.4.1)
    - `tensorflow` (2.11.0)
    - `tidyr` (1.2.1)
    - `tidyverse` (1.3.2)
    - `xtable` (1.8-4)
    - `zoo` (1.8-11)

Portions of the code use bash scripting, which may require Linux.

### Controlled Randomness

- [X] Random seed is set at line 482 of `4a_one_forecast.R`.
- [ ] No Pseudo random generator is used in the analysis described here.

### Memory, Runtime, Storage Requirements

#### Summary

The code is structured in favor of storing nearly all different imputations and 
forecasts at the permo-month level. This was done to minimize the amount of 
pre-processing within each script, at the cost of a large amount of storage
requirements.

Approximate time needed to reproduce the analyses on a standard 2023 desktop machine:

- [ ] <10 minutes
- [ ] 10-60 minutes
- [ ] 1-2 hours
- [ ] 2-8 hours
- [ ] 8-24 hours
- [ ] 1-3 days
- [ ] 3-14 days
- [X] > 14 days
- [X] Not feasible to run on a desktop machine, as described below.

Running the code on a desktop machine is not recommended. Only small chunks
will be feasible.

Approximate storage space needed:

- [ ] < 25 MBytes
- [ ] 25 MB - 250 MB
- [ ] 250 MB - 2 GB
- [ ] 2 GB - 25 GB
- [ ] 25 GB - 250 GB
- [X] > 250 GB

The `raw` data takes up about 22G. The `bcsignals` data with different imputations
takes up about 100G. The PCA forecast results take up 113G (because each stock
has forecasts made separately for k PCs, with k=1,...,J), and the other
forecasting results take up 17G.

#### Details

The code was last run on the Columbia Business School computing grid. 
The grid consists of two types of compute nodes, and batch jobs were allocated
depending on availability. 

The first type consists of 
**28 nodes with 64 cores (2 microprocessors X 32core; 2nd Gen AMD CPUs) and 1TB of RAM**.
The second type consists of 
**30 nodes with 64 cores (2 microprocessors X 32core; 4th Gen AMD CPUs) and 1.5TB of RAM**.

Computation took approximately **72 hours** by running imputations in parallel,
followed by running forecasting jobs in parallel. 

The data download and Box-Cox transformation parts of the code are feasible to
run on a desktop machine.

For the EM imputations, each month is imputed in parallel on one CPU
from the compute nodes above. The EM imputations can take roughly 10 minutes
to a day for one month to converge. Therefore, running
`2a_ar1_em_est.R` on a desktop machine is only feasible for a few select months.
Since each month was run in parallel, this required 
37 years x 12 months = 444 CPUs. The length of time to run the imputations on
a comparable computing cluster is, in practice, the max of runtimes across the
different imputation months, so roughly 24 hours.

The most memory-taxing part of the code is the PCR regressions (`3a_pcr.R`). 
Each forecasting month is run in parallel, and within each forecasting month,
we loop through the different sets of k PCs. Consequently, this takes about
250GB of virtual memory, since R regression output are stored simultaneously. 
However, each forecasting month job runs in less than 5 minutes.
Each job is submitted with 12 CPUs requested, but Columbia maxes out the number
of requested memory per user at 2TB. Therefore, there are 8x12=96 CPUs running
in parallel at any given time. The PCR code therefore takes about 5-6 hours to
run on a comparable computing cluster.

Therefore, to run the main result of the paper (Figure 2), it will take about
30 hours on a computing cluster to (i) download and clean the data, (ii) run
the EM imputations, and (iii) run the PC regressions. Running the 
`5e_pcr_plots.R` code takes under 10 minutes once the data is ready.

## Description of programs/code

The code is structured as a flat directory. Files are prefixed in the order that
they should be run. The core R scripts all have a numeric and alphabetic prefix.
Many of the R scripts take command line arguments to allow
some flexibility for different specifications or parallel computations.
If a file has specific arguments that need to be passed to it, a bash script
is provided that shows the necessary command line arguments.

The shell scripts are provided as examples to show how the files were run on the
Columbia Business School computing grid, but they will likely need to be adapted 
to your computing environment. The grid submission syntax for the CBS grid is

    `Rscript [grid options] mycode.R [script options]` 

The `[script options]` are specific to the program and parsed with the package
`optparse` in R.
So the only thing a user should need to change is the material in `[grid options]`.

The options that we use are the following:

* `--grid_submit=batch`: Run the job in batch mode on the grid.
* `--grid_SGE_TASK_ID=id1-idN`: Runs embarassingly parallel jobs that with
    evironment variable `SGE_TASK_ID=id1,...,idN` passed to R. This is used
    to run imputation years in parallel. 
* `--grid_mem=XG`: requests `X` GB of memory for the job. When paired with
    `--grid_SGE_TASK_ID`, `X` is the amount of memory *per task*.
* `--grid_ncpus=N`: The number `N` of CPUs per job/task.
* `--grid_email=EMAIL`: The user's email for sending job metadata.

The scripts are structured as follows:

* The script `FILEPATHS.R` provides the directory settings for the data and
    output. The user should replace `data_path` with their desired storage
    path for the data and `out_path` with their desired storage path for
    tables and figures. To replicate the paper, leave the `signal_list`
    entry as `YOUR DATA PATH/signals_all_1985.txt`. This option was left
    for earlier drafts that used a smaller set of signals.
* Programs prefixed with 0* will download and prepare the necessary data.
    They should be run, in order, as 

    `Rscript 0a_download_data.R`
    `Rscript 0b_select_signals.R`

    The latter script selects the set of "good" signals to use in imputations. 
    In the final draft of the paper, this is all continuous signals without 
    any totally incomplete months. The first script also constructs the
    necessary data directories.
* The program `1a_box_cox_transform.R` applies the modified Box-Cox transformation
    of Hawkins and Weisberg (2017) to the selected signals. If you are using
    a computing grid with the same syntax as CBS, it can be run with 
    `./1_box_cox_transform.sh`. Alternatively, it can be run with
    `Rscript 1a_box_cox_transform.R --cores_frac=1`.
    This script will create the directory `YOUR DATA PATH/bcsignals/`.
* The programs prefixed `2X` run the different imputation methods. The details
    of the most arguments are shown in `2_run_imputations.sh`. This latter
    script allows the user to choose different imputations with command line
    arguments. They are listed here:
    * Regular EM imputations: `./2_run_imputations em regular`
        runs the script `2a_ar1_em_est.R` with argument for regular EM imputations
        * After running this, run `Rscript 2b_bind_em_years.R --em_type=regular`
            to combine all the separate years together
    * EM imputations on AR(1) residuals: `./2_run_imputations em ar1`
        runs the script `2a_ar1_em_est.R` with argument for AR(1) EM imputations
        * After running this, run `Rscript 2b_bind_em_years.R --em_type=ar1`
            to combine all the separate years together
    * Probabilistic PCA (PPCA) imputations with N principal components: 
        `./2_run_imputations ppca N` runs the script `2c_ppca_est.R` with
        argument for N latent factors. It outputs the file
        `YOUR DATA PATH/bcsignals/bcsignals_ppcaN.csv`
    * Bryzgalova, Lerner, Lettau, and Pelger (BLLP, 2023) imputations with N factors:
        `./2_run_imputations bllp N` runs the script `2d_bllp_imp.R` with 
        argument for N latent factors. Outputs the file
        `YOUR DATA PATH/bcsignals/bcsignals_bllpN.csv`
    * Industry-by-size decile mean imputations: 
        `./2_run_imputations adhoc indsize`
        runs the script `2e_adhoc_imps.R` with flag for industry-by-size
        imputations. Outputs the file
        `YOUR DATA PATH/bcsignals/bcsignals_indsize.csv`
    * Last available value imputations:
        `./2_run_imputations adhoc lastval`
        runs the script `2e_adhoc_imps.R` with flag for last available value
        imputations. Outputs the file
        `YOUR DATA PATH/bcsignals/bcsignals_lastval.csv`
* The script `3a_pcr.R` runs the principal component regressions and forms
    long short portfolios based on predicted returns for a single
    dataset and a given month. Command line arguments are toggled along 
    four dimensions:
    * The type of imputation. This corresponds to one of the suffixes of the
        `bcsignals_IMP.csv` files. In the paper, this would be one of 
        (`none`,`em`,`bllp`).
    * The type of algorithm: by default, this is standard PCA. If you want to
        use the scaled PCA of Huang et al. (2022), you would pass
        `--scaled_PCA=TRUE`
    * The year (1985-2021)
    * The month for the prediction (1-12)

    For replication of the paper, run the bash script `./3_run_pcr.sh`. 
    *Note that this submits many separate batch jobs, so make sure you 
    adapt to your computing environment.*
* The script `4a_one_forecast.R` runs one of the six forecasting algorithms
    described in section 5 of the paper for a given `bcsignals_IMP.csv`
    dataset. To run all the different forecasting algorithms for each dataset,
    run `./4_run_forecasts.sh`. It will also run forecasts for the different
    (micro,small,big) subsets described in the paper. Note that this will
    result in 6x6x3=108 different jobs. The script will create a 
    separate directory `YOUR OUTPUT PATH/forecasts/FORECAST_IMPUTATION_SUBSET`
    to store information on the run and the forecasted returns.
* The script `4b_fcast_table.R` combines the various forecasts-by-imputation
    created by running `./4_run_forecasts.sh`. It then creates the necessary data
    for Table 4 in the paper. 
* The script `5a_miss_plot.R` can be run simply with `Rscript 5a_miss_plot.R`
    once you have downloaded the data. It creates Figure 1 in the paper.
* The script `5b_corr_dist.R` can be run simply with `Rscript 5b_corr_dist.R`.
    once you have applied the Box-Cox transformations and EM imputations.
    It creates the plots for Figure 3 in the paper.
* The script `5c_betas_plots.R` can be run simply with `Rscript 5c_betas_plots.R`.
    once you have applied the Box-Cox transformations and EM imputations.
    It creates the plots for Figure IA.4 in the internet appendix.
* The script `5d_one_signal.R` can be run simply with `Rscript 5d_one_signal.R`.
    once you have applied the Box-Cox transformations and EM imputations.
    It creates the data for Table 3 in the paper.
* The script `5e_pcr_plots.R` can be run simply with `Rscript 5e_pcr_plots.R`
    once you have finished running the PCA regressions with `./3_run_pcr.sh`.
    It creates the plots for Figures 2, 5, 6, A.1, and A.2 in the paper. 

## Instructions to Replicators

- Run `install_packages.R` to install all dependencies.
- Adjust `FILEPATHS.R` to include a `data_path` and `out_path` for storage of
    data and results. _`data_path` must be able to hold 250+ GB of data and
    `out_path` must be able to hold 150+ GB of data._
- Adjust the grid submission commands in the `.sh` scripts to be applicable
    to your server.
- If you have a WRDS login, run `0a_download_data.R` to download the CRSP data
    from WRDS and signal data from Open Source Asset Pricing.
    - It is helpful to have a `~/.pgpass` file set up to log in to WRDS and
        run the script in batch mode.
        If not, there is an option to enter your WRDS username and password
        interactively.
    - If you don't have WRDS login, then you can skip the `0` section of the
        code and use the actual signals file and pseudo CRSP data provided.
- Once the data is downloaded, run `1_box_cox_transform.sh` to apply the
    Hawkins and Weisberg (2017) extension of the Box-Cox transformation to
    the signals.
- Run `2_est_imputations.sh` with the trailing command line arguments as 
    described above.
- Run `3_run_pcr.sh` to generate the PCR portfolios
- Run `4_run_forecasts.sh` to generate the different ML forecasts.
    - Run `Rscript 4b_forecast_table.R` to generate a CSV table of the different
      forecasts.
- Run each of the independent plots and table codes with `Rscript 5X_NAME.R`.
    Once you have run the above code, these can be run in any order (the
    alphabetical ordering is irrelevant).
- You can also run `Rscript 6_sim_correlations.R` to produce the plots for
    Figure IA.2 in the Internet Appendix.

### Details

- These programs were last run at various times in Fall and Winter 2023. 
- When running programs individually, note that ORDER IS IMPORTANT. 

## List of tables and programs

The provided code reproduces:

- [X] All tables and figures in the paper

Below, all tables are outputted to `$out_path/tables/` and all plots are
outputted to `$out_path/plots/`. `$out_path` is the output directory that you
specify in `FILEPATHS.R`.

Without downloading the actual CRSP data, none of the tables and figures
below can be fully reproduced, since the CRSP data is necessary to construct
several signals that are not made publically available by Chen and Zimmerman (2022).
Once you have downloaded the CRSP data from WRDS with `0a_download_data.R`,
however, you can reproduce all of the tables and figures below.

| Figure/Table #    | Program                  | Line Number | Output file                      | Note                            |
|-------------------|--------------------------|-------------|----------------------------------|---------------------------------|
| Table 1           | `0b_select_signals.R`    |175| `long_list.tex` | Makes the full .tex file, which is then curated to produce table.|
| Table 2           | `0b_select_signals.R`    |317 (Panel A); 394 (Panel B)| `pct_obs_one_signal_ptile.csv`; <br>`pct_obs_by_nsignal.csv` | Makes the full CSV files, which are then curated to produce table.|
| Table 3           | `5d_one_signal.R`    |245| `one_signal.csv` | Makes the full CSV file, which is then curated to produce table.|
| Table 4           | `4b_forecast_table`    |110; 111| `fcast_table_ew.csv`;<br> `fcast_table_vw.csv` | Makes the full CSV files, which are then curated to produce table.|
| Figure 1           | `5a_miss_plot.R`    |141 (function called for different months in 158)| `missplotcat_[MONTH].pdf` | Makes plots for June 1985,1990,...,2020. June 1990 is plotted in Figure 1. June 2000 and June 2010 are plotted in Figure IA.1 in internet appendix.|
| Figure 2           | `5e_pcr_plots.R`    |303; 309| `pca_separate_zoomout_expected_rets_main.pdf`;<br> `pca_separate_zoomout_sharpes_main.pdf` | |
| Figure 3           | `5b_corr_dist.R`    |184; 220| `cor_dist_[date].pdf`;<br> `pca_[date].pdf` | Loops through June 1990, 2000, and 2010 for `[date]`, making plots for each, which are taken together to make Figure 3.|
| Figure 4           | `5f_em_errors.R`    |192| `EM_error_rmse.pdf` | |
| Figure 5           | `5e_pcr_plots.R`    |303; 309| `pca_all_zoomout_expected_rets_main.pdf`;<br> `pca_all_zoomout_sharpes_main.pdf` | |
| Figure 6           | `5e_pcr_plots.R`    |303| `spca_separate_zoomout_expected_rets_main.pdf`;<br> `spca_separate_zoomin_expected_rets_main.pdf` | |

## References

Carhart, M. M., 1997. On Persistence in Mutual Fund Performance. Journal of Finance 52, 57-82.

Chen, A., Zimmerman, T., 2022. Open Source Cross-Sectional Asset Pricing. Critical Finance Review 11, 207-264.

Fama, E. F., French, K. R., 2015. A five-factor asset pricing model. Journal of Financial Economics 116, 1-22.

Hawkins, D. M., Weisberg, S., 2017. Combining the Box-Cox power and generalised log transformations to accommodate nonpositive responses in linear and mixed-effects linear models. South African Statistical Journal 51, 317-328. 

Huang, D., Jiang, F., Li, K., Tong, G., Zhou, G., 2022. Scaled PCA: A new approach to dimension reduction. Management Science 68, 1678-1695.


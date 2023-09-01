# options to use different kinds of imputation (or none)

#==============================================================================#
# Packages ----
#==============================================================================#

library(car)
library(doParallel)
library(data.table)
library(foreach)
library(zoo)

source("functions.R")

#==============================================================================#
# Option parsing ----
#==============================================================================#

option_list <- list(
    optparse::make_option(c("--impute_type"),
        type = "character", 
        default = "ind",
        help = "type of imputation: (ind, indsize, lastval)"),  
    optparse::make_option(c("--impute_vec"),
        type = "character",
        default = "../output/signals.txt",
        help = "a comma-separated list of values or .txt file to scan"),
    optparse::make_option(c("--data_path"),
        type = "character",
        default = "../output/bcsignals/",
        help = "where data is located"),
    optparse::make_option(c("--sample_start_year"),
        type = "numeric",
        default = 1985,
        help = "year that sample starts"),
    optparse::make_option(c("--sample_end_year"),
        type = "numeric",
        default = 2020,
        help = "year that sample ends"),
    optparse::make_option(c("--cores_frac"),
        type = "numeric",
        default = 1.0,
        help = "fraction of total cores to use")        
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Check that user selected a valid option
if (!(opt$impute_type %in% c("ind", "indsize", "lastval"))) {
    stop("Option `--impute_type` must be one of (ind, indsize, lastval)\n")
}

# Load file paths
getFilePaths()

# Get the anomalies as a nice vector
impute_vec <- unpackSignalList(FILEPATHS$signal_list)

#==============================================================================#
# Data pull ----
#==============================================================================#

# Read in the signals 
signals <- fread(paste0(FILEPATHS$data_path, "bcsignals/bcsignals_none.csv")) 
signals[, yyyymm := as.yearmon(yyyymm)]

# Read in the other data from CRSP 
crsp_data <- fread(FILEPATHS$data_path, "raw/crsp_data.csv")[, 
    .(permno, yyyymm, hsiccd, me)
]
crsp_data[, ":="(
    yyyymm = as.yearmon(yyyymm),
    sic3 = substr(ifelse(hsiccd < 1000, paste0("0", hsiccd), as.character(hsiccd)), 
        1, 3)
)]

signals <- merge(signals, crsp_data, by = c("permno", "yyyymm"))[, 
    .SD, 
    .SDcols = c("permno", "yyyymm", "sic3", "me", impute_vec)
]

#==============================================================================#
# Make new signals data ----
#==============================================================================#

# Sequence of yearmons to impute
yrmons <- seq(
    as.yearmon(paste0("Jan ", opt$sample_start_year)),
    as.yearmon(paste0("Dec ", opt$sample_end_year)),
    by = 1/12
)

# set up parallel comp
ncores = floor(parallel::detectCores()*opt$cores_frac)
doParallel::registerDoParallel(cores = ncores)

if (opt$impute_type %in% c("ind", "indsize")) {

    signals_new_list <- foreach::"%dopar%"(foreach::foreach(
      i = yrmons,  .packages = c('data.table','zoo')
      ), {
    
        i_file <- gsub("[[:space:]]", "", i)
        
        # Select data for given month
        good <- names(checkMinObs( # Get columns we can use
            signals[yyyymm == i, .SD, .SDcols = impute_vec],
            min_obs = 2
        ))
        
        signals_tmp <- signals[
            yyyymm == i,
            .SD,
            .SDcols = c("permno", "yyyymm", "sic3", "me", good)
        ]
    
        if (opt$impute_type == "ind") {
            # deep copy
            signals_tmp_dup <- copy(signals_tmp)
            
            # Get industry-level and overall cross-sectional means
            signals_tmp[,
                (good) := lapply(.SD, imputeVec),
                .SDcols = good,
                by = .(yyyymm, sic3)
            ]
            signals_tmp_dup[, (good) := lapply(.SD, imputeVec), .SDcols = good]
            
            # If there was no data for industry, replace with overall mean
            signals_tmp_mat <- as.matrix(signals_tmp[, ..good])
            signals_tmp_dup_mat <- as.matrix(signals_tmp_dup[, ..good])
            signals_tmp_mat[is.na(signals_tmp_mat)] <- signals_tmp_dup_mat[is.na(signals_tmp_mat)]
            
            # Final for output
            signals_tmp <- data.table(
                signals_tmp_mat, 
                permno = signals_tmp$permno,
                yyyymm = signals_tmp$yyyymm
            )
        } else if (opt$impute_type == "indsize") {
            # deep copy
            signals_tmp_dup <- copy(signals_tmp)
    
            # Get industry X size, industry-level, and overall cross-sectional means
            signals_tmp[,
                indsize_decile := tryCatch(pctileGroups(me, pctile = 0.1),
                    error = function(e) return(0)), # if can't make unique breaks, all together
                by = .(yyyymm, sic3)
            ][,
                (good) := lapply(.SD, imputeVec),
                .SDcols = good,
                by = .(yyyymm, sic3, indsize_decile)
            ]
            signals_tmp <- signals_tmp[, !c("indsize")] # remove so standard with dup
            signals_tmp_dup[, (good) := lapply(.SD, imputeVec), .SDcols = good]
            
            # If there was no data for industry, replace with overall mean
            signals_tmp_mat <- as.matrix(signals_tmp[, ..good])
            signals_tmp_dup_mat <- as.matrix(signals_tmp_dup[, ..good])
            signals_tmp_mat[is.na(signals_tmp_mat)] <- signals_tmp_dup_mat[is.na(signals_tmp_mat)]
            
            # Final for output
            signals_tmp <- data.table(
                signals_tmp_mat, 
                permno = signals_tmp$permno,
                yyyymm = signals_tmp$yyyymm
            )
        } # end impute if desired
        
        cat("signals_new for", as.character(i), "\n")
    
        # Return the imputed and un-transformed data
        return(signals_tmp)
    })

    # Combine 
    out <- rbindlist(signals_new_list, fill = T)

# Catch the last value imputations if specified (works differently)
} else if (opt$impute_type == "lastval") {
    out <- signals[
        order(permno, yyyymm), # order in advance
        .SD,
        .SDcols = c("permno", "yyyymm", impute_vec)
    ]
    out2 <- copy(out) # for cross-sectional means
    out[, # Use last value for imputation, up to 12 month lag
        (impute_vec) := lapply(.SD, function(x) imputeLastVal(x, yyyymm, k = 12)),
        .SDcols = impute_vec,
        by = .(permno)
    ]
    out2[, # Cross-sectional mean imputation
        (impute_vec) := lapply(.SD, function(x) imputeVec(x)),
        .SDcols = impute_vec,
        by = .(yyyymm)
    ]

    # Ensure same ordering
    out <- out[order(permno, yyyymm)]
    out2 <- out2[order(permno, yyyymm)]

    # If could not impute with last val, use XS mean. Indexing easier with matrix
    out_mat <- as.matrix(out[, .SD, .SDcols = impute_vec])
    out2_mat <- as.matrix(out2[, .SD, .SDcols = impute_vec])
    out_mat[is.na(out_mat)] <- out2_mat[is.na(out_mat)]
    
    # Final for output
    out <- data.table(
        out_mat, 
        permno = out$permno,
        yyyymm = out$yyyymm
    )
}

#==============================================================================#
# Output ----
#==============================================================================#

fwrite(out, 
    paste0(FILEPATHS$data_path, "bcsignals/bcsignals_", opt$impute_type, ".csv"))


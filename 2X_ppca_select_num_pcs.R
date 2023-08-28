#===============================================================================
# NOTE: (Jack 2023 08) The output of this doesn't directly feed into anything.
#   It's just used as a user reference for selecting the number of PCs to
#   use in the pPCA imputation. That's why it is "2X" (2 because it needs the
#   BC transformations to be done first)
#===============================================================================

library(data.table)
library(optparse)
library(pcaMethods)
library(zoo)

source('functions.R')

#===============================================================================
# Option parsing
#===============================================================================

# List of options passed as argments in shell
option_list <- list(
    optparse::make_option(c("--ym_start"),
        type = "character", 
        default = "Jun 1985",
        help = "start year-month of data to use in selection, formatted as 'Mon YYYY'"),
    optparse::make_option(c("--ym_end"),
        type = "character", 
        default = "Jun 1985",
        help = "end year-month of data to use in selection, formatted as 'Mon YYYY'"),
    optparse::make_option(c("--out_path"),
        type = "character", 
        default = "../output/ppca/",
        help = "directory to store output to"),
    optparse::make_option(c("--data_path"),
        type = "character", 
        default = "../data/",
        help = "directory to pull data from"),
    optparse::make_option(c("--impute_vec"),
        type = "character", default = "bm,mom6m",
        help = "a comma-separated list of values or .txt file to scan"),
    optparse::make_option(c("--maxiter"),
        type = "numeric", default = 10000,
        help = "a numeric value for the maximumum number of pPCA iterations"),
    optparse::make_option(c("--n_pcs"),
        type = "numeric", default = 2,
        help = "max number of principal components to estimate")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Get the anomalies as a nice vector
if (grepl("\\.txt", opt$impute_vec)) {
  opt$impute_vec <- scan(opt$impute_vec, character())
} else if (grepl(",", opt$impute_vec)) {
  opt$impute_vec <- trimws(do.call("c", strsplit(opt$impute_vec, ",")))
} else {
  stop("It seems that you did not pass a .txt file or comma-separated list",
    "to the `impute_vec` argument\n")
}

# Years and months to impute
yrmons <- seq(
    as.yearmon(opt$ym_start),
    as.yearmon(opt$ym_end),
    by = 1/12 # monthly interval for yearmon class
)

#===============================================================================
# Read in and set up data
#===============================================================================

# Signals data ====

# Read in the signals we want 
signals <- fread(paste0(opt$data_path, "bcsignals_none.csv")) 
setnames(signals, colnames(signals), tolower(colnames(signals)))

# Convert to proper yearmon format
signals[, yyyymm := as.yearmon(yyyymm)]

# Filter to desired data sample
signals <- signals[yyyymm %in% yrmons]

# A safety check to make sure we don't have any signals with a month of <2 obs
signals_good <- names(which(sapply( # `which` filters out the `FALSE` results
    signals[,
        lapply(.SD, function(x) { # has at least two unique observations
            sum(!is.na(x), na.rm = T) >= 2 & length(unique(x[!is.na(x)])) > 2
        }),
        .SDcols = opt$impute_vec, # the candidate list we passed in
        by = yyyymm # for each year-month
    ][, 
        !c("yyyymm") # don't want to check this column
    ], 
    all # and then check that it holds for all months
)))

# select to those variables
signals <- signals[, .SD, .SDcols = c("permno", "yyyymm", signals_good)]

# For checking in log
cat("There are a total of", length(signals_good), "out of", length(opt$impute_vec), 
    "signals with enough data.\n")

#===============================================================================
# Select the number of PCs
#===============================================================================

# Pre-scale data and select signals
train_data <- scale(signals[, .SD, .SDcols = signals_good])

# get the percentage of variance explained by number of PCs
doParallel::registerDoParallel(cores = detectCores())
pct_var <- rbindlist(foreach::"%dopar%"(foreach::foreach(n = 1:opt$n_pcs), {

    # pPCA with no return data ====

    pc_norets <- pca(train_data[, signals_good],
        method = 'ppca',
        nPcs = n,
        seed = 0, # Consistency
        maxIterations = opt$maxiter,
        scale = 'uv') # "unit variance" scaling (really a double check)

    return(data.table(n = n, pc_num = 1:n, R2cum = pc_norets@R2cum))

}))

# Print for reader
if(any(pct_var$R2cum >= 0.9)) {
    min_pcs <- min(pct_var[R2cum >= 0.9, n], na.rm = TRUE)
    cat("The minimum number of PCs to capture 90% of the variance is", min_pcs)
} else {
    max_var <- max(pct_var$R2cum, na.rm = TRUE)
    cat("There are not enough PCs to capture 90% of the variance in the signals.")
    cat("The max percentage captured is", round(max_var * 100, 2), "with ",
        opt$n_pcs)
}

# Output for reference
dir.create(opt$out_path, showWarnings=FALSE)
fwrite(
    pct_var, 
    paste0(
        opt$out_path,
        "ppca_r2cum_",
        gsub(" ", "", opt$ym_start),
        "_", 
        gsub(" ", "", opt$ym_end),
        ".csv"
    )
)


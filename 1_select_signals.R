
#===============================================================================
# Packages
#===============================================================================

library(data.table)
library(optparse)
library(zoo)

#===============================================================================
# Functions
#===============================================================================

source("functions.R")

#===============================================================================
# Hardcodes
#===============================================================================

option_list <- list(
    optparse::make_option(c("--n_signals", "-n"),
        type = "numeric", default = 10,
        help = "select first n signals in descending order of available observations"),
    optparse::make_option(c("--start_yr"),
        type = "numeric", default = 1980,
        help = "start year of data for selection"),
    optparse::make_option(c("--end_yr"),
        type = "numeric", default = 2018,
        help = "end year of data sample for selection"),
    optparse::make_option(c("--signals_file"),
        type = "character", default = "signals.txt",
        help = "name of files listing signals to use for imputations (output)")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

#===============================================================================
# Read in data
#===============================================================================

# Read in the signals we want 
signals <- fread("signed_predictors_dl_wide.csv") 
setnames(signals, colnames(signals), tolower(colnames(signals)))

# Ease of use
signals[, yyyymm := as.yearmon(as.Date(as.character(yyyymm*10 + 1), "%Y%m%d"))]

# Filter to sample being used
signals <- signals[opt$start_yr <= year(yyyymm) & year(yyyymm) <= opt$end_yr]

crsp_data <- fread("crsp_data.csv")[, .SD, .SDcols = !c("ret", "me")]
crsp_data[, yyyymm := as.yearmon(yyyymm)]

#===============================================================================
# Ensure we have all good variables
#===============================================================================

signals <- merge(signals, crsp_data, by = c("permno", "yyyymm"))

# Only want to use variables that we can work with in each month
bad_cols <- as.character(melt(signals[, # Counts for each month
        lapply(.SD, function(x) length(unique(x[!is.na(x)]))),
        .SDcols = !c("permno", "yyyymm"),
        by = yyyymm
    ][, # Count number of months with leq two unique observations 
        lapply(.SD, function(x) sum(x <= 2))
    ])[
        value >= 1,
        variable
    ])

# Matrix of counts of unique observations, whole sample
counts <- t(signals[, lapply(.SD, function(x) length(unique(x[!is.na(x)]))),
    .SDcols = !c("permno", "yyyymm", bad_cols)])

# Get first N in descending order of counts
signal_names <- rownames(counts)[order(-counts)][1:opt$n_signals]

# Output
writeLines(signal_names, opt$signals_file)



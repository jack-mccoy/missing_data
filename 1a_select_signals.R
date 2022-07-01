
#==============================================================================#
# Packages ----
#==============================================================================#

library(data.table)
library(optparse)
library(zoo)

#==============================================================================#
# Functions ----
#==============================================================================#

source("functions.R")

#==============================================================================#
# Hardcodes ----
#==============================================================================#

option_list <- list(
    optparse::make_option(c("--n_signals", "-n"),
        type = "numeric", default = 10,
        help = "select first n signals in descending order of available observations"),
    optparse::make_option(c("--sample_start_year"),
        type = "numeric", default = 1980,
        help = "start year of sample for analysis"),
    optparse::make_option(c("--sample_end_year"),
        type = "numeric", default = 2020,
        help = "end year of sample for analysis"),
    optparse::make_option(c("--signals_file"),
        type = "character", default = "signals.txt",
        help = "name of files listing signals to use for imputations (output)")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

#==============================================================================#
# Read in data ----
#==============================================================================#

# Read in the signals we want 
signals <- fread("../data/signed_predictors_dl_wide.csv") 
setnames(signals, colnames(signals), tolower(colnames(signals)))

# Ease of use
signals[, ":="(
    yrmon = as.yearmon(as.Date(as.character(yyyymm*10 + 1), "%Y%m%d")),
    yr = as.integer(substr(yyyymm, 1, 4)),
    mon = as.integer(substr(yyyymm, 5, 6))
)]

crsp_data <- fread("../data/crsp_data.csv")[, .SD, .SDcols = !c("ret", "me")]
crsp_data[, yrmon := as.yearmon(yyyymm)]

#==============================================================================#
# Ensure we have all good variables ----
#==============================================================================#

signals <- merge(signals, crsp_data[, !c("yyyymm")], by = c("permno", "yrmon"))

# Only want to use variables that we can work with in each month
# (remove binary signals in an overly-general way)
bad_cols <- as.character(melt(signals[ # Counts for each month
        opt$sample_start_year <= yr & yr <= opt$sample_end_year,
        lapply(.SD, function(x) length(unique(x[!is.na(x)]))),
        .SDcols = !c("permno", "yyyymm", "yrmon", "yr", "mon"),
        by = yyyymm
    ][, # Count number of months with leq two unique observations 
        lapply(.SD, function(x) sum(x <= 2))
    ])[
        value >= 1,
        variable
    ])

# Filter to month specified
# We do this AFTER making sure there are no bad months or it'll cause problems later on
signals_selection <- signals[
    opt$sample_start_year <= yr & yr <= opt$sample_end_year
]

# Matrix of counts of observations, selection month
counts <- t(signals_selection[, lapply(.SD, function(x) sum(!is.na(x))),
    .SDcols = !c("permno", "yyyymm", "yrmon", "yr", "mon", bad_cols)])

# Get first N in descending order of counts
signal_names <- rownames(counts)[order(-counts)][1:opt$n_signals]

# Output
writeLines(signal_names, opt$signals_file)



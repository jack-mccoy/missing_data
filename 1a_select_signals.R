
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
    optparse::make_option(c("--selection_year"),
        type = "numeric", default = 1980,
        help = "year of data for selection"),
    optparse::make_option(c("--selection_month"),
        type = "numeric", default = 6,
        help = "month of data for selection"),
    optparse::make_option(c("--signals_file"),
        type = "character", default = "signals.txt",
        help = "name of files listing signals to use for imputations (output)")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

if (0 >= opt$selection_month | opt$selection_month >= 13) {
    stop("`selection_month` option must be in [1,12]\n")
}

#==============================================================================#
# Read in data ----
#==============================================================================#

# Read in the signals we want 
signals <- fread("../data/signed_predictors_dl_wide.csv") 
setnames(signals, colnames(signals), tolower(colnames(signals)))

# Ease of use
signals[, yyyymm := as.yearmon(as.Date(as.character(yyyymm*10 + 1), "%Y%m%d"))]

crsp_data <- fread("../data/crsp_data.csv")[, .SD, .SDcols = !c("ret", "me")]
crsp_data[, yyyymm := as.yearmon(yyyymm)]

#==============================================================================#
# Ensure we have all good variables ----
#==============================================================================#

signals <- merge(signals, crsp_data, by = c("permno", "yyyymm"))

# Only want to use variables that we can work with in each month
# (remove binary signals in an overly-general way)
bad_cols <- as.character(melt(signals[ # Counts for each month
        opt$sample_start_year <= yyyymm & yyyymm <= opt$sample_end_year,
        lapply(.SD, function(x) length(unique(x[!is.na(x)]))),
        .SDcols = !c("permno", "yyyymm"),
        by = yyyymm
    ][, # Count number of months with leq two unique observations 
        lapply(.SD, function(x) sum(x <= 2))
    ])[
        value >= 1,
        variable
    ])

# Filter to month specified
# We do this AFTER making sure there are no bad months or it'll cause problems later on
signals_selection_month <- signals[
    opt$selection_year == year(yyyymm) &
    opt$selection_month == month(yyyymm)
]

# Matrix of counts of observations, selection month
counts <- t(signals_selection_month[, lapply(.SD, function(x) sum(!is.na(x))),
    .SDcols = !c("permno", "yyyymm", bad_cols)])

# Get first N in descending order of counts
signal_names <- rownames(counts)[order(-counts)][1:opt$n_signals]

# Output
writeLines(signal_names, opt$signals_file)



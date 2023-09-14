
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
    optparse::make_option(c("--start_yr"),
        type = "numeric", 
        default = 1985,
        help = "start year of data to impute"),
    optparse::make_option(c("--end_yr"),
        type = "numeric", 
        default = 2020,
        help = "start year of data to impute"),
    optparse::make_option(c("--maxiter"),
        type = "numeric", default = 10000,
        help = "a numeric value for the maximumum number of EM iterations"),
    optparse::make_option(c("--n_pcs"),
        type = "numeric", default = 2,
        help = "number of principal components to estimate"),
    optparse::make_option(c("--firmset"),
        type = "character", default = "big",
        help = paste0(
            "firm set to impute. one of (micro,small,big,all). ",
            "micro is below 20th ptile NYSE ME, small is 20th to 50th, ",
            "and big is 50th and above."
        ))
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Get file paths
getFilePaths()

# Get the anomalies as a nice vector
impute_vec <- unpackSignalList(FILEPATHS$signal_list)

# Years and months to impute
yrmons <- seq(
    as.yearmon(paste0('Jan ', opt$start_yr)),
    as.yearmon(paste0('Dec ', opt$end_yr)),
    by = 1/12 # monthly interval for yearmon class
)

#===============================================================================
# Read in and set up data
#===============================================================================

# Signals data ====

# Read in the signals we want 
signals <- fread(paste0(FILEPATHS$data_path, "bcsignals/bcsignals_none.csv")) 
setnames(signals, colnames(signals), tolower(colnames(signals)))

# Convert to proper yearmon format
signals[, yyyymm := as.yearmon(yyyymm)]

# Filter to desired data sample
signals <- signals[yyyymm %in% yrmons]

# Returns data ====

crsp_data <- fread(paste0(FILEPATHS$data_path, "raw/crsp_data.csv"))[,
    .(permno, yyyymm, ret, me)
]
crsp_data[, yyyymm := as.yearmon(yyyymm)]

# Combine and filter to firm set ====

signals <- merge(signals, crsp_data, by = c("permno", "yyyymm"))

if (opt$firmset %in% c("micro", "small", "big")) {
    signals <- filterFirms(signals, opt$firmset, FILEPATHS)
} else {
    if (opt$firmset != "all") {
        warning("`--firmset` was not correctly specified as one of (micro,small,big,all)!\n")
    }
    cat("Imputing for all firms...\n")
}

# Safety check here because usable signals will be firmset-specific ====

# A safety check to make sure we don't have any signals with a month of <2 obs
signals_good <- names(which(sapply( # `which` filters out the `FALSE` results
    signals[,
        lapply(.SD, function(x) { # has at least two unique observations
            sum(!is.na(x), na.rm = T) >= 2 & length(unique(x[!is.na(x)])) > 2
        }),
        .SDcols = impute_vec, # the candidate list we passed in
        by = yyyymm # for each year-month
    ][, 
        !c("yyyymm") # don't want to check this column
    ], 
    all # and then check that it holds for all months
)))

# select to those variables
signals <- signals[, .SD, .SDcols = c("permno", "yyyymm", signals_good)]

# For checking in log
cat("There are a total of", length(signals_good), "signals with enough data.\n")

#===============================================================================
# Estimate the model and impute for each month
#===============================================================================

# Register do parallel
doParallel::registerDoParallel(cores = parallel::detectCores())

imp_pc <- foreach::"%dopar%"(foreach::foreach(i = yrmons), {

    cat(paste0('Starting PPCA imputations for ', i, '\n'))
    
    # Subset to this month
    signals_i <- signals[yyyymm == i]
    
    # Estimate the PPCA model
    pc <- pca(as.matrix(signals_i[,.SD, .SDcols = signals_good]),
        method = 'ppca',
        nPcs = opt$n_pcs,
        seed = 0, # Consistency
        maxIterations = opt$maxiter,
        scale = 'uv') # "unit variance" scaling

    # Set up as matrix for convenient NA indexing. Make sure to also scale
    signals_i_mat <- as.matrix(scale(signals_i[, .SD, .SDcols = signals_good]))
    
    # Impute values where missing, otherwise original value
    na_idx <- is.na(signals_i_mat)
    signals_i_mat[na_idx] <- (pc@scores %*% t(pc@loadings))[na_idx]
    
    # Back to data.table
    signals_i_imp <- data.table(
        permno = signals_i$permno,
        yyyymm = signals_i$yyyymm,
        signals_i_mat
    )

    cat(paste0('Finished PPCA imputations for ', i, '\n'))

    # Output
    return(signals_i_imp)

})

#===============================================================================
# Combine the data and output
#===============================================================================

# Combine all the months into one file
out_data <- rbindlist(imp_pc, fill = TRUE)

# Output
fwrite(out_data, 
    paste0(FILEPATHS$data_path, "bcsignals/bcsignals_ppca", opt$n_pcs, "_", opt$firmset, ".csv"))


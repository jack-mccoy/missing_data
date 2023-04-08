
library(data.table)
library(optparse)
library(pcaMethods)
library(zoo)

source('functions.R')

#===============================================================================
# Option parsing
#===============================================================================

# Check if on cluster
on_cluster <- Sys.getenv('SGE_TASK_ID') != ""

# List of options passed as argments in shell
option_list <- list(
    optparse::make_option(c("--impute_yr"),
        type = "numeric", 
        default = ifelse(on_cluster, as.integer(Sys.getenv("SGE_TASK_ID")), 1990),
        help = "year of data to impute"),
    optparse::make_option(c("--out_path"),
        type = "character", 
        default = "../output/impute_ests/",
        help = "directory to store output to"),
    optparse::make_option(c("--impute_vec"),
        type = "character", default = "bm,mom6m",
        help = "a comma-separated list of values or .txt file to scan"),
    optparse::make_option(c("--maxiter"),
        type = "numeric", default = 10000,
        help = "a numeric value for the maximumum number of EM iterations"),
    optparse::make_option(c("--n_pcs"),
        type = "numeric", default = 2,
        help = "number of principal components to estimate")
    #optparse::make_option(c("--tol"),
    #    type = "numeric", default = 1e-6,
    #    help = "a numeric value for the convergence check")
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
yrmons <- zoo::as.yearmon(paste0(month.abb, " ", opt$impute_yr))

#===============================================================================
# Read in and set up data
#===============================================================================

# Read in the signals we want 
signals <- fread(paste0(opt$out_path, "bcsignals_none.csv")) 
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
cat("There are a total of", length(signals_good), "signals with enough data.\n")
print(signals_good)

print(head(signals))

#===============================================================================
# Estimate the model and impute for each month
#===============================================================================

# Register do parallel
doParallel::registerDoParallel(cores = parallel::detectCores())

imp_pc <- foreach::"%dopar%"(foreach::foreach(i = yrmons), {
    
    # Subset to this month
    signals_i <- signals[yyyymm == i]
    
    # Estimate the PPCA model
    pc <- pca(signals_i[,.SD, .SDcols = opt$impute_vec],
        method = 'ppca',
        nPcs = opt$n_pcs,
        seed = 0, # Consistency
        maxIterations = opt$maxiter,
        scale = 'uv') # "unit variance" scaling

    # Set up as matrix for convenient NA indexing. Make sure to also scale
    signals_i_mat <- as.matrix(scale(signals_i[, .SD, .SDcols = opt$impute_vec]))
    
    # Impute values where missing, otherwise original value
    na_idx <- is.na(signals_i_mat)
    signals_i_mat[na_idx] <- (pc@scores %*% t(pc@loadings))[na_idx]
    
    # Back to data.table
    signals_i_imp <- data.table(
        signals_i$permno,
        signals_i$yyyymm,
        signals_i_mat
    )

    # Output
    return(signals_i_imp)

})

#===============================================================================
# Combine the data and output
#===============================================================================

out_data <- rbindlist(imp_pc, fill = TRUE)

fwrite(out_data, paste0(opt$out_path, "bcsignals_ppca_", opt$impute_yr, ".csv"))


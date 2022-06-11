
#===============================================================================
# Packages
#===============================================================================

library(car) # Box-Cox stuff
library(data.table) # Standard data manipulation
library(optparse) # Python-like option parsing
library(RPostgres) # SQL query to WRDS
library(zoo) # yearmon convention is nice to work with here

#===============================================================================
# Option parsing
#===============================================================================

option_list <- list(
    optparse::make_option(c("--impute_yr"),
        type = "numeric", default = as.integer(Sys.getenv("SGE_TASK_ID")),
        help = "year of data to impute"),
    optparse::make_option(c("--boxcox"),
        type = "logical", default = FALSE, action = "store_true",
        help = "logical to indicate if Box-Cox transformations should be done before imputing"),
    optparse::make_option(c("--out_path"),
        type = "character", 
        default = "./",
        help = "directory to store output to"),
    optparse::make_option(c("--impute_vec"),
        type = "character", default = "bm,mom6m",
        help = "a comma-separated list of values or .txt file to scan"),
    optparse::make_option(c("--maxiter"),
        type = "numeric", default = 10000,
        help = "a numeric value for the maximumum number of EM iterations"),
    optparse::make_option(c("--tol"),
        type = "numeric", default = 1e-6,
        help = "a numeric value for the convergence check"),
    optparse::make_option(c("--force_convergence", "-f"),
        type = "logical", default = FALSE, action = "store_true",
        help = "logical to override maxiter and run the EM procedure until it converges to tol")
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

#===============================================================================
# Hardcodes
#===============================================================================

# Months to impute
yrmons <- zoo::as.yearmon(paste0(month.abb, " ", opt$impute_yr))

#===============================================================================
# Functions
#===============================================================================

source("functions.R")

#===============================================================================
# Read in data
#===============================================================================

# Read in the signals we want 
signals <- fread("../data/signed_predictors_dl_wide.csv") 
setnames(signals, colnames(signals), tolower(colnames(signals)))

# Ease of use
signals[, yyyymm := as.yearmon(as.Date(as.character(yyyymm*10 + 1), "%Y%m%d"))]

# Filter to test data
signals <- signals[yyyymm %in% yrmons]

crsp_data <- fread("../data/crsp_data.csv")
crsp_data[, yyyymm := as.yearmon(yyyymm)]

#===============================================================================
# Ensure we have all good variables
#===============================================================================

signals <- merge(signals, crsp_data, by = c("permno", "yyyymm"))

# Find vars even worth using (at least 2 different obs for variance) 
signals_good <- names(which(sapply(signals[, .SD, .SDcols = opt$impute_vec],
  function(x) sum(!is.na(x), na.rm = T) >= 2 & length(unique(x[!is.na(x)])) > 2)))

# select to those variables
signals <- signals[, .SD, .SDcols = c("permno", "yyyymm", signals_good)]

# Memory
rm(crsp_data)

#===============================================================================
# Box-Cox transformations and scaling
#===============================================================================

# Get vector of time periods to run through 
# (ensuring that we're only using yrmons in dataset)
yrmons <- unique(signals$yyyymm)
yrmons <- yrmons[order(yrmons)]

# Run the Box-Cox transformations ----

start_b <- Sys.time()

doParallel::registerDoParallel(cores = parallel::detectCores())
bctrans <- foreach::"%dopar%"(foreach::foreach(i = yrmons), {

    cat("Box-Cox transformations for", as.character(i), "\n")

    # Edit data for given month
    transformed <- signals[yyyymm == i][,
        (signals_good) := lapply(.SD, function(x) {
          tryCatch(
            {
                if (opt$boxcox) {
                    x <- winsorize(x, tail = 0.005) # 1% symmetric winsorization
                    params <- car::powerTransform(x, family = "bcnPower")
                    x <- scale(car::bcnPower(x, 
                        lambda = params$lambda, gamma = params$gamma))
                    attributes(x) <- c(attributes(x), 
                          list(
                              'bcn:lambda' = params$lambda,
                              'bcn:gamma' = params$gamma)
                          )
                } else {
                    x <- scale(winsorize(x, tail = 0.005))
                    attributes(x) <- c(attributes(x), 
                        list('bcn:lambda' = NA, 'bcn:gamma' = NA))
                }
                return(x)
            }, # if can't transform, just scale and winsor
            error = function(e) {
                warning("A column did not get BC-transformed for ",
                    as.character(i), " but BC transformation was specified")
                x <- scale(winsorize(x, tail = 0.005))
                attributes(x) <- c(attributes(x), 
                    list('bcn:lambda' = NA, 'bcn:gamma' = NA))
                return(x)
            }
          )
        }),
        .SDcols = signals_good
    ][, # Remove the permno and yyyymm columns
        .SD,
        .SDcols = signals_good
    ]

    # Getting the parameters to CSV to undo Box-Cox and scaling later ----

    params <- transformed[, lapply(.SD, function(x) {
        do.call("c", attributes(x))[c( # to ensure order consistency
            "scaled:center",
            "scaled:scale", 
            "bcn:lambda", 
            "bcn:gamma"
        )]
    })][, param := c( # Label the parameters in the data
        "scaled:center",
        "scaled:scale", 
        "bcn:lambda", 
        "bcn:gamma"
    )]

    # Each month gets its own CSV with 3 columns: param, variable, and value
    fwrite(melt(params, id.vars = "param"),
        paste0(opt$out_path, "bcn_scale_", gsub("[[:space:]]", "", 
            as.character(i)), ".csv"))

    return(transformed)

})

names(bctrans) <- as.character(yrmons)

# Timing for the log
end_b <- Sys.time()
bc_time <- difftime(end_b, start_b, unit = "mins")
cat("Box-Cox run time:", bc_time, "minutes\n")

# Memory
rm(signals)

#===============================================================================
# Imputations
#===============================================================================

# Timing for the log
start_i <- Sys.time()

# Getting started on the imputations in parallel ====

doParallel::registerDoParallel(cores = 28)
imp_par <- foreach::"%dopar%"(foreach::foreach(i = as.character(yrmons)), {

    cat("Starting imputations for", i, "\n")

    # Data prepping ----

    good <- names(checkMinObs(bctrans[[i]], min_obs = 2)) # Get cols we can use
    na_sort <- do.call("order", as.data.frame(-is.na(bctrans[[i]]))) #Sort by NA
    raw_i <- as.matrix(bctrans[[i]][na_sort, .SD, .SDcols = good]) # Final mat

    # Intial mean and covariance matrix
    E0 <- colMeans(raw_i, na.rm = T)
    R0 <- cov(raw_i, use = "pairwise.complete.obs")
    id <- diag(nrow = nrow(R0), ncol = ncol(R0))

    # Error checking. Catch missing means and covariance
    if (any(is.na(E0))) E0[which(is.na(E0))] <- 0
    if (any(is.na(R0))) R0[is.na(R0)] <- id[is.na(R0)]

    # Ensure covariance matrix is positive semi-definite ----

    eig <- eigen(R0)
    if (any(eig$values < 0)) {
        j <- 1
        while (any(eig$values < 0)) {
            R0 <- structure(
                (
                    eig$vectors %*% 
                    diag(ifelse(eig$values <= 0, 0, eig$values)) %*% 
                    t(eig$vectors)
                ),
                dimnames = list(rownames(R0), colnames(R0))
            )
            eig <- eigen(R0)
            if (j == 100) break
        }
        raw_i <- raw_i * matrix( # make it so data has same var as new cov mat
            rep(sqrt(diag(R0)), nrow(raw_i)), 
            byrow = T,
            ncol = ncol(raw_i)
        )
    } 

    # Imputation algorithm ----

    em_out <- mvn_emf(raw_i, E0, R0, maxiter = opt$maxiter, tol = opt$tol)

    # Ensure that imputation converged
    if (em_out$maxiter >= opt$maxiter & opt$force_convergence) {
      while (em_out$maxiter >= opt$maxiter) {
        cat("Imputations for", i, "did not converge. Trying again...\n")
        em_out <- mvn_emf(raw_i, 
            em_out$estE, em_out$estR, 
            maxiter = opt$maxiter, 
            tol = opt$tol)
      }
    }

    # Write the files ----

    # The estimated means
    write.csv(setNames(em_out$estE, good),
        paste0(opt$out_path, "estE_",
            gsub("[[:space:]]", "", as.character(i)), ".csv"))
    # The estimated covariances
    write.csv(structure(em_out$estR, dimnames = list(good, good)),
        paste0(opt$out_path, "estR_", 
            gsub("[[:space:]]", "", as.character(i)), ".csv"))

    cat("Finished imputations for", i, "\n")

    return(TRUE) # Confirmation of convergence, but not really necessary

})

# Timing thing to print to the output log for occassional check ====

end_i <- Sys.time()
imp_time <- difftime(end_i, start_i, units = "mins")
cat("Imputations for", opt$impute_yr, "ran in", imp_time, "minutes")


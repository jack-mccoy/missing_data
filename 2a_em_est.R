
#==============================================================================#
# Packages ----
#==============================================================================#

library(car) # Box-Cox stuff
library(data.table) # Standard data manipulation
library(optparse) # Python-like option parsing
library(RPostgres) # SQL query to WRDS
library(zoo) # yearmon convention is nice to work with here
library(tidyverse) # sorry Jack

#==============================================================================#
# Option parsing ----
#==============================================================================#

# check if on cluster
on_cluster = Sys.getenv('SGE_TASK_ID') != ""

option_list <- list(
    optparse::make_option(c("--impute_yr"),
        type = "numeric", 
        default = ifelse(on_cluster, as.integer(Sys.getenv("SGE_TASK_ID")), 1990),
        help = "year of data to impute"),
    optparse::make_option(c("--boxcox"),
        type = "logical", default = FALSE, action = "store_true",
        help = "logical to indicate if Box-Cox transformations should be done before imputing"),
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
    optparse::make_option(c("--tol"),
        type = "numeric", default = 1e-6,
        help = "a numeric value for the convergence check"),
    optparse::make_option(c("--force_convergence", "-f"),
        type = "logical", default = FALSE, action = "store_true",
        help = "logical to override maxiter and run the EM procedure until it converges to tol"),
    optparse::make_option(c("--cores_frac"),
        type = "numeric", default = 1.0,
        help = "fraction of total cores to use")    
)

  
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

    # debug
    opt$impute_vec = '../output/signals_20.txt'

# Get the anomalies as a nice vector
if (grepl("\\.txt", opt$impute_vec)) {
  opt$impute_vec <- scan(opt$impute_vec, character())
} else if (grepl(",", opt$impute_vec)) {
  opt$impute_vec <- trimws(do.call("c", strsplit(opt$impute_vec, ",")))
} else {
  stop("It seems that you did not pass a .txt file or comma-separated list",
    "to the `impute_vec` argument\n")
}

#==============================================================================#
# Setup ----
#==============================================================================#

# Months to impute
yrmons <- zoo::as.yearmon(paste0(month.abb, " ", opt$impute_yr))

source("functions.R")

# make output folder
dir.create(opt$out_path, showWarnings = F)

closeAllConnections()

#==============================================================================#
# Read in data ----
#==============================================================================#

# Read in the bc-transformed signals 
bcsignals = fread('../output/bcsignals/bcsignals_none.csv') 
bcsignals[ , yyyymm := as.yearmon(yyyymm)]  # careful with reading yearmon format from csv!

# Filter to test data
bcsignals <- bcsignals[yyyymm %in% yrmons, .SD, .SDcols = c("permno", "yyyymm", opt$impute_vec)]


#==============================================================================#
# AR1 Residuals ----
#==============================================================================#

ar1_grouping = c('signalname') # grouping for ar1 model

# reshape to long
bclong = bcsignals %>% 
  pivot_longer(
    cols = !c('permno','yyyymm'), names_to = 'signalname'
  ) %>% 
  arrange(signalname, permno, yyyymm, ) %>% 
  setDT()

bclong[
  order(signalname, permno , yyyymm)
  , lag_value := shift(value, type = 'lag')
  , by = c('signalname', 'permno')
]

# estimate AR1
ar1_param = bclong[ 
  , list(
    ar1_slope = coef(lm(value ~ 0 + lag_value))[1]
  )
  , by = ar1_grouping
]

# find residuals
bclong = bclong %>% left_join(ar1_param, by = ar1_grouping) %>% 
  mutate(
    ar1_pred = lag_value*ar1_slope
    , ar1_resid = value - ar1_pred
  )


# reshape back to wide
bcsignals = bclong %>% select(permno,yyyymm,signalname,ar1_resid) %>% 
  pivot_wider(names_from = signalname, values_from = ar1_resid) %>% 
  setDT()


#==============================================================================#
# Imputation parameter estimates ----
#==============================================================================#
  
  # DEBUG
  yrmons = max(bcsignals$yyyymm)
  i = yrmons

# Timing for the log
start_i <- Sys.time()

# imp_par <- foreach::"%dopar%"(foreach::foreach(i = as.character(yrmons)), {
  imp_par <- foreach::"%do%"(foreach::foreach(i = as.character(yrmons)), {  

    cat("Starting imputations for", i, "\n")

    ## Data prepping ----
    bcsmall = bcsignals[ yyyymm == i , .SD, .SDcols = opt$impute_vec ]

    good <- names(checkMinObs(bcsmall, min_obs = 2)) # Get cols we can use
    na_sort <- do.call("order", as.data.frame(-is.na(bcsmall))) #Sort by NA
    raw_i <- as.matrix(bcsmall[na_sort, .SD, .SDcols = good]) # Final mat

    # Initialize mean and cov matrix
    #   enforce 0s on diagonal, as in the norm2 package
    E0 = colMeans(raw_i, na.rm = T)
    R0 = diag(diag(cov(raw_i, use = "pairwise.complete.obs"))) 

    # Error checking. Catch missing means and covariance
    if (any(is.na(E0))) stop(paste0('A signal has zero obs this month ', i))
    if (any(is.na(R0))) stop(paste0('A signal has < 2 obs this month ', i))

    # Imputation algorithm ----
    em_out <- mvn_emf(raw_i, E0, R0, maxiter = opt$maxiter, tol = opt$tol,
        update_estE = FALSE)

    # Ensure that imputation converged
    if (em_out$maxiter >= opt$maxiter & opt$force_convergence) {
      while (em_out$maxiter >= opt$maxiter) {
        cat("Imputations for", i, "did not converge. Trying again...\n")
        em_out <- mvn_emf(raw_i, 
            em_out$estE, em_out$estR, 
            maxiter = opt$maxiter, 
            tol = opt$tol,
            update_estE = FALSE)
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


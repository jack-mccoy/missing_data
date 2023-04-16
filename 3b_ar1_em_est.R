
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



option_list <- list(
  optparse::make_option(c("--impute_yr"),
                        type = "numeric", 
                        default = 2015,
                        help = "year of data to impute"),
  optparse::make_option(c("--impute_vec"),
                        type = "character", default = "bm,mom6m",
                        help = "a comma-separated list of values or .txt file to scan"),
  optparse::make_option(c("--ar1_sample_length"),
                        type = "numeric", 
                        default = 60,
                        help = "Months to use in AR1 estimate"),  
  optparse::make_option(c("--maxiter"),
                        type = "numeric", default = 2000,
                        help = "a numeric value for the maximumum number of EM iterations"),
  optparse::make_option(c("--tol"),
                        type = "numeric", default = 1e-5,
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
opt$impute_vec = '../output/signals_75.txt'
opt$cores_frac = 0.5

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

source("functions.R")

# Months to impute
yrmon_list <- zoo::as.yearmon(paste0(month.abb, " ", opt$impute_yr))

closeAllConnections()

# Read in the bc-transformed signals 
bcsignals = fread('../output/bcsignals/bcsignals_none.csv') 
bcsignals[ , yyyymm := as.yearmon(yyyymm)]  # careful with reading yearmon format from csv!

# keep only selected signals and months needed for ar1 samples
bcsignals <- bcsignals[
  yyyymm >= as.yearmon(paste0(opt$impute_yr, '-01')) - opt$ar1_sample_length + 1 
  & yyyymm <= as.yearmon(paste0(opt$impute_yr, '-12'))
  , .SD, .SDcols = c("permno", "yyyymm", opt$impute_vec)
]

# reshape to long
bclong = bcsignals %>% 
  pivot_longer(
    cols = !c('permno','yyyymm'), names_to = 'signalname'
  ) %>% 
  arrange(signalname, permno, yyyymm, ) %>% 
  setDT()

# find signal update frequency
signaldoc = fread('../data/SignalDoc.csv') %>% 
  filter(Cat.Signal == 'Predictor') %>% 
  transmute(signalname = tolower(Acronym)
            , port_start = `Start Month`, port_period = `Portfolio Period`) %>% 
  mutate(
    port_start = if_else(is.na(port_start), 6, port_start)
    , port_period = if_else(is.na(port_period), 1, port_period)
    , port_period = pmin(port_period, 12)
  ) %>% 
  select(-port_start) # not used right now, let's be stingy

bclong = bclong %>% left_join(signaldoc, by = c('signalname'))

# memory
rm(bcsignals)

#==============================================================================#
# AR1 Residuals Function ----
#==============================================================================#


ar1_grouping = c('signalname') # grouping for ar1 model

# prob don't want to parallelize this since it takes so much ram
# bclong is roughly 5-10 gigs with 100 signals, so you really can't parallelize much anyway

ar1long <- foreach::"%do%"(foreach::foreach(
  i = yrmon_list,  .packages = c('data.table','zoo','tidyverse')
  , .combine = rbind
), {
  
  print(paste0('running ar1 model for ', i))
  
  # don't look ahead
  bclong2 = bclong[yyyymm >= i - opt$ar1_sample_length + 1  & yyyymm <= i]
  
  # deal with signal update frequency
  bclong2 = bclong2 %>% mutate(
    mon = round(12*(yyyymm - floor(yyyymm))+1)
  ) %>% 
    filter(
      port_period == 1 | (
        port_period == 3 & mon %in%  (month(i) + 3*(-12:12))
      ) | (
        port_period == 6 & mon %in%  (month(i) + 6*(-12:12))
      ) | (
        port_period == 12 & mon %in%  (month(i) + 12*(-12:12))
      )
    ) %>% 
    select(permno,yyyymm,signalname,value)
  
  # lag
  bclong2[
    order(signalname, permno , yyyymm)
    , lag_value := shift(value, type = 'lag')
    , by = c('signalname', 'permno')
  ]
  
  # estimate AR1
  ar1_param = bclong2[ 
    , list(
      ar1_slope = coef(lm(value ~ 0 + lag_value))[1]
    )
    , by = ar1_grouping
  ]
  
  # make xs for current month i
  xs = bclong2[yyyymm == i] %>% left_join(ar1_param, by = ar1_grouping) %>% 
    mutate(
      pred = lag_value*ar1_slope
      , resid = value - pred
    ) %>% 
    # if there's no data, the model says E(signal) = 0
    mutate(
      pred = if_else(is.na(pred), 0, pred)
    ) %>% 
    select(permno, yyyymm, signalname, value, pred, resid)
  
  return(xs)

})

# memory: bclong is ar1_sample_length/12 times as big as ar1long
rm(bclong)
rm(bclong2)


#==============================================================================#
# EM on Residuals ----
#==============================================================================#


# Timing for the log
start_i <- Sys.time()


ncores = floor(parallel::detectCores()*opt$cores_frac)
doParallel::registerDoParallel(cores = ncores)


bcsignals_emar1 <- foreach::"%dopar%"(foreach::foreach(
  i = as.character(yrmon_list),  .packages = c('data.table','zoo','tidyverse')
  , .combine = rbind
), {
  
  cat("Starting imputations for", as.character(i), "\n")
  
  # turn ar1long into cross-sections
  xs = list()
  
  xs$resid = ar1long[yyyymm == i] %>% select(permno,yyyymm,signalname,resid) %>% 
    pivot_wider(names_from = signalname, values_from = resid) %>% 
    arrange(permno) %>% 
    setDT()
  
  xs$pred = ar1long[yyyymm == i] %>% select(permno,yyyymm,signalname,pred) %>% 
    pivot_wider(names_from = signalname, values_from = pred) %>% 
    arrange(permno) %>%   
    setDT()
  
  xs$obs = ar1long[yyyymm == i] %>% select(permno,yyyymm,signalname,value) %>% 
    pivot_wider(names_from = signalname, values_from = value) %>% 
    arrange(permno) %>%   
    setDT()

  # Data prep --
  
  # sort by missingness, arrange columns
  na_sort <- do.call("order", as.data.frame(-is.na(xs$resid))) 
  xs$resid = xs$resid[ na_sort,   ] %>% select(c(permno,yyyymm, opt$impute_vec))
  xs$pred = xs$pred[na_sort, ] %>% select(c(permno,yyyymm, opt$impute_vec))
  xs$obs = xs$obs[na_sort, ] %>% select(c(permno,yyyymm, opt$impute_vec))  
  
  # make matrices (note they're sorted above)
  mat_resid <- xs$resid %>% select(-c(permno,yyyymm)) %>% as.matrix
  mat_pred  <- xs$pred %>% select(-c(permno,yyyymm)) %>% as.matrix
  mat_obs  <- xs$obs %>% select(-c(permno,yyyymm)) %>% as.matrix
  
  # Initialize mean and cov matrix
  #   force means to be zero, since we have demeaned
  #   enforce 0s on diagonal, as in the norm2 package
  E0 = 0*colMeans(mat_resid, na.rm = T)
  R0 = diag(diag(cov(mat_resid, use = "pairwise.complete.obs"))) 
  
  # Error checking. Catch missing means and covariance
  if (any(is.na(E0))) stop(paste0('A signal has zero obs this month ', i))
  if (any(is.na(R0))) stop(paste0('A signal has < 2 obs this month ', i))
  
  # EM impute residuals ---
  em_out <- mvn_emf(mat_resid, E0, R0, maxiter = opt$maxiter, tol = opt$tol,
                    update_estE = FALSE)
  
  # output log file if no convergence
  if (em_out$maxiter >= opt$maxiter){
    sink('3b_ar1_em_est.err')
    
    print('ar1_em_est.R error: divergence')
    i
    opt$impute_vec
    print('iter, tolerance')
    em_out$maxiter
    em_out$tol
    print('E0')
    E0
    
    sink()
  }
  
  # Ensure that imputation converged
  if (em_out$maxiter >= opt$maxiter & opt$force_convergence) {
    while (em_out$maxiter >= opt$maxiter) {
      cat("Imputations for", i, "did not converge. Trying again...\n")
      em_out <- mvn_emf(mat_resid, 
                        em_out$estE, em_out$estR, 
                        maxiter = opt$maxiter, 
                        tol = opt$tol,
                        update_estE = FALSE)
    }
  }
  
  # Combine EM residuals with predictions
  mat_imp = mat_pred + em_out$Ey
  
  # replace with observed if observed
  mat_obs_zeros = mat_obs
  mat_obs_zeros[which(is.na(mat_obs_zeros))] = 0
  mat_imp = as.numeric(is.na(mat_obs)) * mat_imp + as.numeric(!is.na(mat_obs)) * mat_obs_zeros
  
  
  # add labels, make back into data table
  imp = data.table(
    permno = xs$resid$permno,
    yyyymm = xs$resid$yyyymm,
    mat_imp
  )
  

  cat("Finished imputations for",  as.character(i), "\n")
  
  return(imp)
  
})

# Timing thing to print to the output log for occasional check ===

end_i <- Sys.time()
imp_time <- difftime(end_i, start_i, units = "mins")
cat("Imputations for", opt$impute_yr, "ran in", imp_time, "minutes")


#==============================================================================#
# Save ----
#==============================================================================#

dir.create('../output/emar1_intermediate', showWarnings = F)
fwrite(bcsignals_emar1, paste0('../output/emar1_intermediate/bcsignals_emar1_',opt$impute_yr, '.csv' ))

#==============================================================================#
# Sum stats for the console ----
#==============================================================================#


bcsignals_emar1 %>% 
  pivot_longer(cols = !c(permno,yyyymm)) %>% 
  group_by(name, yyyymm) %>% 
  summarize(
    nobs = sum(!is.na(value))
    , sd = sd(value)
    , m = mean(value)
  ) %>% 
  group_by(name) %>% 
  summarize(
    across(c(nobs,sd,m), ~ mean(.x))
  ) %>% 
  print(n = 200)

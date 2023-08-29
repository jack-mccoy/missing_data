# input: bcsignals_none.csv
# output: many bcsignals_emar1_YYYY.csv files (temporary)

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
    optparse::make_option(c("--em_type"),
        type = "character", default = "regular",
        help = "one of (regular, ar1)"),
    optparse::make_option(c("--output_ts_pred"),
        type = "logical", default = FALSE,
        help = "output the prediction of time-series model?"),  
    optparse::make_option(c("--out_path"), # Need this flexibility. Files are too big to store in non-scratch
        type = "character", default = "../output/",
        help = "directory including input and output files"),
    optparse::make_option(c("--impute_yr"),
       type = "numeric", 
       default = ifelse(on_cluster, as.integer(Sys.getenv("SGE_TASK_ID")), 2015),
       help = "year of data to impute"),
    optparse::make_option(c("--impute_months"),
       type = "character", default = '1,2,3,4,5,6,7,8,9,10,11,12',
       help = "months of data to impute"),    
    optparse::make_option(c("--impute_vec"),
       type = "character", default = "bm,mom6m",
       help = "a comma-separated list of values or .txt file to scan"),
    optparse::make_option(c("--ar1_sample_length"),
       type = "numeric", 
       default = 5,
       help = "Years to use in AR1 estimate"),  
    optparse::make_option(c("--winsor_p"),
       type = "numeric", default = 0,
       help = "Winsorize the extreme p fraction of residuals (by month-signalname)"),      
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
       type = "numeric", default = 0.5,
       help = "fraction of total cores to use")    
)


opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Check that user passed valid option
if (!(opt$em_type %in% c("regular", "ar1"))) {
    stop("Option `--em_type` must be one of (regular, ar1)\n")
}

# for consistency
if (opt$em_type == 'regular'){
  opt$ar1_sample_length = 1
}


# Get the anomalies as a nice vector
if (grepl("\\.txt", opt$impute_vec)) {
  opt$impute_vec <- scan(opt$impute_vec, character())
} else if (grepl(",", opt$impute_vec)) {
  opt$impute_vec <- trimws(do.call("c", strsplit(opt$impute_vec, ",")))
} else {
  stop("It seems that you did not pass a .txt file or comma-separated list",
       "to the `impute_vec` argument\n")
}

# convert impute_months to numeric
opt$impute_months = as.numeric(strsplit(opt$impute_months, ',')[[1]])

#==============================================================================#
# Setup ----
#==============================================================================#

source("functions.R")

closeAllConnections()

# Read in the bc-transformed signals 
bcsignals = fread(paste0(opt$out_path, "bcsignals/bcsignals_none.csv")) 
bcsignals[ , yyyymm := as.yearmon(yyyymm)]  # careful with reading yearmon format from csv!

# keep only selected signals and months needed for ar1 samples
bcsignals <- bcsignals[
  yyyymm >= as.yearmon(paste0(opt$impute_yr, '-01')) - opt$ar1_sample_length + 1/12 
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
    port_start = ifelse(is.na(port_start), 6, port_start)
    , port_period = ifelse(is.na(port_period), 1, port_period)
    , port_period = pmin(port_period, 12)
  ) %>% 
  select(-port_start) # not used right now, let's be stingy

bclong = bclong %>% left_join(signaldoc, by = c('signalname'))

# memory
rm(bcsignals)

#==============================================================================#
# AR1 Residuals, or Not ----
#==============================================================================#
# makes bclong2

# Months to impute
yrmon_list <- as.yearmon(paste(opt$impute_yr, opt$impute_months, sep = '-'))

if (opt$em_type == 'regular'){
  # if EM type is regular, it's just an AR1 model with zero persistence.
  
  bclong2 = bclong %>% select(-port_period) %>% 
    mutate(pred = 0, resid = value) 
  
} else if (opt$em_type == 'ar1'){
  # if EM type is ar1, run ar1 model to decompose value into pred and resid
  
  # prob don't want to parallelize this since it takes so much ram
  # bclong is roughly 5-10 gigs with 100 signals, so you really can't parallelize much anyway
  bclong2 <- foreach::"%do%"(foreach::foreach(
    i = yrmon_list,  .packages = c('data.table','zoo','tidyverse')
    , .combine = rbind
  ), {
    
    ar1_grouping = c('signalname') # grouping for ar1 model, hardcode for now
    
    print(paste0('running ar1 model for ', i))
    
    # don't look ahead
    bctemp = bclong[yyyymm >= i - opt$ar1_sample_length + 1/12  & yyyymm <= i]
    
    # deal with signal update frequency
    bctemp = bctemp %>% mutate(
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
    bctemp[
        order(signalname, permno , yyyymm),
        lag_value := shift(value, type = 'lag'),
        by = c('signalname', 'permno')
    ]
    
    # estimate AR1
    ar1_param = bctemp[ 
      , list(
        ar1_slope = tryCatch(coef(lm(value ~ 0 + lag_value))[1], error = function(e) 0)
      )
      , by = ar1_grouping
    ]
    
    # make xs for current month i
    xs = bctemp[yyyymm == i] %>% left_join(ar1_param, by = ar1_grouping) %>% 
      mutate(
        pred = lag_value*ar1_slope
        , resid = value - pred
      ) %>% 
      # if there's no data, the model says E(signal) = 0
      mutate(
        pred = ifelse(is.na(pred), 0, pred)
      ) %>% 
      select(permno, yyyymm, signalname, value, pred, resid)
    
    return(xs)
    
  })
  
  rm(bctemp)
  
} # end if opt$em_type

# memory: bclong can be really big if ar1_sample_length is long
rm(bclong)


#==============================================================================#
# Winsorize Residuals ----
#==============================================================================#
# ar1 model can make kurtosis even more extreme

if (opt$winsor_p > 0){
  bclong2[
    , by = c('yyyymm','signalname')
    , resid := winsorize(resid, tail = opt$winsor_p/2)
  ][
    , value := pred + resid
  ]
  
} # end if opt$winsor_p > 0


#==============================================================================#
# EM on Residuals ----
#==============================================================================#


# Timing for the log
start_i <- Sys.time()

if (opt$cores_frac > 0){
  ncores = floor(parallel::detectCores()*opt$cores_frac)
  doParallel::registerDoParallel(cores = ncores)
}

bcsignals_emar1 <- foreach::"%dopar%"(foreach::foreach(
  i = as.character(yrmon_list),  .packages = c('data.table','zoo','tidyverse')
  , .combine = rbind
), {
  
  cat("Starting imputations for", as.character(i), "\n")
  
  # turn bclong2 into cross-sections
  xs = list()
  
  xs$resid = bclong2[yyyymm == i] %>% select(permno,yyyymm,signalname,resid) %>% 
    pivot_wider(names_from = signalname, values_from = resid) %>% 
    arrange(permno) %>% 
    setDT()
  
  xs$pred = bclong2[yyyymm == i] %>% select(permno,yyyymm,signalname,pred) %>% 
    pivot_wider(names_from = signalname, values_from = pred) %>% 
    arrange(permno) %>%   
    setDT()
  
  xs$obs = bclong2[yyyymm == i] %>% select(permno,yyyymm,signalname,value) %>% 
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
    Sys.time()
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

dir.create(paste0(opt$out_path, 'em_intermediate'), showWarnings = F)
if (opt$em_type == "ar1") {
    filename <- paste0(opt$out_path, 'em_intermediate/bcsignals_emar1_',opt$impute_yr, '.csv' )
} else {
    filename <- paste0(opt$out_path, 'em_intermediate/bcsignals_em_',opt$impute_yr, '.csv' )
}
fwrite(bcsignals_emar1, filename)
sink(paste0(opt$out_path, "em_intermediate/readme.log"))
Sys.time()
opt
sink()

if (opt$output_ts_pred){
    ts_pred = bclong2 %>% select(permno,yyyymm,signalname,pred) %>% 
        pivot_wider(names_from = signalname, values_from = pred) %>% 
        arrange(permno,yyyymm) 
    
    fwrite(ts_pred,
        paste0(opt$out_path, "em_intermediate/ts_prediction_",opt$impute_yr, ".csv"))
}

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

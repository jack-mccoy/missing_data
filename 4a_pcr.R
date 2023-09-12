# runs pcr for (year, month) =  (opt$iter_year, opt$iter_month)

#==============================================================================#
# Packages ----
#==============================================================================#

library(doParallel)
library(data.table)
library(foreach)
library(stringr)
library(zoo)

closeAllConnections()

source("functions.R")

#==============================================================================#
# Option parsing ----
#==============================================================================#

# check if on cluster
on_cluster = Sys.getenv('SGE_TASK_ID') != ""

option_list <- list(
    #optparse::make_option(c("--out_path"),
    #    type = "character", default = "../output/pca_returns/em/",
    #    help = "directory to store output to"),
    #optparse::make_option(c("--signals_keep"),
    #    type = "character", default = "../output/signals_10.txt",
    #    help = "a comma-separated list of values or .txt file to scan"),
    optparse::make_option(c("--signal_file"),
        type = "character", default = "bcsignals_em.csv",
        help = "path of imputed file"),
    optparse::make_option(c("--iter_year"),
        type = "numeric", default = 1990,
        help = "year for making prediction (numeric)"),
    optparse::make_option(c("--iter_month"),
        type = "numeric", 
        default = ifelse(on_cluster, as.integer(Sys.getenv("SGE_TASK_ID")), 5),
        help = "month for making prediction (numeric)"),
    optparse::make_option(c("--n_yrs"),
        type = "numeric", default = 7,
        help = "number of years to run for Fama-Macbeth principal component regressions"),
    optparse::make_option(c("--scaled_pca"),
        type = "logical", default = FALSE,
        help = "logical to indicate if Huang et al 2022 scaled pca should be used"),  
    optparse::make_option(c("--scaled_pca_weight"),
        type = "character", default = "ew",
        help = "ew or vw: allows value-weighted scaled pca"),
    optparse::make_option(c("--quantile_prob"),
        type = "numeric", default = 0.2,
        help = paste0("the ratio out of 1 used to form long-short portfolios ",
            "(i.e., 0.2 means quintile portfolios)")),
    optparse::make_option(c("--n_pcs"),
        type = "numeric", default = 25,
        help = "maximum number of principal components to use in PCRs"),
    optparse::make_option(c("--cores_frac"),
        type = "numeric", default = 1.0,
        help = "fraction of total cores to use")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Get the file paths
getFilePaths()

# Signals to use
signals_keep <- unpackSignalList(FILEPATHS$signal_list)

#==============================================================================#
# Get prediction month from array iteration map ----
#==============================================================================#

pred_mon <- as.yearmon(opt$iter_year + opt$iter_month/12 - 1/12)

# Log file
cat("PC regressions for ", as.character(pred_mon), "\n")

#==============================================================================#
# Data pull ----
#==============================================================================#

# Read in data and merge
signals <- fread(paste0(FILEPATHS$data_path, "bcsignals/", opt$signal_file))
setnames(signals, colnames(signals), tolower(colnames(signals))) # ensure lower

crsp_data <- fread(paste0(FILEPATHS$data_path, "raw/crsp_data.csv"))[, 
    .(permno, yyyymm, ret, me)
]
signals <- merge(signals, crsp_data, by = c("permno", "yyyymm"))

# Memory
rm(crsp_data)

# Nice yearmons format
signals[, yyyymm := as.yearmon(yyyymm)]

# Filter down to months we are actually using for PCA
signals <- signals[
    (pred_mon - opt$n_yrs) <= yyyymm &
    yyyymm <= pred_mon
]

# Simple mean imputations. Shouldn't matter for EM-imputed data
# Easiest to just catch it all here for non-imputed data
signals[,
    (signals_keep) := lapply(.SD, imputeVec, na.rm = T),
    .SDcols = signals_keep,
    by = .(yyyymm)
]

# Complete cases are necessary for principal components. Ensuring here
signals <- signals[complete.cases(
  signals[, .SD, .SDcols = c("permno", "yyyymm", "ret", "me", signals_keep)]
)]

#==============================================================================# 
# Timing adjustments ----
#==============================================================================#

# Align signals with next month's returns
signals[
    order(permno, yyyymm),
    bh1m := shift(ret, 1, type = "lead"), # we want to predict next month's ret
    by = permno
]

# Define 'time_avail_m': the month of return prediction where signal was available
signals[, time_avail_m := yyyymm + 1/12]

#==============================================================================# 
# Huang et al. Scaling (optional) ----
#==============================================================================#

if (opt$scaled_pca){
  
    # set up weights
    if (opt$scaled_pca_weight == 'ew'){
        signals$tempw = 1
    } else if (opt$scaled_pca_weight == 'vw') {
        signals$tempw = signals$me
    } else {
        stop('opt$scaled_pca_weight invalid')
    }
    
    # regress bh1m on each signal, dropping current month's bh1m
    temp = lapply(
        signals_keep,
        function(signalname){
            summary(lm(
                paste0('bh1m ~ ', signalname), signals[time_avail_m < pred_mon],
                weights = tempw
            ))$coefficients[signalname, 'Estimate']
        }
    )
    slopedat <- data.table(
        signalname = signals_keep,
        slope = as.numeric(temp)
    )
    
    # re-scale with slope
    slopemat <- matrix(1, dim(signals)[1], 1) %*% t(as.matrix(slopedat$slope))
    
    signals <- cbind(
        signals[ , .SD, .SDcols = signals_keep] * slopemat,
        signals[ , .SD, .SDcols = !signals_keep] 
    )
    pc <- prcomp(signals[, .SD, .SDcols = signals_keep], 
        center = FALSE, scale = FALSE)
} else {
    pc <- prcomp(signals[, .SD, .SDcols = signals_keep], 
        center = TRUE, scale = TRUE)
}

#==============================================================================# 
# Run regressions ----
#==============================================================================#

# Get principal components and form regression data
reg_data <- data.table(signals[, .(permno, time_avail_m, bh1m, me)], pc$x)

# Max number of PCs we want to run
n_pcs <- min(opt$n_pcs, ncol(pc$x), na.rm = T)

rm(signals, pc) # Memory

# Regressions in parallel to speed things up

ncores <- floor(parallel::detectCores()*opt$cores_frac)
doParallel::registerDoParallel(cores = ncores)

pcr_pred <- foreach::"%dopar%"(foreach::foreach(
  j = seq(1,n_pcs,1), .packages = c('data.table','zoo')
), {

  # Need EW and VW regressions as separate models
  mod_ew <- lm(paste0("bh1m ~ ", paste(paste0("PC", 1:j), collapse = "+")), 
    reg_data[time_avail_m < pred_mon])
  mod_vw <- lm(paste0("bh1m ~ ", paste(paste0("PC", 1:j), collapse = "+")), 
    reg_data[time_avail_m < pred_mon], 
    weights = me) # Weight by the market cap

  # Get separate predictions
  pred_ew <- predict(mod_ew, reg_data[time_avail_m == pred_mon])
  pred_vw <- predict(mod_vw, # weighted predictions need to remove any NA weights
    reg_data[time_avail_m == pred_mon & !is.na(me)])
  wgts <- reg_data[time_avail_m == pred_mon & !is.na(me), me]

  # Remove because they take up a lot of memory
  rm(mod_ew, mod_vw)

  # Actual returns
  rets_ew <- reg_data[time_avail_m == pred_mon, bh1m]
  rets_vw <- reg_data[time_avail_m == pred_mon & !is.na(me), bh1m]

  # Get upper and lower quantiles for long-short portfolio formation
  upper_ew <- quantile(pred_ew, probs = 1 - opt$quantile_prob, na.rm = T)
  lower_ew <- quantile(pred_ew, probs = opt$quantile_prob, na.rm = T)
  upper_vw <- quantile(pred_vw, probs = 1 - opt$quantile_prob, na.rm = T)
  lower_vw <- quantile(pred_vw, probs = opt$quantile_prob, na.rm = T)

  cat("Finished PC", j, "\n")

  # Return 1-row data table of long-short returns
  data.table(
    ew_ls = mean(rets_ew[pred_ew > upper_ew], na.rm = T) - 
      mean(rets_ew[pred_ew < lower_ew], na.rm = T),
    vw_ls = weighted.mean(rets_vw[pred_vw > upper_vw], 
            w = wgts[pred_vw > upper_vw], na.rm = T) - 
        weighted.mean(rets_vw[pred_vw < lower_vw], 
            w = wgts[pred_vw < lower_vw], na.rm = T),
    pc = j # Labelling
  )

})

#==============================================================================# 
# Output ----
#==============================================================================#

if (opt$scaled_pca) {
    if (opt$scaled_pca_weight == "ew") {
        fcast <- "spca1"
    } else {
        fcast <- "spca2"
    }
} else {
    fcast <- "pca"
}

# Automatic output folder based on dataset
outfolder <- paste0( # Name based on forecast and imputation
    fcast, "_", 
    str_remove(str_remove(basename(opt$signal_file), '.csv'), "bcsignals_")
)

dir.create(paste0(FILEPATHS$out_path, "pca_returns/"), showWarnings = FALSE)
dir.create(paste0(FILEPATHS$out_path, "pca_returns/", outfolder),
    showWarnings = FALSE)

fwrite(
    rbindlist(pcr_pred)[, n_signals := length(signals_keep)], 
    paste0(
        FILEPATHS$out_path, 
        "pca_returns/",
        outfolder,
        '/ret_pc_', format(pred_mon, '%Y_%m'), '.csv'
    )
)


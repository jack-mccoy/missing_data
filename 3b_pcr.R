
#===============================================================================
# Packages
#===============================================================================

library(doParallel)
library(data.table)
library(foreach)
library(zoo)

#===============================================================================
# Option parsing
#===============================================================================

option_list <- list(
  optparse::make_option(c("--out_path"),
    type = "character", default = "./",
    help = "directory to store output to"),
  optparse::make_option(c("--signals_keep"),
    type = "character", default = "size,bm",
    help = "a comma-separated list of values or .txt file to scan"),
  optparse::make_option(c("--data_file"),
    type = "character", default = "mn_imputed_tmp.csv",
    help = "path of imputed file"),
  optparse::make_option(c("--iter_year"),
    type = "numeric", default = 1990,
    help = "year for making prediction (numeric)"),
  optparse::make_option(c("--iter_month"),
    type = "numeric", default = as.integer(Sys.getenv("SGE_TASK_ID")),
    help = "month for making prediction (numeric)"),
  optparse::make_option(c("--prefix"),
    type = "character", default = "pcr_",
    help = "prefix of output files"),
  optparse::make_option(c("--n_yrs"),
    type = "numeric", default = 7,
    help = "number of years to run for Fama-Macbeth principal component regressions"),
  optparse::make_option(c("--quantile_prob"),
    type = "numeric", default = 0.2,
    help = paste0("the ratio out of 1 used to form long-short portfolios ",
        "(i.e., 0.2 means quintile portfolios)")),
  optparse::make_option(c("--n_pcs"),
    type = "numeric", default = 25,
    help = "maximum number of principal components to use in PCRs")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Get the anomalies as a nice vector
if (grepl("\\.txt", opt$signals_keep)) {
  opt$signals_keep <- scan(opt$signals_keep, character())
} else if (grepl(",", opt$signals_keep)) {
  opt$signals_keep <- trimws(do.call("c", strsplit(opt$signals_keep, ",")))
} else {
  stop("It seems that you did not pass a .txt file or comma-separated list ",
    "to the `signals_keep` argument\n")
}

# Good directory ending
if (substr(opt$out_path, nchar(opt$out_path), nchar(opt$out_path)) != "/") {
  opt$out_path <- paste0(opt$out_path, "/")
}

#===============================================================================
# Functions
#===============================================================================

source("functions.R")

#===============================================================================
# Get prediction month from array iteration map
#===============================================================================

pred_mon <- as.yearmon(opt$iter_year + opt$iter_month - 1/12)

#===============================================================================
# Data pull
#===============================================================================

# Read in data and merge
signals <- fread(opt$data_file)
setnames(signals, colnames(signals), tolower(colnames(signals))) # ensure lower

crsp_data <- fread("crsp_data.csv")[, .(permno, yyyymm, ret, me)]
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
    (opt$signals_keep) := lapply(.SD, imputeVec, na.rm = T),
    .SDcols = opt$signals_keep,
    by = .(yyyymm)
]

# Complete cases are necessary for principal components. Ensuring here
signals <- signals[complete.cases(
  signals[, .SD, .SDcols = c("permno", "yyyymm", "ret", "me", opt$signals_keep)]
)]

#=============================================================================== 
# Timing adjustments
#===============================================================================

# Align signals with next month's returns
signals[
    order(permno, yyyymm),
    bh1m := shift(ret, 1, type = "lead"), # we want to predict next month's ret
    by = permno
]

# Define 'time_avail_m': the month of return prediction where signal was available
signals[, time_avail_m := yyyymm + 1/12]


#=============================================================================== 
# Run regressions
#===============================================================================

# Get principal components and form regression data
pc <- prcomp(signals[, .SD, .SDcols = opt$signals_keep], 
    center = TRUE, scale = TRUE)
reg_data <- data.table(signals[, .(permno, time_avail_m, bh1m, me)], pc$x)

# Max number of PCs we want to run
n_pcs <- min(opt$n_pcs, ncol(pc$x), na.rm = T)

rm(signals, pc) # Memory

# Regressions in parallel to speed things up
doParallel::registerDoParallel(cores = parallel::detectCores())
pcr_pred <- foreach::"%dopar%"(foreach::foreach(j = 1:n_pcs), {

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

#=============================================================================== 
# Output
#===============================================================================

fwrite(rbindlist(pcr_pred)[, n_signals := length(opt$signals_keep)], 
  paste0(opt$out_path, opt$prefix, 
    gsub("[[:space:]]", "", as.character(pred_mon)), ".csv"))



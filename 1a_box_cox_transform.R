# extended box-cox transforms signals
# outputs to opt$out_path/bcsignals_none.csv

#==============================================================================#
# Setup ----
#==============================================================================#

library(dplyr)
library(tidyr) # sadly
library(zoo)

source("functions.R")

closeAllConnections()

#==============================================================================#
# Option parsing ----
#==============================================================================#

option_list <- list(
    optparse::make_option(c("--cores_frac"),
        type = "numeric", 
        default = 0.5,
        help = "fraction of total cores to use"),
    optparse::make_option(c("--impute_vec"),
        type = "character", 
        default = "../data/signals_best150_1985.txt",
        help = "a comma-separated list of values or .txt file to scan, by default top 150 continuous and gapless"),
    optparse::make_option(c("--data_path"),
        type = "character",
        default = "../data/",
        help = "path for input data"),
    optparse::make_option(c("--out_path"),
        type = "character",
        default = "../data/",
        help = "path for output data")
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

#==============================================================================#
# Read in data ----
#==============================================================================#

# Read in the signals we want 
signals <- fread(paste0(opt$data_path, "signed_predictors_dl_wide_filtered.csv")) 
setnames(signals, colnames(signals), tolower(colnames(signals)))

# Ease of use
signals[, yyyymm := as.yearmon(yyyymm)]

# merge on crsp predictors
crsp_data <- fread(paste0(opt$data_path, "crsp_data.csv"))
crsp_data[, yyyymm := as.yearmon(yyyymm)]

signals <- merge(signals, crsp_data, by = c("permno", "yyyymm"))

# keep only those requested in impute_vec (which should be called something more helpful)
signals <- signals[, .SD, .SDcols = c("permno", "yyyymm", opt$impute_vec)]

# Memory
rm(crsp_data)

#==============================================================================#
# Box-Cox transformations and scaling ----
#==============================================================================#

# Get vector of time periods to run through 
# (ensuring that we're only using yrmons in dataset)
yrmons <- unique(signals$yyyymm)
yrmons <- yrmons[order(yrmons)]

# Run the Box-Cox transformations 

start_b <- Sys.time()

ncores = floor(parallel::detectCores()*opt$cores_frac)
doParallel::registerDoParallel(cores = ncores)
bctrans <- foreach::"%do%"(foreach::foreach(
  i = yrmons, .packages = c('data.table','zoo'), .combine = rbind
), {
  
    # A safety check to make sure we don't have any signals with a month of <2 obs
    signals_good <- names(which(sapply( # `which` filters out the `FALSE` results
        signals[
            yyyymm == i,
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
  
    # For checking in log
    # cat("There are a total of", length(signals_good), "signals with enough data.\n")
    # print(signals_good)
    # print(head(signals))

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
  ]
  
  return(transformed)
  
})

# Timing for the log
end_b <- Sys.time()
bc_time <- difftime(end_b, start_b, unit = "mins")
cat("Box-Cox run time:", bc_time, "minutes\n")

#==============================================================================#
# Output ----
#==============================================================================#

dir.create(opt$out_path, showWarning=FALSE)
fwrite(bctrans, paste0(opt$out_path, "bcsignals_none.csv"))

#==============================================================================#
# Checks to console ----
#==============================================================================#

sum_sd = bctrans[ , lapply(.SD, function(x) {sd(x, na.rm=T)})
                  , .SDcols = opt$impute_vec
                  , by = yyyymm
]

# show instances where SD is not 1
# it seems like it's pretty rare, looks good
sum_sd %>% 
  pivot_longer(cols = !c(yyyymm)) %>% 
  filter(!is.na(value)) %>% 
  filter(abs(value-1) > 0.1) %>% 
  distinct(name,value, .keep_all = T) %>% 
  print(n=200)

# bizarrely, herf has sd = 0 for a lot of the sample 
# it turns out the herf = 0 for many permnos before 1952
signals[!is.na(herf)  & (yyyymm - floor(yyyymm)) == 0.5
        , .(permno,yyyymm,herf) ]  %>%  
  arrange(yyyymm) %>% 
  print(topn = 20)

signals[!is.na(herf)  & (yyyymm - floor(yyyymm)) == 0.5, .(permno,yyyymm,herf) ][
  , .(n_zero = sum(herf == 0))
  , by = yyyymm
] %>% 
  arrange(yyyymm) %>% 
  print(topn= 10)

# check number of observations
n_bc = bctrans[ , lapply(.SD, function(x) {sum(!is.na(x))})] %>% 
  pivot_longer(cols = everything(), values_to = 'n_boxcox') %>% 
  print(n= 100)

n_raw = signals[ , lapply(.SD, function(x) {sum(!is.na(x))})] %>% 
  pivot_longer(cols = everything(), values_to = 'n_raw') %>% 
  print(n= 100)

n_raw %>% left_join(n_bc) %>% 
  mutate(diff = n_raw - n_boxcox) %>% 
  print(n=100)

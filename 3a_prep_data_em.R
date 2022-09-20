
#==============================================================================#
# Packages
#==============================================================================#

library(car)
library(doParallel)
library(data.table)
library(foreach)
library(zoo)

#==============================================================================#
# Option parsing
#==============================================================================#

option_list <- list(
    optparse::make_option(c("--impute_vec"),
        type = "character", default = "signals.txt",
        help = "a comma-separated list of values or .txt file to scan"),
    optparse::make_option(c("--sample_start_year"),
        type = "numeric", default = 1985,
        help = "year that sample starts"),
    optparse::make_option(c("--sample_end_year"),
        type = "numeric", default = 2020,
        help = "year that sample ends"),
    optparse::make_option(c("--tmp_file_imp"),
        type = "character", 
        default = "imputed_tmp.csv",
        help = "name of temporary output file for imputed dataset"),
    optparse::make_option(c("--tmp_file_bc"),
        type = "character", 
        default = "bc_tmp.csv",
        help = "name of temporary output file for imputed dataset"),
    optparse::make_option(c("--params_path"),
        type = "character", 
        default = "./",
        help = "directory holding parameter estimates")
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

# Good directory ending
if (substr(opt$params_path, nchar(opt$params_path), nchar(opt$params_path)) != "/") {
  opt$params_path <- paste0(opt$params_path, "/")
}

#==============================================================================#
# Functions
#==============================================================================#

source("functions.R")

#==============================================================================#
# Data pull
#==============================================================================#

# Read in the signals 
signals <- fread("../data/signed_predictors_dl_wide.csv") 
setnames(signals, colnames(signals), tolower(colnames(signals)))

# Ease of use
signals[, yyyymm := as.yearmon(as.Date(as.character(yyyymm*10 + 1), "%Y%m%d"))]

# Read in the other data from CRSP 
crsp_data <- fread("../data/crsp_data.csv")
crsp_data[, yyyymm := as.yearmon(yyyymm)]

signals <- merge(signals, crsp_data, by = c("permno", "yyyymm"))[, 
    .SD, 
    .SDcols = c("permno", "yyyymm", opt$impute_vec)
]

#==============================================================================#
# Imputation
#==============================================================================#

# Sequence of yearmons to impute
yrmons <- seq(
    as.yearmon(paste0("Jan ", opt$sample_start_year)),
    as.yearmon(paste0("Dec ", opt$sample_end_year)),
    by = 1/12
)

doParallel::registerDoParallel(cores = parallel::detectCores())
imputed <- foreach::"%dopar%"(foreach::foreach(i = yrmons), {

  i_file <- gsub("[[:space:]]", "", i)

  # Get the Box-Cox and scaling parameters ----

  params <- fread(paste0(opt$params_path, "bcn_scale_", i_file, ".csv"))

  # Get the imputation parameters ----

  estE <- as.matrix(read.csv(paste0(opt$params_path, "estE_", i_file, ".csv"),
    row.names = 1))[, 1][opt$impute_vec] # Need as vector with names, not 1 column matrix
  estR <- as.matrix(read.csv(paste0(opt$params_path, "estR_", i_file, ".csv"),
    row.names = 1))[opt$impute_vec, opt$impute_vec]

  # Select data for given month
  good <- names(checkMinObs( # Get columns we can use
    signals[yyyymm == i, .SD, .SDcols = opt$impute_vec],
    min_obs = 2
  ))
  tmp <- signals[yyyymm == i]
  
  # Perform the Box-Cox transformation and scaling ----

  tmp[,
    (good) := foreach::"%do%"(foreach::foreach(j = good), {

      lambda_j <- params[variable == j & param == "bcn:lambda", value]
      gamma_j <- params[variable == j & param == "bcn:gamma", value]
      center_j <- params[variable == j & param == "scaled:center", value]
      scale_j <- params[variable == j & param == "scaled:scale", value]

      # Winsorize 
      x <- winsorize(tmp[, get(j)], tail = 0.005)
      
      # BC-transform if possible
      if (!is.na(lambda_j) & !is.na(gamma_j)) {
        x <- car::bcnPower(x, lambda = lambda_j, gamma = gamma_j)
      }

      # returning scaled value exactly as we had before
      (x - center_j) / scale_j
    }),
    .SDcols = good
  ]

  # Deal with negative eigenvalues here ----

  R0 <- cov(tmp[, .SD, .SDcols = good], use = "pairwise.complete.obs")
  id <- diag(nrow = nrow(R0), ncol = ncol(R0))

  # Error checking. Catch missing means and covariance
  if (any(is.na(R0))) R0[is.na(R0)] <- id[is.na(R0)]

  # Ensure covariance matrix is positive semi-definite ----

  eig <- eigen(R0)
  if (any(eig$values < 0)) {
    k <- 1
    while (any(eig$values < 0)) {
      R0 <- structure(
        eig$vectors %*% diag(ifelse(eig$values <= 0, 0, eig$values)) %*% t(eig$vectors),
        dimnames = list(rownames(R0), colnames(R0))
      )
      eig <- eigen(R0)
      if (k == 100) break
    }
    tmp[, # Now make it so that data has same variance as new cov mat
      (good) := foreach::"%do%"(foreach::foreach(j = good), {
        tmp[, get(j)] * R0[j, j]
      })
    ]
  } 

  # Impute the data from estimated parameters ----

  imp_i <- mvn_emf(as.matrix(tmp[, .SD, .SDcols = good]), 
    E0 = estE[good], R0 = estR[good, good], 
    tol = 1e-4, maxiter = 1, update_estE = FALSE)$Ey # Most of these should be converging quickly  

  colnames(imp_i) <- good # set the correct names

  imp_i <- data.table(imp_i, permno = tmp$permno, yyyymm = tmp$yyyymm)

  cat("Imputed for", as.character(i), "\n")

  # Return the imputed and un-transformed data
  list(imp = imp_i, bc = tmp[, .SD, .SDcols = c("permno", "yyyymm", good)])

})

#==============================================================================#
# Output
#==============================================================================#

fwrite(rbindlist(lapply(imputed, function(x) return(x$imp)), fill = T), 
    opt$tmp_file_imp)
fwrite(rbindlist(lapply(imputed, function(x) return(x$bc)), fill = T),
    opt$tmp_file_bc)


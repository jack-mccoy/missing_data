# inputs: box-cox and (if desired) mean and cov estimates 
# outputs: big permno-month-signal datasets
# options to use different kinds of imputation (or none)

#==============================================================================#
# Packages ----
#==============================================================================#

library(car)
library(doParallel)
library(data.table)
library(foreach)
library(zoo)

#==============================================================================#
# Option parsing ----
#==============================================================================#

option_list <- list(
    optparse::make_option(c("--impute_type"),
        type = "character", 
        default = "em",
        help = "type of imputation: (none, em, availcase)"),  
    optparse::make_option(c("--impute_vec"),
        type = "character",
        default = "../output/signals.txt",
        help = "a comma-separated list of values or .txt file to scan"),
    optparse::make_option(c("--sample_start_year"),
        type = "numeric",
        default = 1985,
        help = "year that sample starts"),
    optparse::make_option(c("--sample_end_year"),
        type = "numeric",
        default = 2020,
        help = "year that sample ends"),
    optparse::make_option(c("--bcsignals_filename"),
        type = "character", 
        default = "../output/bcsignals/bcsignals_em.csv",
        help = "name of temporary output file for dataset"),
    optparse::make_option(c("--params_path"),
        type = "character", 
        default = "../output/impute_ests",
        help = "directory holding parameter estimates"),
    optparse::make_option(c("--cores_frac"),
        type = "numeric",
        default = 1.0,
        help = "fraction of total cores to use")        
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

# default output folder 
dir.create('../output/bcsignals/', showWarnings = F)

#==============================================================================#
# Functions ----
#==============================================================================#

source("functions.R")

#==============================================================================#
# Data pull ----
#==============================================================================#

# Read in the signals 
signals <- fread("../data/signed_predictors_dl_wide.csv") 
setnames(signals, colnames(signals), tolower(colnames(signals)))

# Ease of use
signals[, yyyymm := as.yearmon(as.Date(as.character(yyyymm*10 + 1), "%Y%m%d"))]

# Read in the other data from CRSP 
crsp_data <- fread("../data/crsp_data.csv")
crsp_data[, ":="(
    yyyymm = as.yearmon(yyyymm),
    sic3 = substr(ifelse(hsiccd < 1000, paste0("0", hsiccd), as.character(hsiccd)), 
        1, 3)
)]

signals <- merge(signals, crsp_data, by = c("permno", "yyyymm"))[, 
    .SD, 
    .SDcols = c("permno", "yyyymm", "sic3", "me", opt$impute_vec)
]

#==============================================================================#
# Make new signals data ----
#==============================================================================#

# Sequence of yearmons to impute
yrmons <- seq(
    as.yearmon(paste0("Jan ", opt$sample_start_year)),
    as.yearmon(paste0("Dec ", opt$sample_end_year)),
    by = 1/12
)

# set up parallel comp
ncores = floor(parallel::detectCores()*opt$cores_frac)
doParallel::registerDoParallel(cores = ncores)

signals_new_list <- foreach::"%dopar%"(foreach::foreach(
  i = yrmons,  .packages = c('data.table','zoo')
  ), {

    i_file <- gsub("[[:space:]]", "", i)
    
    # Select data for given month
    good <- names(checkMinObs( # Get columns we can use
        signals[yyyymm == i, .SD, .SDcols = opt$impute_vec],
        min_obs = 2
    ))
    
    signals_tmp <- signals[
        yyyymm == i,
        .SD,
        .SDcols = c("permno", "yyyymm", "sic3", "me", good)
    ]

    ## Box-Cox transform ----
    
    # Get the Box-Cox and scaling parameters 
    params <- fread(paste0(opt$params_path, "bcn_scale_", i_file, ".csv"))

    # Perform the Box-Cox transformation and scaling 
    signals_tmp[,
        (good) := foreach::"%do%"(foreach::foreach(j = good), {
            
            lambda_j <- params[variable == j & param == "bcn:lambda", value]
            gamma_j <- params[variable == j & param == "bcn:gamma", value]
            center_j <- params[variable == j & param == "scaled:center", value]
            scale_j <- params[variable == j & param == "scaled:scale", value]

            # Winsorize 
            x <- winsorize(signals_tmp[, get(j)], tail = 0.005)
            
            # BC-transform if possible
            if (!is.na(lambda_j) & !is.na(gamma_j)) {
                x <- car::bcnPower(x, lambda = lambda_j, gamma = gamma_j)
            }

            # returning scaled value exactly as we had before
            (x - center_j) / scale_j
        }),
        .SDcols = good
    ]

    # Imputation ----
    if (opt$impute_type %in% c('em', 'availcase')){
        
        # Get the imputation parameters (depending on em or availcase)
        if (opt$impute_type == 'em'){
            # em imputation settings
            estE <- as.matrix(read.csv(paste0(opt$params_path, "estE_", i_file, ".csv"),
                                       row.names = 1))[, 1][opt$impute_vec] # Need as vector with names, not 1 column matrix
            estR <- as.matrix(read.csv(paste0(opt$params_path, "estR_", i_file, ".csv"),
                                       row.names = 1))[opt$impute_vec, opt$impute_vec]
            
            # maxiter = 1 should work, but for sanity 10 is better due to rounding issues
            maxiter = 10 
        } else if (opt$impute_type == 'availcase'){
            # availcase settings
            # this is just: impute using sample moments (1 iteration)
            tmpmat = signals_tmp[ , ..good]
            estE = colMeans(tmpmat, na.rm=T)
            estR = cov(tmpmat, use = 'pairwise.complete.obs')
            
            # use Highham's method to ensure PSD
            # e.g. https://github.com/cran/pracma/blob/c1688b374d201c13fb40b4dda2d2a89e34b94ec6/R/nearest_spd.R
            eig <- eigen(estR)
            if (any(eig$values < 0)) {
                k <- 1
                while (any(eig$values < 0)) {
                    estR <- structure(
                        eig$vectors %*%
                            diag(ifelse(eig$values <= 0, 0, eig$values)) %*%
                            t(eig$vectors),
                        dimnames = list(rownames(estR), colnames(estR))
                    )
                    eig <- eigen(estR)
                    if (k == 100) break
                }
            } 
            
            # avail case is just 1 iter of EM
            maxiter = 1
        }
        
        # Impute the data from estimated parameters ----
        imp_i <- mvn_emf(as.matrix(signals_tmp[, .SD, .SDcols = good]), 
            E0 = estE[good], R0 = estR[good, good], 
            tol = 1e-4, maxiter = maxiter, update_estE = FALSE)$Ey # Most of these should be converging quickly  
        
        colnames(imp_i) <- good # set the correct names
        
        imp_i <- data.table(
            imp_i, 
            permno = signals_tmp$permno,
            yyyymm = signals_tmp$yyyymm
        )
        
        # overwrite signals_tmp
        signals_tmp = imp_i
        
    } else if (opt$impute_type == "ind") {
        # deep copy
        signals_tmp_dup <- copy(signals_tmp)
        
        # Get industry-level and overall cross-sectional means
        signals_tmp[,
            (good) := lapply(.SD, imputeVec),
            .SDcols = good,
            by = .(yyyymm, sic3)
        ]
        signals_tmp_dup[, (good) := lapply(.SD, imputeVec), .SDcols = good]
        
        # If there was no data for industry, replace with overall mean
        signals_tmp_mat <- as.matrix(signals_tmp[, ..good])
        signals_tmp_dup_mat <- as.matrix(signals_tmp_dup[, ..good])
        signals_tmp_mat[is.na(signals_tmp_mat)] <- signals_tmp_dup_mat[is.na(signals_tmp_mat)]
        
        # Final for output
        signals_tmp <- data.table(
            signals_tmp_mat, 
            permno = signals_tmp$permno,
            yyyymm = signals_tmp$yyyymm
        )
    } else if (opt$impute_type == "indsize") {
        # deep copy
        signals_tmp_dup <- copy(signals_tmp)

        # Get industry X size, industry-level, and overall cross-sectional means
        signals_tmp[,
            indsize_decile := tryCatch(pctileGroups(me, pctile = 0.1),
                error = function(e) return(0)), # if can't make unique breaks, all together
            by = .(yyyymm, sic3)
        ][,
            (good) := lapply(.SD, imputeVec),
            .SDcols = good,
            by = .(yyyymm, sic3, indsize_decile)
        ]
        signals_tmp <- signals_tmp[, !c("indsize")] # remove so standard with dup
        signals_tmp_dup[, (good) := lapply(.SD, imputeVec), .SDcols = good]
        
        # If there was no data for industry, replace with overall mean
        signals_tmp_mat <- as.matrix(signals_tmp[, ..good])
        signals_tmp_dup_mat <- as.matrix(signals_tmp_dup[, ..good])
        signals_tmp_mat[is.na(signals_tmp_mat)] <- signals_tmp_dup_mat[is.na(signals_tmp_mat)]
        
        # Final for output
        signals_tmp <- data.table(
            signals_tmp_mat, 
            permno = signals_tmp$permno,
            yyyymm = signals_tmp$yyyymm
        )
    } # end impute if desired
    
    cat("signals_new for", as.character(i), "\n")

    # Return the imputed and un-transformed data
    return(signals_tmp)
})

# Combine 
out <- rbindlist(signals_new_list, fill = T)

# Catch the last value imputations if specified (works differently)
if (opt$impute_type == "lastval") {
    out <- out[order(permno, yyyymm)] # order in advance
    out2 <- copy(out) # for cross-sectional means
    out[, # Use last value for imputation, up to 12 month lag
        (opt$impute_vec) := lapply(.SD, function(x) imputeLastVal(x, yyyymm, k = 12)),
        .SDcols = opt$impute_vec,
        by = .(permno)
    ]
    out2[, # Cross-sectional mean imputation
        (opt$impute_vec) := lapply(.SD, function(x) imputeVec(x)),
        .SDcols = opt$impute_vec,
        by = .(yyyymm)
    ]

    # Ensure same ordering
    out <- out[order(permno, yyyymm)]
    out2 <- out2[order(permno, yyyymm)]

    # If could not impute with last val, use XS mean. Indexing easier with matrix
    out_mat <- as.matrix(out[, .SD, .SDcols = opt$impute_vec])
    out2_mat <- as.matrix(out2[, .SD, .SDcols = opt$impute_vec])
    out_mat[is.na(out_mat)] <- out2_mat[is.na(out_mat)]
    
    # Final for output
    out <- data.table(
        out_mat, 
        permno = out$permno,
        yyyymm = out$yyyymm
    )
}

#==============================================================================#
# Output ----
#==============================================================================#

fwrite(out, opt$bcsignals_filename)


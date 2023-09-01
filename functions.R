#===============================================================================#
# Packages
#===============================================================================#

library(car)          # necessary in Box-Cox transformations
library(data.table)   # base of all the data manipulation in project
library(lubridate)    # helpful date/time formatting in project
library(doParallel)   # fast split/apply/combine analysis
library(ggplot2)      # missingness plot

#===============================================================================#
# Functions
#===============================================================================#

#-------------------------------------------------------------------------------#
# Wrapper for mean-variance EM in Fortran
#-------------------------------------------------------------------------------#
mvn_emf <- function(y, E0, R0, update_estE = TRUE, 
    tol = 1e-6, maxiter = 1e3) {
  
  if (!is.loaded('mvn_emf')) {
    if(.Platform$OS.type == "unix") {
      dyn.load("mvn_emf.so")
    } else {
      dyn.load("mvn_emf_win.so")
    }
  }
  
  N <- nrow(y)
  K <- ncol(y)
  nan.index <- is.na(y)
  y[nan.index] <- 0
  
  fvals <- .Fortran(
    "mvn_emf", Ey = as.double(y),
    estE = as.double(E0), estR = as.double(R0),
    nan.index = as.logical(nan.index),
    tol = as.double(tol), maxiter = as.integer(maxiter),
    K = as.integer(K), N = as.integer(N),
    update_estE = update_estE
  )
  
  fvals$Ey <- matrix(fvals$Ey,N,K)
  fvals$estR <- matrix(fvals$estR,K,K)
  fvals$nan.index <- matrix(fvals$nan.index,N,K)

  return(fvals)
  
}

#-------------------------------------------------------------------------------#
# "not in" binary
#-------------------------------------------------------------------------------#
"%ni%" <- Negate("%in%")

#-------------------------------------------------------------------------------#
# Pctile group classification 
#-------------------------------------------------------------------------------#
pctileGroups <- function(vec, pctile = 0.2) {
  return(cut(vec,
    breaks = unname(quantile(vec ,probs = seq(0, 1, by = pctile), na.rm = T)),
    include.lowest = T, labels = F
  ))
}

#-------------------------------------------------------------------------------#
# Checking for consistent matrix columns 
#-------------------------------------------------------------------------------#
sameMatCol <- function(mat_list) {
  if (any(do.call("cbind", lapply(mat_list, class)) != "matrix")) {
    stop("All elements of `mat_list` must be matrices.\n")
  }
  if (length(unique((do.call("rbind", lapply(mat_list, ncol))))) != 1) {
    stop(paste0(
      "All matrices in `mat_list` must have same number of columns (variables). ",
      "Also ensure that each matrix has named columns or columns ",
      "following the same order. The sequential imputation assumes this ",
      "to assist in error checking\n."
    ))
  }
}

#-------------------------------------------------------------------------------#
# checkMinObs 
#-------------------------------------------------------------------------------#
checkMinObs <- function(mat, min_obs) {
  return(which(apply(mat, MARGIN = 2, FUN = function(x) {
    length(x[!is.na(x)]) > min_obs
  })))
}

#-------------------------------------------------------------------------------#
# Simple imputation of a vector 
#-------------------------------------------------------------------------------#
imputeVec <- function(x, na.rm = T) {
  x[is.na(x)] <- mean(x, na.rm = na.rm)
  return(x)
}

#-------------------------------------------------------------------------------#
# Winsorizing data
#-------------------------------------------------------------------------------#
winsorize <- function(vec, tail = 0.005, upper = 1-tail, lower = tail) {
  if((upper + lower) > 1){
    stop("Sum of tails greater than 1")
  }
  lower_val <- unname(quantile(vec, probs = lower, na.rm = T))
  upper_val <- unname(quantile(vec, probs = upper, na.rm = T))
  vec[which(vec < lower_val)] <- lower_val
  vec[which(vec > upper_val)] <- upper_val
  return(vec)
}

#-------------------------------------------------------------------------------#
# Function: getBetas
#  Takes covariance matrix and empty dataset
#  Calculates EMF betas for each observation
#  applying to each row of observed data
#  we select the missing and not-missing columns
#  we then calculate betas
#  using the partitioned covariance matrix
#-------------------------------------------------------------------------------#
getBetas <- function(covs, empties){
  beta_list <- apply(
    empties, FUN = function(x) {
      missing <- which(x)
      observed <- which(!x)
      sigma_12 <- covs[missing,observed]
      sigma_22 <- covs[observed,observed]
      sigma_22_inv <- solve(sigma_22)
      betas <- sigma_12 %*% sigma_22_inv
      return(betas)
    },
    MARGIN = 1
  )
}

#-------------------------------------------------------------------------------#
# Function: wtdMean
#  Enables weighed mean with na.rm for `x` and `w`
#-------------------------------------------------------------------------------#
wtdMean <- function(x, w, ..., na.rm = F) {
  x_obs <- !is.na(x)
  w_obs <- !is.na(w)
  if (na.rm) {
    out <- weighted.mean(
      x[which(x_obs & w_obs)], w[which(x_obs & w_obs)],
      ..., na.rm = na.rm
    )
  } else {
    out <- weighted.mean(x, w, ...)
  }
  return(out)
}

#-------------------------------------------------------------------------------#
# Function: na2zero
#  Replaces NAs in vector with 0
#-------------------------------------------------------------------------------#
na2zero <- function(x) {
  return(ifelse(is.na(x), 0, x))
}

#-------------------------------------------------------------------------------#
# Function: lsos
#  lists objects by size
#-------------------------------------------------------------------------------#

# improved list of objects
.ls.objects <- function(pos = 1, pattern, order.by,
  decreasing = FALSE, head = FALSE, n = 5
) {
  napply <- function(names, fn) {
    sapply(names, function(x) fn(get(x, pos = pos)))
  }
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.prettysize <- napply(names, function(x) {
    format(utils::object.size(x), units = "auto") 
  })
  obj.size <- napply(names, object.size)
  obj.dim <- t(napply(names, function(x) as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
  names(out) <- c("Type", "Size", "PrettySize", "Length/Rows", "Columns")
  if (!missing(order.by)) {
    out <- out[order(out[[order.by]], decreasing=decreasing), ]
  }
  if (head) {
    out <- head(out, n)
  }
  return(out)
}

# shorthand
lsos <- function(..., n = 10) {
  .ls.objects(..., order.by = "Size", decreasing=TRUE, head=TRUE, n=n)
}

#-------------------------------------------------------------------------------#
# Read in a specific file from a zip using function of your choice
#-------------------------------------------------------------------------------#
readZipFile <- function(zipfile, tmp_dir, datafile, FUN, prefix = "", ...) {
  if (!any(system(paste0("unzip -l ", zipfile), intern = T) %like% datafile)) {
    stop(datafile, " not found in ", zipfile, "\n")
  }
  if (!("function" %in% class(FUN)) & class(FUN) == "character") {
    if (!exists(FUN)) stop("`", FUN, "` not found in workspace\n")
  }
  system(paste0("unzip -p ", zipfile, " ",  datafile,
    " > ", tmp_dir, prefix, datafile))
  dat <- do.call(FUN, list(paste0(tmp_dir, prefix, datafile), ...))
  system(paste0("rm -f ", tmp_dir, prefix, datafile))
  return(dat)
}

#-------------------------------------------------------------------------------#
# Convert stata YYYYmM format to zoo's yearmon format
#-------------------------------------------------------------------------------#
stataMon2yearmon <- function(statamon) {
  if (!all(statamon %like% "m")) { 
    stop("`statamon` must be Stata-formatted year-month, i.e. '1996m1', '2019m12'\n")
  }
  zoo::as.yearmon(sapply(base::strsplit(statamon, "m"), function(x) {
    as.numeric(x[1]) + (as.numeric(x[2]) - 1) / 12 
  }))
}

#-------------------------------------------------------------------------------#
# Nearest PSD matrix a la N.J. Highham thanks to pracma
# https://github.com/cran/pracma/blob/c1688b374d201c13fb40b4dda2d2a89e34b94ec6/R/nearest_spd.R
#-------------------------------------------------------------------------------#

nearest_spd <- function(A) {
  stopifnot(is.numeric(A), is.matrix(A))
  eps <- .Machine$double.eps
  
  m <- nrow(A); n <- ncol(A)
  if (m != n) {
    stop("Argument 'A' must be a square matrix.")
  } else if (n == 1 && A <= 0)
    return(as.matrix(eps))
  
  B <- (A + t(A)) / 2                 # symmetrize A
  svdB <- svd(B)                      # H is symmetric polar factor of B
  H <- svdB$v %*% diag(svdB$d) %*% t(svdB$v)
  
  Ahat <- (B + H) / 2
  Ahat <- (Ahat + t(Ahat)) / 2
  
  # Test that Ahat is in fact positive-definite;
  # if it is not so, then tweak it just a bit.
  k <- 0; not_pd <- TRUE
  while (not_pd) {
    k <- k + 1
    try_R <- try(chol(Ahat), silent = TRUE)
    if(inherits(try_R, "try-error")) {
      mineig <- min(eigen(Ahat, symmetric = TRUE, only.values = TRUE)$values)
      Ahat = Ahat + (-mineig*k^2 + eps(mineig)) * diag(1, n)
    } else
      not_pd <- FALSE
  }
  Ahat
}

# Figures Setup -----------------------------------------------------------------

MATBLUE = rgb(0,0.4470,0.7410)
MATRED = rgb(0.8500, 0.3250, 0.0980)
MATYELLOW = rgb(0.9290, 0.6940, 0.1250)

NICEBLUE = "#619CFF"
NICEGREEN = "#00BA38"
NICERED = "#F8766D"

chen_theme =   theme_minimal() +
  theme(
    text = element_text(family = "Palatino Linotype")
    , panel.border = element_rect(colour = "black", fill=NA, size=1)    
    
    # Font sizes
    , axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26),
    axis.text.x = element_text(size = 26),
    axis.text.y = element_text(size = 26),
    legend.text = element_text(size = 22),
    
    # Tweaking legend
    legend.position = c(0.7, 0.8),
    legend.text.align = 0,
    legend.background = element_rect(fill = "white", color = "black"),
    legend.margin = margin(t = 5, r = 2, b = 5, l = 5), 
    legend.key.size = unit(1.8, "cm"), 
    legend.spacing.x = unit(0.001, 'cm'),
    legend.title = element_blank()    
  ) 

imputeLastVal <- function(x, yrmons, k = 12) {
    if (length(x) != length(yrmons)) {
        stop('`x` and `yrmons` must have same length\n')
    }
    if (length(x) == 1) {
        return(x)
    } else {
        if (any((yrmons - shift(yrmons))[2:length(yrmons)] <= 0)) {
            stop('`yrmons` must be strictly increasing\n')
        }
        y <- x
        for (i in 1:length(x)) {
            if (is.na(y[i])) {
                x_prev <- x[which(
                    ((yrmons[i] - k * (1/12)) < yrmons) &
                    (yrmons < yrmons[i])
                )]
                if (any(!is.na(x_prev))) {
                    non_nas <- x_prev[which(!is.na(x_prev))]
                    y[i] <- last(non_nas)
                }
            }
        }
        return(y)
    }
}

getFilePaths <- function(signal_list_exists=TRUE) {
    pathnames <- c("data_path", "out_path", "signal_list")
    # Load file
    if (file.exists("FILEPATHS.R")) {
        source('FILEPATHS.R')
        # Make sure it has a list `FILEPATHS`
        if (exists("FILEPATHS")) {
            # Make sure it has the right names
            if (all(pathnames %in% names(FILEPATHS))) {
                # If desired, make sure signal list exists (created after a few scripts)
                if (signal_list_exists) {
                    if (!file.exists(FILEPATHS$signal_list)) {
                        stop("signal_list `", FILEPATHS$signal_list, "` does not exist\n", sep="")
                    }
                }
                for (path in pathnames[1:2]) {
                    if (!dir.exists(FILEPATHS[[path]])) {
                        cat("Creating directory `", FILEPATHS[[path]], "`\n", sep = "")
                        dir.create(FILEPATHS[[path]], showWarnings = FALSE)
                    }
                }
            } else {
                stop(
                    "`FILEPATHS` list does not contain all of (", 
                    paste(pathnames, collapse = ", "), ")\n",
                    sep = ""
                )
            }
        } else {
            stop(
                "File `FILEPATHS.R` must have list `FILEPATHS` and entries (",
                paste(pathnames, collapse = ", "), ")\n",
                sep = ""
            )
        }
    } else {
        stop(
            "Must make `FILEPATHS.R` file with list `FILEPATHS` and entries (",
            paste(pathnames, collapse = ", "), ")\n",
            sep = ""
        )
    }
}

unpackSignalList <- function(signal_list) {
    if (grepl("\\.txt", signal_list)) {
        impute_vec <- scan(signal_list, character())
    } else if (grepl(",", signal_list)) {
        impute_vec <- trimws(do.call("c", strsplit(signal_list, ",")))
    } else {
        stop("It seems that `signal_list` is not a .txt file or comma-separated list\n")
    }
    return(impute_vec)
}

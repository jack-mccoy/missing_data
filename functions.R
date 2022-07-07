#==============================================================================#
#  
# Functions for anomalies project w/ Andrew Chen & Fabian Winkler
# Created: August 2018
# Updated: all the time
#
#==============================================================================#

#==============================================================================#
# Packages
#==============================================================================#

# trying to limit required packages
library(car)          # necessary in Box-Cox transformations
library(data.table)   # base of all the data manipulation in project
library(lubridate)    # helpful date/time formatting in project
library(doParallel)   # fast split/apply/combine analysis
library(ggplot2)      # missingness plot
#library(gridExtra)    # missingness plot
#library(kableExtra)   # LaTeX tables
#library(xlsx)         # for header info

#==============================================================================#
# Functions
#==============================================================================#

#-------------------------------------------------------------------------------
# Function: getGitRepo
#-------------------------------------------------------------------------------
if (!(exists("getGitRepo"))) {
  getGitRepo <- function() {
    tmp <- system("git rev-parse --show-toplevel", intern = TRUE)
    if (!is.null(attr(tmp, "status"))) {
      if (attr(tmp, "status") == 128) stop("You're not in your git repo.\n") 
    }
    return(tmp)
  }
}

#-------------------------------------------------------------------------------
# Function: makePos
#-------------------------------------------------------------------------------
makePos <- function(x, sd_mult = 1, na.rm = TRUE, 
  cond = !(length(which(is.na(x))) == length(x))
) {  
  if (eval(substitute(cond))) {
    x <- x - min(x, na.rm = na.rm)
    return(x + min(x[x > sd_mult*sd(x, na.rm = na.rm)], na.rm = na.rm))
  } else { 
    return(x)
  }
}

#-------------------------------------------------------------------------------
# Function: boxCoxTransform
#-------------------------------------------------------------------------------
boxCoxTransform <- function(x, fmla = x~1, param = "lambda", 
  cond = !(length(which(is.na(x))) == length(x))
) {
  if (eval(substitute(cond))) {
    lam <- unname(powerTransform(as.formula(fmla))[[param]])
    return(bcPower(x, lam))
  } else {
    return(x)
  }
}

#-------------------------------------------------------------------------------
# Wrapper for mean-variance EM in Fortran
#-------------------------------------------------------------------------------
mvn_emf <- function(y, E0, R0, tol = 1e-6, maxiter = 1e3) {
  
  if (!is.loaded('mvn_emf')) {
    dyn.load("mvn_emf.so")
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
    K = as.integer(K), N = as.integer(N)
  )
  
  fvals$Ey <- matrix(fvals$Ey,N,K)
  fvals$estR <- matrix(fvals$estR,K,K)
  fvals$nan.index <- matrix(fvals$nan.index,N,K)

  return(fvals)
  
}

#-------------------------------------------------------------------------------
# "not in" binary
#-------------------------------------------------------------------------------
"%ni%" <- Negate("%in%")

#-------------------------------------------------------------------------------
# simple moving average function
#-------------------------------------------------------------------------------
movavg <- function(x, window)
{
  output <- x

  ## starting in first obs with enough previous data, calculate moving average
  for(i in window:length(x)) {
    output[i] <- mean(x[(i - window + 1):i], na.rm = TRUE)
  }

  ## subset to only include observations with enough data available
  output <- output[window:length(output)]
  
  ## return the moving average
  return(output)
}

#-------------------------------------------------------------------------------
# Get complete columns
#-------------------------------------------------------------------------------
completeCols <- function(data) {
  check <- do.call("rbind", lapply(data, function(x) !any(is.na(x))))
  singular <- do.call("rbind", lapply(data, function(x) {
    length(unique(x[!is.na(x)])) == 1 
  }))
  goodset <- rownames(check)[which(check & !singular)]
  return(goodset)
}

#-------------------------------------------------------------------------------
# Regressions by group subset
#-------------------------------------------------------------------------------
regressions <- function(data, dep_var, rhs_vars, group_var,
  date_adj = 1 ,weight_var = NULL
) {
  
  ## this function uses data.table syntax, so must be formatted appropriately
  if(!is.data.table(data)){
    data <- as.data.table(data)
  }

  ## order data by regression group and save row names
  row_names <- as.matrix(unique(data[,..group_var][order(get(group_var))]))

  ## regress model for each date subset
  registerDoParallel(cores = 14)
  coefs <- foreach(dat = 1:(length(row_names) - date_adj)) %dopar% {
    
    this_dat <- row_names[dat]
    temp <- data[get(group_var) == as.Date(this_dat)]
    
    ## find complete columns that we want to use in regressions for given month
    ## essentially prioritizes dropping variables over dropping observations
    goodcols <- completeCols(temp)
    colidx <- which(rhs_vars %in% goodcols) 
    these_rhs_vars <- rhs_vars[colidx]

    reg_spec <- paste0(dep_var, "~", paste(these_rhs_vars, collapse = "+"))
    
    ## regression
    if(!is.null(weight_var)){
      this_lm <- lm(as.formula(reg_spec),data = temp,weights = get(weight_var))
    } else {
      this_lm <- lm(as.formula(reg_spec),data = temp)
    }
    
    ## extract coefficients
    this_coef <- matrix(0,nrow=1,ncol=length(rhs_vars)+1)
    this_coef[1,c(1,colidx+1)] <- coef(this_lm)
    colnames(this_coef) <- c("Intercept",rhs_vars)

    ## label with the date of regression
    ## we want regressions to show data that we had available at time of predictions
    ## hence the "+ date_adj" when necessary to avoid look-ahead bias
    rownames(this_coef) <- row_names[dat + date_adj]

    ## return
    return(this_coef)

  }

  ## summary statistics to report
  coefmat <- do.call("rbind", coefs)
  
  ## for some subsets there are not enough DoF for all variables to have coefficients
  coefmat[which(is.na(coefmat))] <- 0

  ## return
  return(coefmat)

}

#-------------------------------------------------------------------------------
# Rolling FM betas 
#-------------------------------------------------------------------------------
fmBetas <- function(parms_df, window, parms, intercept = TRUE) { 
  
  ## including intercept if it is in data but not parameter vector
  if(intercept & ("(Intercept)" %ni% parms | "Intercept" %ni% parms)) {
    parms <- c("Intercept", parms)    
  }

  ## to filter out where we have enough data
  min_row <- window
  max_row <- length(rownames(parms_df))

  ## storing date labels for ease of lapply below
  beta_dates <- rownames(parms_df)[min_row:max_row]
  parms_df <- as.data.frame(parms_df)

  ## creating rolling average parameters dataset and housekeeping
  parms_df <- do.call("cbind", lapply(parms_df, function (x) {
    movavg(x, window = window)
  }))
  
  ## some tidying up
  if(intercept) { 
    
    ## make data table
    parms_df <- as.data.table(parms_df)

    ## renaming for ease
    if("(Intercept)" %in% colnames(parms_df)){
      setnames(parms_df, "(Intercept)", "Intercept")
    }
    ## select what we want and add date
    parms_df[ , date := as.Date(beta_dates)]

  } else { 
    
    ## make data table
    parms_df <- as.data.table(parms_df)
    
    ## adding in dummy intercept
    parms <- c("Intercept", parms)

    ## select what we want and add date
    parms_df[ ,..parms ][ , ":="( date = as.Date(beta_dates), Intercept = 0)]

  }

  ## return
  return(parms_df)

}


#-------------------------------------------------------------------------------
# Predictions function 
#-------------------------------------------------------------------------------
predictions <- function(parms_df, obs_df, parms, id_var, group_var, dep_var) {

  ## check
  if(!is.data.table(parms_df)){
    parms_df2 <- copy(as.data.table(parms_df))
    if(group_var %ni% colnames(parms_df2)){
      parms_df2[,(group_var) := rownames(parms_df)]
      if(class(obs_df[,get(group_var)]) == "Date"){
        parms_df2[,(group_var) := as.Date(get(group_var))]
      }
    }
  } else {
    parms_df2 <- copy(parms_df)
  }
  if(!is.data.table(obs_df)){
    obs_df2 <- copy(as.data.table(obs_df))
  } else {
    obs_df2 <- copy(obs_df)
  }
  
  ## adding in "Intercept" to our variables
  if("Intercept" %in% colnames(parms_df2) & "Intercept" %ni% parms){
    parms <- c("Intercept", parms)
  }
  if("Intercept" %in% colnames(parms_df2) & "Intercept" %ni% colnames(obs_df2)){
    obs_df2[,Intercept := 1] 
  }
  
  ## ids
  id_vars <- c(id_var, group_var)
  dep_vars <- c(dep_var, id_vars)

  ## save regular return variable
  rets <- copy(obs_df2)[,..dep_vars]

  ## it is cleaner to do the prediction multiplication with long-format dataset
  parms_df2 <- 
    data.table::melt(
      parms_df2
      ,id.vars = group_var
      ,measure.vars = parms
      ,value.name = "value.parms"
    )
  obs_df2 <- 
    data.table::melt(
      obs_df2
      ,id.vars = id_vars
      ,measure.vars = parms
      ,value.name = "value.obs"
    )
  
  ## merge together
  all_df  <- 
    merge(
      obs_df2
      ,parms_df2
      ,by = c("variable", group_var)
      ,all.x = F
    )

  ## value multiplication for predictions
  all_df <- all_df[
    ,value.pred := value.parms * value.obs ## multiply
  ][
    ,.(pred_rets = sum(value.pred, na.rm = T)) ## calculate predicted returns                                       
    ,by = setNames(
      list(                         ## by id variables (generally permno and date) 
        get(id_var)
        ,get(group_var)
      )
      ,c(
        id_var
        ,group_var
      )
    ) 
  ]

  ## merge back regular return variables
  all_df <- 
    merge(
      x = all_df
      ,y = rets
      ,by = id_vars
      ,all.x = T
    )

  return(all_df)

}


#-------------------------------------------------------------------------------
# Pctile group classification 
#-------------------------------------------------------------------------------
pctileGroups <- function(vec, pctile = 0.2) {
  return(cut(vec,
    breaks = unname(quantile(vec ,probs = seq(0, 1, by = pctile), na.rm = T)),
    include.lowest = T, labels = F
  ))
}

#-------------------------------------------------------------------------------
# Long/short portfolio returns 
#-------------------------------------------------------------------------------
longShort <- function(data, pctile_var = "tiles", group_var = "date",
  ret_var = "bh1m", weight_var = NULL
) {
  
  ## this way we do not affect the inputted data, just a local copy
  data2 <- copy(data)

  ## making min pctile of predicted rets negative for "short" 
  rets <- 
    data2[
      ,.(
        ret_var = 
          ifelse(
            !is.null(weight_var)
            ,weighted.mean(get(ret_var),get(weight_var))
            ,mean(get(ret_var))
          )
      )
      ,by = setNames(list(get(pctile_var), get(group_var)), c(pctile_var, group_var))
    ][
      get(pctile_var) == min(get(pctile_var))
      ,ret_var := -ret_var
    ]
  
  ## constructing l/s portfolio returns off values
  longshort <- 
    rets[
      get(pctile_var) == min(get(pctile_var))
      | get(pctile_var) == max(get(pctile_var))
    ][
      ,.(ret_var = sum(ret_var))
      ,by = setNames(list(get(group_var)), c(group_var))
    ]
  
  ## returning names
  setnames(longshort, "ret_var", "longshort")

  ## return
  return(longshort)

}

#-------------------------------------------------------------------------------
# Long/short portfolio returns
# Takes vectors as arguments
#-------------------------------------------------------------------------------
longShort.small <- function(
  x
  ,portfolio.id
  ,long.id = max(portfolio.id)
  ,short.id = min(portfolio.id)
  ,weight = NULL
)
{
  longs <- x[which(portfolio.id == long.id)]
  shorts <- x[which(portfolio.id == short.id)]
  if(!is.null(weight)){
    longs.weights <- weight[which(portfolio.id == long.id)]
    shorts.weights <- weight[which(portfolio.id == short.id)]
    longs.mean <- weighted.mean(x = longs, w = longs.weights) 
    shorts.mean <- weighted.mean(x = shorts, w = shorts.weights) 
  } else {
    longs.mean <- mean(x = longs) 
    shorts.mean <- mean(x = shorts) 
  }
  long_short <- longs.mean - shorts.mean
  return(long_short)
}

#-------------------------------------------------------------------------------
# Sharpe ratio summary 
#-------------------------------------------------------------------------------
sharpe <- function(data, var_name = "longshort", new_name = "sharpe")
{
  
  data <- as.data.frame(data)

  ## reassign variable name as specified
  colnames(data)[which(colnames(data) == var_name)] <- new_name
  
  Er    <- mean(data[,new_name], na.rm = TRUE)
  sigma <- sd(data[,new_name], na.rm = TRUE)
  sharpe <- Er/sigma
  
  summary <- rbind(Er, sigma, sharpe)
  rownames(summary) <- c("Expected_return", "Sigma", "Sharpe")
  colnames(summary) <- new_name

  return(summary)

}

#-------------------------------------------------------------------------------
# Rank transform dataset 
#-------------------------------------------------------------------------------
rankTransform <- function(
  data
  ,id_var
  ,group_var
  ,informat = "wide"
  ,outformat = "wide"
)
{

  data2 <- copy(data)

  if(informat == "wide"){
    data2 <- 
      melt(
        data2
        ,id.vars = c(id_var, group_var)
      )
  }
  
  data2 <-
    data2[
      ,value := frank(value,na.last = "keep",ties.method = "min")/length(value)
      ,by = .(get(group_var), variable)
    ]

  if(outformat == "wide"){
    data2 <- 
      data.table::dcast(
        data2
        ,as.formula(paste0(paste(c(id_var, group_var), collapse = "+")," ~ variable"))
        ,value.var = "value"
      )
  }

  return(data2)

}

#-------------------------------------------------------------------------------
# Checking for consistent matrix columns 
#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
# Imputing multiple matrices in a list 
#-------------------------------------------------------------------------------
seqImputeMatList <- function(mat_list, E0_init = colMeans(mat_list[[1]], na.rm = T),
  R0_init = cov(mat_list[[1]], use = "pairwise.complete.obs"), min_obs = 100,
  maxiter = 10000, tol = 1e-4, print_progress = FALSE
) {
  # checking for most likely errors ----
  if (class(mat_list) != "list") stop("`mat_list` must be a list.\n")
  sameMatCol(mat_list)
  if (any(do.call("rbind", lapply(mat_list, function(x) is.null(colnames(x)))))) {
    warning(paste0(
      "n >= 1 matrices in `mat_list` do not have column names, ",
      "so seqImputeMatList will assign dummy names to all.\n"
    ))
    mat_list <- lapply(mat_list, function(x) {
      colnames(x) <- paste0("x", 1:ncol(x))
      return(x)
    })
  }
  # initialize and loop through imputations ----
  E0 <- list(E0_init)
  R0 <- list(R0_init)
  for (i in 1:length(mat_list)) {
    raw_i <- mat_list[[i]]
    good <- names(checkMinObs(raw_i, min_obs))
    na_sort <- do.call("order", as.data.frame(-is.na(raw_i)))
    raw_i <- scale(raw_i[na_sort, good])
    if (i != 1) {
      if (any(good %ni% good_p) | any(good_p %ni% good)) {
        E0[[i]] <- colMeans(raw_i, na.rm = T)
        R0[[i]] <- cov(raw_i, use = "pairwise.complete.obs")
      } else {
        E0[[i]] <- E0[[i - 1]]
        R0[[i]] <- R0[[i - 1]]
      }
    }
    if (any(is.na(E0[[i]]))) E0[[i]][which(is.na(E0[[i]]))] <- 0
    if (any(is.na(R0[[i]]))) R0[[i]] <- diagRowCol(R0[[i]])
    if (print_progress) {
      if (!is.null(names(mat_list)) & i <= length(names(mat_list))) {
        cat(paste0(
          "Running imputations for ", names(mat_list)[i], 
          " (", i, "/", length(mat_list), ")\n"
        ))
      } else {
        cat(paste0(
          "Running imputations for list element ",
          i, "/", length(mat_list), "\n"
        ))
      }
    }
    # Actual Fortran imputation algorithm ----
    em_out <- mvn_emf(raw_i, 
      E0[[i]][good], R0[[i]][good, good], 
      maxiter = maxiter, tol = tol
    )
    E0[[i]] <- em_out$estE
    R0[[i]] <- em_out$estR
    names(E0[[i]]) <- good
    colnames(R0[[i]]) <- good
    rownames(R0[[i]]) <- good
    good_p <- good
  }
  if (!is.null(names(mat_list))) {
    names(E0) <- names(mat_list)
    names(R0) <- names(mat_list)
  }
  return(list(estE = E0, estR = R0))
}

#-------------------------------------------------------------------------------
# diagRowCol 
#-------------------------------------------------------------------------------
diagRowCol <- function(mat) {
  if (nrow(mat) != ncol(mat)) stop("Must input square matrix to diagRowCol.")
  bad <- apply(mat, FUN = function(x) length(x[is.na(x)]), MARGIN = 2)
  idx <- (1:ncol(mat))[order(-bad)]
  for(i in idx){
    if (any(is.na(mat)) == FALSE) {
      break
    } else {
      # replace with identity
      mat[ ,i] <- 0 # columns
      mat[i, ] <- 0 # rows
      mat[i,i] <- 1 # diagonal
    }
  }
  return(mat)
}

#-------------------------------------------------------------------------------
# getImpList
#-------------------------------------------------------------------------------
getImpList <- function(raw_data, E0, R0, maxiter = 100, tol = 1e-4) {
  if (length(unique(length(raw_data), length(E0), length(R0))) != 1) {
    stop("Must have covariance and means to impute with for each dataset.")
  }
  imp <- list()
  for (i in 1:length(raw_data)) {
    cols <- colnames(raw_data[[i]])[
      colnames(raw_data[[i]]) %in% names(E0[[i]]) &
      colnames(raw_data[[i]]) %in% colnames(R0[[i]])
    ]
    na_sort <- do.call("order", as.data.frame(-is.na(raw_data[[i]])))
    imp_i <- mvn_emf(raw_data[[i]][na_sort, cols], 
      E0[[i]][cols], R0[[i]][cols, cols],
      maxiter = maxiter, tol = tol
    )
    imp[[i]] <- imp_i$Ey
    colnames(imp[[i]]) <- cols
  }
  if (!is.null(names(raw_data))) {
    names(imp) <- names(raw_data)
  }
  return(imp)
}

#-------------------------------------------------------------------------------
# checkMinObs 
#-------------------------------------------------------------------------------
checkMinObs <- function(mat, min_obs) {
  return(which(apply(mat, MARGIN = 2, FUN = function(x) {
    length(x[!is.na(x)]) > min_obs
  })))
}

#-------------------------------------------------------------------------------
# Simple imputation of a vector 
#-------------------------------------------------------------------------------
imputeVec <- function(x, na.rm = T) {
  x[is.na(x)] <- mean(x, na.rm = na.rm)
  return(x)
}

#-------------------------------------------------------------------------------
# Missingness plot 
#-------------------------------------------------------------------------------
missPlot <- function(data, rhs_vars, xlab = "Stock i.d.", ylab = "Anomaly i.d.",
  title = "Missingness Map"
) {
  data2 <- copy(as.data.table(data))

  nas <- is.na(data2[,..rhs_vars])

  sort_by_na <- do.call("order", as.data.table(is.na(data2)))

  nas <- as.data.table(nas)[sort_by_na]

  col_miss <- colSums(nas)
  col_ind <- order(-col_miss) 
  nas <- nas[, ..col_ind]
  signal_id <- data.table(variable = colnames(nas), id = ncol(nas):1) 
  nas[, firm_id := 1:nrow(nas)]

  ratio <- nrow(nas)*1/(ncol(nas) - 1)
  nas_long <- merge(melt(nas, id = "firm_id"), signal_id, by = "variable") 

  miss_plot <- ggplot(data = nas_long, aes(x = firm_id, y = id)) + 
    geom_raster(aes(fill = value)) +
    labs(x = xlab, y = ylab, color = "") + 
    scale_fill_manual(
      name = "",
      labels = c("Present","Missing"),
      values = c("white", "indianred3")
    ) + 
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      legend.position = "bottom",
      legend.background = element_rect(fill = "lightgrey")
    ) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    ggtitle(title)# + 
    #coord_fixed(ratio)
  
  return(miss_plot)
}

#-------------------------------------------------------------------------------
# Sharpe table 
#-------------------------------------------------------------------------------
sharpeTable <- function(
  table_list
  ,row_groups
  ,group_len
  ,col_heads
  ,title = "Sharpe Ratios of Selected Specifications"
  ,digits = 2
)
{
  cols_align <- rep("c", length(col_heads)) 

  for(i in 1:length(table_list)){
    
    table_list[[i]] <- do.call("cbind", table_list[[i]])
    table_list[[i]] <- formatC(table_list[[i]], digits = digits, format = "f")
    
    if(ncol(table_list[[i]]) < length(col_heads)){
      for(j in 1:(length(col_heads) - ncol(table_list[[i]]))){
        table_list[[i]] <- cbind(table_list[[i]],"")
      }
    }
  
    colnames(table_list[[i]]) <- NULL

  }

  table_df <- do.call("rbind", table_list)
  
  col_widths <- (6 - 2) / ncol(table_df)
  col_widths <- paste0(round(col_widths,2), "in")
  
  out_latex <- 
    kable(
      table_df
      ,format = "latex"
      ,col.names = col_heads
      ,row.names = TRUE
      ,caption = title
      ,booktabs= TRUE
      ,align = cols_align
    ) %>%    
    column_spec(1, border_right = T, width = "2in") %>%
    column_spec(2:(ncol(table_df)+1), width = col_widths)


  for(i in 1:length(row_groups)){
    out_latex <- 
      group_rows(
        out_latex
        ,row_groups[i]
        ,1 + (i-1)*group_len
        ,i*group_len
      )
  }

  return(out_latex)
}

#-------------------------------------------------------------------------------
# Principal components custom 
# Returns principal components weighted data
# Arguments:
#   data: the data to get PCs for
#   covmat: covariance matrix of data. 
#     If blank, will calculate covariance of data from data argument
#   scale: if we want to scale the data first
#   n.components: number of components to use
#-------------------------------------------------------------------------------
prcompCustom <- function(data, covmat = cov(data, use = "pairwise.complete.obs"),
  scale = TRUE, n.components = max(ncol(data))
) {

  if(scale){
    data <- scale(data, center = T, scale = T)
  }

  ## if we put in more components than data 
  ## i.e. if we had to eliminate missing variable columns
  if (n.components > ncol(data)) {
    n.components <- ncol(data)
  }

  eigens <- eigen(covmat)

  eigenvectors <- eigens$vectors[,1:n.components]
  rownames(eigenvectors) <- colnames(data)
  colnames(eigenvectors) <- paste0("PC",1:n.components)

  eigenvalues <- eigens$values[1:n.components]
  names(eigenvalues) <- paste0("PC",1:n.components)

  prcomps <- data %*% eigenvectors
  colnames(prcomps) <- paste0("PC",1:n.components)

  out <- list(
    eigenvalues = eigenvalues,
    rotation = eigenvectors,
    transformed = prcomps
  )

  return(out)

}

#-------------------------------------------------------------------------------
# Winsorizing data
#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
# Function: getBetas
#  Takes covariance matrix and empty dataset
#  Calculates EMF betas for each observation
#  applying to each row of observed data
#  we select the missing and not-missing columns
#  we then calculate betas
#  using the partitioned covariance matrix
#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
# Function: stopQuietly
#  Stops R without error message
#-------------------------------------------------------------------------------
stopQuietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

#-------------------------------------------------------------------------------
# Function: wtdMean
#  Enables weighed mean with na.rm for `x` and `w`
#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
# Function: na2zero
#  Replaces NAs in vector with 0
#-------------------------------------------------------------------------------
na2zero <- function(x) {
  return(ifelse(is.na(x), 0, x))
}

#-------------------------------------------------------------------------------
# Function: lsos
#  lists objects by size
# I believe this comes from jwutil package
# (a stackoverflow user)
#-------------------------------------------------------------------------------

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


# Skew of a single variable. Not used right now
skew <- function(x, na.rm = FALSE) {
  if (na.rm) x <- x[!is.na(x)]
  n <- length(x)
  return((sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2))
}

# Read in a specific file from a zip using function of your choice
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

# Convert stata YYYYmM format to zoo's yearmon format
stataMon2yearmon <- function(statamon) {
  if (!all(statamon %like% "m")) { 
    stop("`statamon` must be Stata-formatted year-month, i.e. '1996m1', '2019m12'\n")
  }
  zoo::as.yearmon(sapply(base::strsplit(statamon, "m"), function(x) {
    as.numeric(x[1]) + (as.numeric(x[2]) - 1) / 12 
  }))
}

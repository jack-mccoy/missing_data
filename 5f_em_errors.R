# 2023 10 Andrew: made from ppca-emp-crossval.R
# find imputation errors by masking 1/num_folds of the observed data
# and running MVN EM

# Setup -------------------------------------------------------------------
rm(list = ls())
library(data.table)
library(optparse)
library(pcaMethods)
library(zoo)
library(tidyverse)
library(doParallel)

source('functions.R')

# file paths
getFilePaths()
dir.create(paste0(FILEPATHS$out_path, 'plots'), showWarnings = FALSE)

# settings ----------------------------------------------------------------
maxiter <- 100000
tol = 1e-4
num_folds = 10
yearm_list <- c('Jun 1990','Jun 2000','June 2010')
seed <- 864
ncores = round(0.75*detectCores())
fastrun = FALSE # TRUE for testing

# read data ----------------------------------------------------------------
signals <- fread(paste0(FILEPATHS$data_path, 'bcsignals/bcsignals_none.csv'))
signals[, yyyymm := as.yearmon(yyyymm)] # careful with reading yearmon format from csv!
signals_good <- names(signals %>% select(-c("permno", "yyyymm")))

crsp = fread(paste0(FILEPATHS$data_path, 'raw/crsp_data.csv'))
crsp <- crsp[, yyyymm := as.yearmon(yyyymm)]
# declare function ---------------------------------------------------------

# function for finding prediction error for a fold
fold_err = function(itest){
    # example: itest = ifold[[1]]
    itrain = setdiff(iobs, itest)

    # force test ids to missing
    train_data <- xsignal
    train_data[itest] <- NA
    test_data <- xsignal
    test_data[itrain] <- NA

    # remove stocks with no obs
    badstock <- apply(is.na(train_data), 1, all)
    train_data <- train_data[!badstock, ]
    test_data <- test_data[!badstock, ]

    # remove signals with no obs
    badsignal <- apply(is.na(train_data), 2, all)
    train_data <- train_data[, !badsignal]
    test_data <- test_data[, !badsignal]

    # prepare for EM, sort by missingness for speed
    na_sort = do.call("order", as.data.frame(-is.na(train_data))) 
    train_data = train_data[na_sort, ] 
    test_data = test_data[na_sort, ]

    # run EM
    tic = Sys.time()
    E0 = 0*colMeans(train_data, na.rm = TRUE)
    R0 = diag(diag(cov(train_data, use = "pairwise.complete.obs")))     
    em_out = mvn_emf(train_data, E0, R0, maxiter = maxiter, tol = tol, 
        update_estE = FALSE)
    test_pred = em_out$Ey
    rownames(test_pred) = rownames(test_data)
    colnames(test_pred) = colnames(test_data)
    toc = Sys.time()
    difftime(toc, tic, units = "mins")

    # output error for test data in long format
    errlong <- test_data - test_pred %>% as.data.frame() 
    errlong <- errlong %>%
        mutate(permno = as.integer(row.names(test_data))) %>%
        pivot_longer(cols = -permno, names_to = 'signalname', values_to = "err")   %>% 
        filter(!is.na(err))

    return(errlong)

} # end function fold_err

# find errors for each year ------------------------------------------------

# loop over years
xerrall = foreach (i=1:length(yearm_list), .combine = 'rbind') %do% {
    # testing: yearm_cur = yearm_list[1]
    yearm_cur = yearm_list[i]

    print(paste0("EM errors for ", yearm_cur))

    # make signals matrix (no identifiers, it's a nice matrix)
    xsignal <- signals[yyyymm == yearm_cur] %>% 
        select(-c(yyyymm,permno))  %>% 
        as.matrix()
    rownames(xsignal) <- signals[yyyymm == yearm_cur]$permno

    # fast run for checking
    if (fastrun == TRUE){
        colselect = floor(seq(1, ncol(xsignal), length.out = 20))
        xsignal = xsignal[ , colselect]
    }

    # shuffle observed stock-signal indexes
    iobs <- which(!is.na(xsignal))
    set.seed(seed)
    iobs <- sample(iobs, length(iobs))

    # split into ifolds list
    ifold <- split(iobs, ceiling(seq_along(iobs) / length(iobs) * num_folds))

    cl = makePSOCKcluster(ncores)
    registerDoParallel(cl)

    print('finding EM errors for a whole cross section')
    print(paste0('time started = ', Sys.time()))
    tic = Sys.time()
    xerrlong = foreach(
        i = 1:num_folds,
        .combine = rbind,
        .packages = c("tidyverse","data.table")
    ) %dopar% {
        print(paste0("fold = ", i, " of ", num_folds))
        errfold = fold_err(ifold[[i]])
        return(errfold)
    }
    toc = Sys.time()
    print(difftime(toc, tic, units = "mins"))

    xerrlong$yyyymm = yearm_cur

    return(xerrlong)

} # end for yearm_cur

# add market equity
xerrall2 = xerrall  %>% 
    mutate(yyyymm = as.yearmon(yyyymm))  %>%
    filter(!is.na(err))  %>% 
    left_join(crsp[ , .(permno, yyyymm, me)], by = c('permno','yyyymm')) 

# save just in case -------------------------------------------------------
fwrite(xerrall2, paste0(FILEPATHS$out_path, 'temp_EM_errors.csv'))


# make plots ----------------------------------------------------------------
num_size_bins = 10

# add market equity bins
plotme = xerrall2  %>% 
    filter(!is.na(me))  %>% 
    group_by(yyyymm)  %>%
    mutate(bin_me = ntile(me, num_size_bins)
    , bin_me = as.factor(bin_me)
    , yyyymm = as.factor(yyyymm))  %>% 
    arrange(yyyymm, bin_me)   %>% 
    ungroup()

# plot rmse by group
tempp = plotme  %>% 
    group_by(yyyymm, bin_me)  %>% 
    summarize(rmse = sqrt(mean(err^2)))  %>%
    ggplot(aes(x = bin_me, y = rmse, group = yyyymm)) +
    geom_hline(yintercept = 1.0, color = 'black') +
    geom_hline(yintercept = 0.0, color = 'black') +    
    geom_line(aes(color = yyyymm, linetype = yyyymm)) +
    geom_point(aes(color = yyyymm, shape = yyyymm)) +
    theme_bw() +
    theme(
            legend.position = c(8,8)/10,
            legend.background = element_blank(),
            legend.key = element_blank(),
            legend.spacing.y = unit(0.005, 'cm'),
            legend.spacing.x = unit(0.1, 'cm'),
            legend.box = 'horizontal',
            legend.key.width = unit(0.55, 'cm'),
            legend.key.height = unit(0.38, 'cm'),
            legend.text = element_text(size = 7),
            legend.title = element_blank()
    ) + 
    scale_linetype_manual(values = c('solid', 'longdash', 'dotted')) +
    scale_size_manual(values = c(0.8, 0.8, 0.5)) +
    scale_color_manual(values = c(MATRED, MATBLUE, MATYELLOW)) +
    xlab('Market Equity Decile') +
    ylab('Out-of-Sample RMSE') +
    scale_y_continuous(limits = c(0.5, 1.0))

ggsave(tempp, 
    filename = paste0(FILEPATHS$out_path, 'plots/EM_error_rmse.pdf'),
    width = 8, height = 5, unit = "in", scale = 0.55
)


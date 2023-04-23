# Imputes data following Bryzgalova, Lettau, Lerner, Pelger 2023
# The Backward-cross-section (B-XS) model, based on Appendix A

# email from Markus Pelger April 17 2023
# "All empirical results including the IPCA results are based on the
# implementation described in Appendix A."

# input: bcsignals_none.csv
# output: bcsignals_bllp.csv

# AC: Eq A.1 in BLLP is equivalent to a complete case regression of
# signals on loadings across signalnames (by permno)
# Seems too easy so I kept the sanity check in comments below

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
  optparse::make_option(c("--num_PCs"),
                        type = "numeric", default = 6,
                        help = "number of PCs of AC Cov to use"),
  optparse::make_option(c("--out_path"), # Need this flexibility. Files are too big to store in non-scratch
                        type = "character", default = "../output/",
                        help = "directory including input and output files"),
  optparse::make_option(c("--out_name"), 
                        type = "character", default = "bcsignals_bllp.csv",
                        help = "directory including input and output files"),  
  optparse::make_option(c("--impute_vec"),
                        type = "character", default = "bm,mom6m",
                        help = "a comma-separated list of values or .txt file to scan"),
  optparse::make_option(c("--cores_frac"),
                        type = "numeric", default = 1.0,
                        help = "fraction of total cores to use")    
)


opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)


# ## debug -------------------------------------------------------------------
# opt$impute_vec = '../data/signals_best125_1985.txt'


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
bcsignals = fread(paste0(opt$out_path, "bcsignals/bcsignals_none.csv")) 
bcsignals <- bcsignals[, .SD, .SDcols = c("permno", "yyyymm", opt$impute_vec)]
bcsignals[ , yyyymm := as.yearmon(yyyymm)]  # careful with reading yearmon format from csv!


#==============================================================================#
# impute_one_month function ----
#==============================================================================#




impute_one_month = function(yearm_cur){
  

  # set up ------------------------------------------------------------------

  
  # current signal matrix
  signalmat = bcsignals[yyyymm == yearm_cur] %>% select(-yyyymm) %>% select(-c(permno)) %>% 
    as.matrix()
  rownames(signalmat) = bcsignals[yyyymm == yearm_cur]$permno
  
  cov_ac = cov(signalmat, use = 'pairwise.complete.obs')
  
    # # debug
    # ibad = which(!complete.cases(cov_ac))[2]
    # rownames(cov_ac)[ibad]
    # which(is.na(cov_ac[ibad,]))
    # 
    # bcsignals[yyyymm == yearm_cur, .(permno,yyyymm,betatailrisk,firmagemom,ageipo)] %>% 
    #   print(topn = 20)
    # 
    # bcsignals[yyyymm == yearm_cur, .(permno,yyyymm,betatailrisk,firmagemom,ageipo)] %>% 
    #   summarize(cov(betatailrisk, firmagemom))
    
  
  
  # find preliminary loadings Lambda (num_signals x num_PCs)
  Lamtil = eigen(cov_ac)$vectors[ , 1:opt$num_PCs]
  rownames(Lamtil) = colnames(signalmat)
  colnames(Lamtil) = paste0('load', 1:opt$num_PCs)
  
  
  # make long data with permno, name, value, lag, loadings
  temp1 = bcsignals[yyyymm == yearm_cur] %>% 
    select(-yyyymm) %>% 
    pivot_longer(cols = !c('permno')) %>% 
    setDT() 
  
  temp2 = bcsignals[yyyymm == yearm_cur-1/12] %>% 
    select(-yyyymm) %>% 
    pivot_longer(cols = !c('permno')) %>% 
    rename(lag = value) %>% 
    setDT()
  
  longdat = temp1 %>% 
    left_join(temp2, by = c('permno','name')) %>% 
    rename(signalname = name) %>% 
    left_join(
      data.table(
        signalname = rownames(Lamtil), Lamtil
      )    
      , by = 'signalname'
    ) 
  
  
  # find factors ---------------------------------------------------
  
  # https://stackoverflow.com/questions/29420674/fit-model-by-group-using-data-table-package
  
  # run complete case regression of signals on factor loadings by permno
  form = paste0("value ~ 0 + ", paste(paste0("load", 1:opt$num_PCs), collapse = "+")) %>% as.formula
  formulas <- list(form)
  temp1 = longdat[
    , lapply(formulas, function(x) list(coef(lm(x, data=.SD))))
    , by=permno
  ]
  
  # temp is data frame of lists, we need to unpack it and clean up
  temp2 = apply(
    temp1, 1, function(x) {unlist(x)}
  ) %>% t() %>% 
    as.data.table()
  names(temp2) = sub('^V1.load', 'fac', names(temp2))
  
  stockfac = temp2
  
  # # Sanity check  
  # # this check uses notation from BLLP
  # stocki = 5
  # 
  # L = ncol(signalmat)
  # K = opt$num_PCs
  # 
  # # make matrix calc
  # signalzero = signalmat
  # signalzero[which(is.na(signalmat))] = 0
  # 
  # obs = matrix(0, nrow(signalmat), ncol(signalmat))
  # obs[which(!is.na(signalmat))] = 1
  # 
  # 
  # # sum version
  # x = matrix(0,K,K)
  # for (l in 1:L){
  #   x = x +  obs[stocki,l] * matrix(Lamtil[l, ], nrow = K) %*% matrix(Lamtil[l, ], nrow = 1)
  # }
  # x = solve(x)
  # 
  # y = matrix(0,1,K)
  # for (l in 1:L){
  #   y = y + obs[stocki,l] * matrix(Lamtil[l,], nrow = 1) * signalzero[stocki, l]
  # }
  # 
  # Fhatcheck = x %*% t(y) %>% as.numeric
  # 
  # # Compare
  # Fhatcheck
  # Fhatlong[permno == rownames(signalzero)[stocki] ]
  
  
  
  
  # stage 2: add lagged info ------------------------------------------------
  
  # BLLP: B-XS model (no eqn number)
  
  temp = longdat[ , .(permno,signalname,value,lag)] %>% 
    left_join(stockfac, by = 'permno')
  
  # run (complete case) regression of signals on factors and lagged signal by signalname
  form = paste0("value ~ 0 + lag + ", paste(paste0("fac", 1:opt$num_PCs), collapse = "+")) %>% as.formula
  formulas <- list(form)
  temp2 = temp[
    , lapply(formulas, function(x) list(coef(lm(x, data=.SD))))
    , by=signalname
  ]
  
  # clean up, create dt of ts model coefficients
  tscoef = apply(
    temp2, 1, function(x) {unlist(x)}
  ) %>% t() %>% 
    as.data.table()
  names(tscoef) = sub('^V1.', '', names(tscoef))
  tscoef = tscoef %>% 
    mutate(across(!signalname, as.numeric))
  
  # make predictions valuehat
  longfac = copy(longdat)[
    , .(permno, signalname, lag)
  ] %>% 
    left_join(stockfac, by = 'permno') %>% 
    mutate(
      lag = ifelse(is.na(lag), 0, lag) # if lag is missing, assume zero, not sure how to do this
    )
  
  longcoef = copy(longdat)[
    , .(permno, signalname)
  ] %>% 
    left_join(tscoef, by = 'signalname')
  
  temp1 = (longfac %>% select(-c(permno,signalname)) %>% as.matrix) *
    (longcoef %>% select(-c(permno,signalname)) %>% as.matrix)
  
  longdat$valuehat = rowSums(temp1)
  
  # make final data valueimp
  # replace with NA if missing (not sure how to handle this)
  # can come up if a stock has only one observed signal
  longdat = longdat %>% 
    mutate(
      valueimp = case_when(
        !is.na(value) ~ value
        , is.na(value) & !is.na(valuehat) ~ valuehat
        , T ~ 0
      )
    )
  
  
  impcur = longdat %>% 
    select(permno, signalname, valueimp) %>% 
    pivot_wider(id_cols = permno, names_from = signalname, values_from = valueimp) %>% 
    mutate(yyyymm = yearm_cur) %>% 
    select(permno, yyyymm, everything()) %>% 
    setDT()
  
  return(impcur)
  
  
} # end imput_one_month





#==============================================================================#
# loop over months ----
#==============================================================================#


# hard code
yearm_min = as.yearmon('1985-01')
yearm_max = as.yearmon('2020-12')


ymlist = bcsignals$yyyymm %>% unique()
ymlist = ymlist[
  ymlist >= yearm_min & ymlist <= yearm_max
]


tic = Sys.time()
bcsignals_bllp = foreach(yearm_cur = ymlist, .combine = rbind) %do% {
  print(paste0('bllp imputing yearm = ', yearm_cur))
  impcur = impute_one_month(yearm_cur)
  toc = Sys.time()
  print(toc - tic)
  
  return(impcur)
}



#==============================================================================#
# write to disk ----
#==============================================================================#

fwrite(bcsignals_bllp, paste0(opt$out_path, 'bcsignals/', opt$out_name))


#==============================================================================#
# check ----
#==============================================================================#


names(bcsignals_bllp) 

signallist = names(bcsignals_bllp)
signallist = signallist[3:length(signallist)]


bcsignals_bllp[
  , lapply(.SD, sd)
  , .SDcols = signallist
]



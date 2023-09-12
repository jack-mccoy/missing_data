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

# 12 seconds per month, about 90 minutes for the full data

#==============================================================================#
# Packages ----
#==============================================================================#

library(car) # Box-Cox stuff
library(data.table) # Standard data manipulation
library(optparse) # Python-like option parsing
library(RPostgres) # SQL query to WRDS
library(zoo) # yearmon convention is nice to work with here
library(tidyverse) # sorry Jack

source("functions.R")

#==============================================================================#
# Option parsing ----
#==============================================================================#

# check if on cluster
on_cluster = Sys.getenv('SGE_TASK_ID') != ""

option_list <- list(
  optparse::make_option(c("--num_PCs"),
                        type = "numeric", default = 6,
                        help = "number of PCs of AC Cov to use"),
  #optparse::make_option(c("--out_path"), # Need this flexibility. Files are too big to store in non-scratch
  #                      type = "character", default = "../output/",
  #                      help = "directory including input and output files"),
  #optparse::make_option(c("--out_name"), 
  #                      type = "character", default = "bcsignals_bllp.csv",
  #                      help = "directory including input and output files"),  
  #optparse::make_option(c("--impute_vec"),
  #                      type = "character", default = "../data/signals_best125_1985.txt",
  #                      help = "a comma-separated list of values or .txt file to scan"),
  optparse::make_option(c("--lag_max_years"),
                        type = "numeric", default = 2,
                        help = "max num years to look back for 2nd stage")
)


opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Load file paths
getFilePaths()

# Get the anomalies as a nice vector
impute_vec <- unpackSignalList(FILEPATHS$signal_list)

if (length(impute_vec) <= opt$num_PCs){
  stop('error num signals < num pcs ')
}

#==============================================================================#
# Setup ----
#==============================================================================#

# Months to impute
yrmon_list <- zoo::as.yearmon(paste0(month.abb, " ", opt$impute_yr))

closeAllConnections()

# Read in the bc-transformed signals 
bcsignals = fread(paste0(FILEPATHS$data_path, "bcsignals/bcsignals_none.csv")) 
bcsignals <- bcsignals[, .SD, .SDcols = c("permno", "yyyymm", impute_vec)]
bcsignals[ , yyyymm := as.yearmon(yyyymm)]  # careful with reading yearmon format from csv!

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



#==============================================================================#
# impute_one_month function ----
#==============================================================================#




impute_one_month = function(yearm_cur){
  

  # set up real-time sample ------------------------------------------------------------------
  
  # crop months
  longdat = bcsignals[yyyymm <= yearm_cur & yyyymm >= yearm_cur - opt$lag_max_years]%>% 
    pivot_longer(cols = !c('permno','yyyymm'), names_to = 'signalname') %>%
    arrange(permno,signalname,-yyyymm)
  
  # keep only permnos with some data this month 
  permnolist = longdat %>% 
    filter(yyyymm == yearm_cur & !is.na(value))  %>% 
    pull(permno)
  longdat = longdat %>% filter(permno %in% permnolist)
  
  # keep only relevant rebalancing periods
  longdat = longdat %>% 
    left_join(signaldoc, by = 'signalname') %>% mutate(
      mon = round(12*(yyyymm - floor(yyyymm))+1)
    ) %>% 
    filter(
      port_period == 1 | (
        port_period == 3 & mon %in%  (month(yearm_cur) + 3*(-12:12))
      ) | (
        port_period == 6 & mon %in%  (month(yearm_cur) + 6*(-12:12))
      ) | (
        port_period == 12 & mon %in%  (month(yearm_cur) + 12*(-12:12))
      )
    )
  
  # drop lags if missing (but keep all current rows, even if missing, this is important)
  longdat = longdat %>% mutate(current_val = yyyymm==yearm_cur)  %>% 
    filter(current_val | (!current_val & !is.na(value))) %>% 
    arrange(permno,signalname,-yyyymm)
  
  # keep only the most recent 1 lag (hardcode 1 for now)
  setDT(longdat) # group_by is slow
  longdat[ , lagnum := seq(1,.N), by = c('permno','signalname','current_val')]
  longdat = longdat[current_val | lagnum == 1]
  longdat = longdat %>% 
    pivot_wider(id_cols = c('permno','signalname'), names_from = 'current_val'
                , values_from = value, names_prefix = 'current_') %>% 
    transmute(permno,signalname, value = current_TRUE, lag = current_FALSE)
  
  # checks (see how often lags are observed if current is missing)
  # might use this to look back further than 2 years
  # longdat %>% mutate(obs = !is.na(value)) %>% group_by(obs) %>% summarize(frac_obs_lag = mean(!is.na(lag)))

  # find prelim loadings ----------------------------------------------------

  # current signal matrix
  signalmat = longdat %>% select(permno, signalname, value) %>%   
    pivot_wider(id_cols = 'permno', names_from = signalname, values_from = value) %>%
    select(-permno) %>% 
    as.matrix()

  cov_ac = cov(signalmat, use = 'pairwise.complete.obs')
  
  # error handling 
  # if no overlapping obs, replace with zero
  if (any(is.na(cov_ac))){
    cov_ac[is.na(cov_ac)] = 0
  }
  
  # find preliminary loadings Lambda (num_signals x num_PCs)
  Lamtil = eigen(cov_ac)$vectors[ , 1:opt$num_PCs]
  rownames(Lamtil) = colnames(signalmat)
  colnames(Lamtil) = paste0('load', 1:opt$num_PCs)
  Lamtil = data.table(signalname = rownames(Lamtil), Lamtil)    
  
  # merge onto longdat
  longdat = longdat %>% 
    left_join(Lamtil, by = 'signalname')
  
  # find factors ---------------------------------------------------
  
  # https://stackoverflow.com/questions/29420674/fit-model-by-group-using-data-table-package
  
  # run complete case regression of signals on factor loadings by permno
  form = paste0("value ~ 0 + ", paste(paste0("load", 1:opt$num_PCs), collapse = "+")) %>% as.formula
  formulas <- list(form)
  setDT(longdat)
  temp1 = longdat[
    !is.na(value)
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
  
  # merge on factor data
  #   drops permnos if no factor data (less than two signals observed)
  temp = longdat[ , .(permno,signalname,value,lag)] %>% 
    inner_join(stockfac, by = 'permno') 
  
  # run (complete case) regression of signals on factors and lagged signal by signalname
  form = paste0("value ~ 0 + lag + ", paste(paste0("fac", 1:opt$num_PCs), collapse = "+")) %>% as.formula
  formulas <- list(form)
  temp2 = temp[
    !is.na(value) & !is.na(lag)
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
      lag = ifelse(is.na(lag), 0, lag) # if lag is missing, assume lag is zero (not sure how handle this)
    )
  
  longcoef = copy(longdat)[
    , .(permno, signalname)
  ] %>% 
    left_join(tscoef, by = 'signalname')
  
  # multiply factors and coefficients
  temp1 = (longfac %>% select(-c(permno,signalname)) %>% as.matrix) *
    (longcoef %>% select(-c(permno,signalname)) %>% as.matrix)
  
  # sum
  longdat$valuehat = rowSums(temp1)
  
  # check
  # ggplot(longdat[1:1000, ], aes(x=valuehat,y=value)) +geom_point() + geom_abline(slope = 1)
  
  # make final data valueimp
  # replace with 0 if value and valuehat are missing (not sure how to handle this)
  # can come up if a stock has only one observed signal
  # (occurs about 10 percent of the time)
  longdat = longdat %>% 
    mutate(
      valueimp = case_when(
        !is.na(value) ~ value
        , is.na(value) & !is.na(valuehat) ~ valuehat
        , T ~ 0
      )
    )

  # make wide
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
yearm_max = as.yearmon('2020-11')

ymlist = bcsignals$yyyymm %>% unique()
ymlist = ymlist[
  ymlist >= yearm_min & ymlist <= yearm_max
]


tic = Sys.time()
bcsignals_bllp = foreach(yearm_cur = ymlist, .combine = rbind) %do% {
  
  print(paste0('bllp imputing yearm = ', yearm_cur))
  
  impcur = impute_one_month(yearm_cur)
  toc = Sys.time()
  
  min_per_month = as.numeric(toc-tic, units = 'mins') / which(yearm_cur == ymlist) 
  min_remain = (length(ymlist) - which(yearm_cur == ymlist) ) * min_per_month
  print(paste0('min remain = ', min_remain))
  
  return(impcur)
}

#==============================================================================#
# write to disk ----
#==============================================================================#

fwrite(bcsignals_bllp, 
    paste0(FILEPATHS$data_path, 'bcsignals/bcsignals_bllp', opt$num_PCs, ".csv"))

#==============================================================================#
# check ----
#==============================================================================#


names(bcsignals_bllp) 

signallist = names(bcsignals_bllp)
signallist = signallist[3:length(signallist)]


tempfun = function(x) sum(!is.na(x))
tempsum = bcsignals_bllp[
  , by = yyyymm
  , lapply(.SD, tempfun)
  , .SDcols = signallist
] 

temptab = tempsum %>% 
  mutate(year = as.numeric(floor(yyyymm))) %>% 
  group_by(year) %>% 
  summarize(across(all_of(signallist), mean)) %>% 
  arrange(year) %>% 
  pivot_longer(cols = -year)


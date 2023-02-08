# 2023 02 trying ml stuff on chen-Zimmermann data

#==============================================================================#
# Setup ----
#==============================================================================#
rm(list = ls())
library(ranger)
library(data.table)
library(tidyverse)
library(ggplot2)
library(zoo)
library(janitor)

opt = list(
  data_file = '../output/bcsignals/bcsignals_none.csv'
  , yearm_begin = '1995-06'
  , refit_period = 12 # months
  , train_months = 12*10 # months
  , yearm_end = '2020-12'
  , cores_frac = 0.5
)

#==============================================================================#
# Data pull ----
#==============================================================================#

# Read in data 
dat <- fread(opt$data_file)
setnames(dat, colnames(dat), tolower(colnames(dat))) # ensure lower
dat[, yyyymm := as.yearmon(yyyymm)]
crsp_data <- fread("../data/crsp_data.csv")[, .(permno, yyyymm, ret)] # me is in dat
crsp_data[ , yyyymm := as.yearmon(yyyymm)]

# adjust timing and merge
crsp_data[ , yyyymm := yyyymm - 1/12] # lead returns
setnames(crsp_data, 'ret', 'bh1m')
dat <- merge(dat, crsp_data, by = c("permno", "yyyymm"))

# Memory
rm(crsp_data)

# Simple mean imputations. Shouldn't matter for EM-imputed data
# Easiest to just catch it all here for non-imputed data
signals_list = names(dat) %>% setdiff(c('permno','yyyymm','sic3','me'))

imputeVec <- function(x, na.rm = T) {
  x[is.na(x)] <- mean(x, na.rm = na.rm)
  return(x)
}
dat[,
        (signals_list) := lapply(.SD, imputeVec, na.rm = T),
        .SDcols = signals_list,
        by = .(yyyymm)
]

# keep only complete cases
dat <- dat[complete.cases(
  dat[, .SD, .SDcols = c("permno", "yyyymm", "bh1m", "me", signals_list)]
)]


#==============================================================================# 
# Loop over forecast periods ----
#==============================================================================#

# create list of forecast periods
forecast_list = data.frame(
  train_end = seq(as.yearmon(opt$yearm_begin), as.yearmon(opt$yearm_end), opt$refit_period/12)
) %>% 
  mutate(
    train_begin = train_end - opt$train_months/12
    , test_begin = train_end + 1/12
    , test_end = lead(train_end) - 1/12
  )

# loop
outlist = list()
for (fi in 1:dim(forecast_list)[1]){
  fcur = forecast_list[fi, ]
  print(paste0('forecasting in ', fcur$train_end))
  
  # fit = ranger(form, data=train, num.trees = 20,  num.threads = floor(parallel::detectCores()*opt$cores_frac))
  # pred.ra = predict(ra$fit, test) # annoying middle step for ranger
  # Ebh1m = pred.ra$predictions
  
  fit = lm(form, data=train)
  Ebh1m = predict(fit, test)
  
  outlist[[fi]] = test[,.(permno,yyyymm)][
    , Ebh1m := Ebh1m
  ]
}

# bind output
forecast = rbindlist(outlist)


#==============================================================================# 
# Portfolios ----
#==============================================================================#


temp = merge(dat[, .(permno,yyyymm,bh1m)], forecast, all.x=T, by = c('permno','yyyymm') )
temp %>% 
  group_by(yyyymm) %>% 
  mutate(port = ntile(Ebh1m,10)) %>% 
  group_by(port) %>% 
  summarize(
    rbar = mean(bh1m), vol = sd(bh1m)
  ) 


port = temp %>% 
  group_by(yyyymm) %>% 
  mutate(port = ntile(Ebh1m,10)) %>% 
  group_by(port,yyyymm) %>% 
  summarize(
    bh1m = mean(bh1m, na.rm=T)
  ) 

port %>% filter(port %in% c(1,10))   %>% 
  pivot_wider(id_cols = 'yyyymm', names_from = 'port', values_from = 'bh1m', names_prefix = 'port') %>% 
  mutate(portLS = port10 - port1) %>% 
  filter(!is.na(portLS)) %>% 
  pivot_longer(
    cols = -yyyymm
  ) %>% 
  group_by(name) %>% 
  summarize(
    rbar = mean(value), vol = sd(value), SRann = rbar/vol*sqrt(12)
  )








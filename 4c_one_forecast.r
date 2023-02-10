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

library(tensorflow)
library(keras)

opt = list(
  data_file = '../output/bcsignals/bcsignals_none.csv'
  , yearm_begin = '1995-06'
  , refit_period = 12 # months
  , insamp_months = 12*10 # months
  , yearm_end = '2020-12'
  , cores_frac = 0.5
  , model = 'keras' # 'ranger' or 'lm' or 'keras'
)

# output names
customname = NULL
outroot = '../output/forecast/'

# settings for specific algos
# max.depth = NULL dramatically slows things down
optranger = list( 
  num.tree = 500, max.depth = 3, mtry = 30
)

optkeras = list(
  lambda1 = 1e-3, epochs = 100, batch_size = 10000, learning_rate = 0.001, patience = 5
)

# keep this subset for testing, use NULL for keep all
signals_keep = NULL # e.g. c('size','bm','mom12m')

# auto setup
if (is.null(customname)) {
  outfolder = paste(
    str_remove(basename(opt$data_file), '.csv')
    , opt$model, sep = '-'
  )
} else {
  outfolder = customname
}
dir.create(paste0(outroot), showWarnings = F)
dir.create(paste0(outroot,outfolder), showWarnings = F)


#==============================================================================#
# Data pull ----
#==============================================================================#

# Read in data 
dat <- fread(opt$data_file)
setnames(dat, colnames(dat), tolower(colnames(dat))) # ensure lower
dat[, yyyymm := as.yearmon(yyyymm)]

# bcsignals_none.csv seems to have crsp info stuff that we should remove
dat[ , ':=' (me = NULL)]

crsp_data <- fread("../data/crsp_data.csv")[, .(permno, yyyymm, me, ret)] 
crsp_data[ , yyyymm := as.yearmon(yyyymm)]

# split crsp into known and unknown
crsp_info = crsp_data %>% select(permno, yyyymm, me)
crsp_bh1m = crsp_data %>% select(permno, yyyymm, ret) %>% 
  mutate(yyyymm := yyyymm - 1/12) %>% 
  rename(bh1m = ret)

dat <- merge(dat, crsp_info, by = c("permno", "yyyymm"), all.x = T)
dat <- merge(crsp_bh1m, dat, by = c("permno", "yyyymm"), all.x = T)

# Memory
rm(crsp_data, crsp_info, crsp_bh1m)

# Simple mean imputations. Shouldn't matter for EM-imputed data
# Easiest to just catch it all here for non-imputed data
# should be a cleaner way to do this
signals_list = names(dat) %>% setdiff(c('permno','yyyymm','sic3','me','bh1m'))

# subset to select ones for debugging
if (!is.null(signals_keep)){
  signals_list = intersect(signals_keep, signals_list)
}

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
# fit_and_forecast function ----
#==============================================================================#

# function for choosing forecast
fit_and_forecast = function(modelname){
  if (modelname == 'lm'){
    form = paste0('bh1m ~ ', paste(signals_list, collapse = '+'))  
    fit = lm(form, data=insamp)
    Ebh1m = predict(fit, oos)
    
  } else if (modelname == 'ranger'){
    
    fit = ranger(x = as.matrix(insamp %>% select(all_of(signals_list)))
                 , y = as.matrix(insamp %>% select(bh1m))
                 , num.trees = optranger$num.tree
                 , num.threads = floor(parallel::detectCores()*opt$cores_frac)
                 , max.depth = optranger$max.depth
                 , mtry = optranger$mtry
    )
    
    pred.ra = predict(fit, oos) 
    Ebh1m = pred.ra$predictions
    
  } else if (modelname == 'keras') {
    
    # create norm layer, fit the state of layer to data
    #   seems to be necessary, or something like it that uses adapt
    normalizer <- layer_normalization(axis = -1L)
    normalizer %>% adapt(as.matrix(insamp %>% select(all_of(signals_list))))
    
    # build model
    keras_model = keras_model_sequential() %>% 
      normalizer() %>% 
      layer_dense(32, activation = 'relu'
                  , kernel_regularizer = regularizer_l1(optkeras$lambda1)) %>%
      # layer_dense(64, activation = 'relu') %>%
      layer_dense(1)
    
    keras_model %>%  compile(
      loss = 'mean_squared_error',
      optimizer = optimizer_adam(
        learning_rate = optkeras$learning_rate
      )
    )
    
    # fit
    # ref: https://tensorflow.rstudio.com/guides/keras/basics.html#callbacks
    history <- keras_model %>% fit(
      as.matrix(insamp %>% select(all_of(signals_list))),
      as.matrix(insamp %>% select(bh1m)),
      verbose = 1,
      epochs = optkeras$epochs,
      batch_size = optkeras$batch_size,
      callbacks = callback_early_stopping(patience = optkeras$patience)
    )
    
    Ebh1m = predict(keras_model, as.matrix(oos %>% select(all_of(signals_list))))
    
  }
  
  
  return(Ebh1m)
}  



#==============================================================================# 
# Loop over sample periods ----
#==============================================================================#

# create list of forecast periods
sample_list = data.frame(
  insamp_end = seq(as.yearmon(opt$yearm_begin), as.yearmon(opt$yearm_end), opt$refit_period/12)
) %>% 
  mutate(
    insamp_begin = insamp_end - opt$insamp_months/12
    , oos_begin = insamp_end + 1/12
    , oos_end = lead(insamp_end) - 1/12
  ) %>% 
  filter(row_number() < n()) # don't oos after opt$yearm_end


# loop
tic = Sys.time()
outlist = list()
file.remove(paste0(outroot,outfolder,'/forecast-loop.log'))
for (fi in 1:dim(sample_list)[1]){
  fcur = sample_list[fi, ]
  
  print(paste0('forecasting in ', fcur$insamp_end))  
  
  # create subsamples
  insamp = dat[yyyymm >= fcur$insamp_begin & yyyymm <= fcur$insamp_end]
  oos  = dat[yyyymm >= fcur$oos_begin & yyyymm <= fcur$oos_end]
  
  # fit and forecast
  Ebh1m = fit_and_forecast(opt$model)
  
  # save
  outlist[[fi]] = oos[,.(permno,yyyymm)][
    , Ebh1m := Ebh1m
  ]
  
  # feedback for sanity
  min_elapsed = round(as.numeric(Sys.time() - tic, units = 'mins'), 1)
  log.text <- paste0(
    Sys.time()
    , " insamp_end = ", fcur$insamp_end
    , " min elapsed = ", min_elapsed
    , " min per fit = ", round(min_elapsed / fi, 2)
    , " min remaining = ", round(min_elapsed / fi * (dim(sample_list)[1]-fi+1), 1)
  )
  write.table(log.text, paste0(outroot,outfolder, '/forecast-loop.log')
              , append = TRUE, row.names = FALSE, col.names = FALSE)
}

# bind output
forecast = rbindlist(outlist)


#==============================================================================# 
# Portfolios ----
#==============================================================================#

temp = merge(dat[, .(permno,yyyymm,bh1m)], forecast
             ,by = c('permno','yyyymm') ) 

port = temp %>% 
  group_by(yyyymm) %>% 
  mutate(port = ntile(Ebh1m,10)) %>% 
  mutate(port = sprintf('%02.0f',port)) %>% 
  group_by(port,yyyymm) %>% 
  summarize(nstock = sum(!is.na(bh1m)), bh1m = mean(bh1m, na.rm=T))

LSport = port %>% filter(port %in% c('01','10'))   %>% 
  pivot_wider(
    id_cols = 'yyyymm', names_from = 'port'
    , values_from = c('bh1m','nstock'), names_prefix = 'port'
  ) %>% 
  mutate(
    bh1m = bh1m_port10 - bh1m_port01
    , nstock = nstock_port10 + nstock_port01
  ) %>% 
  transmute(port = 'LS', yyyymm, bh1m, nstock)

port = port %>% rbind(LSport)


#==============================================================================# 
# Save to disk ----
#==============================================================================#


# save forecasts
fwrite(forecast, paste0(outroot,outfolder,'/permno-month-forecast.csv'))

# save portfolios
fwrite(port, paste0(outroot,outfolder,'/portsort.csv'))

# save settings and summary
sink(paste0(outroot,outfolder,'/summary and settings.txt'))
print('================================')
print('port sumstats')
port %>% group_by(port) %>% 
  summarize(
    rbar = mean(bh1m), vol = sd(bh1m), SRann = rbar/vol*sqrt(12)
    , nstock = mean(nstock)
    , nmonth = n()
  ) %>% 
  filter(port != 'NA')
cat('\n\n')

print('================================')
print('opt')
format(opt)
cat('\n\n')

print('================================')
print('optranger')
format(optranger)
cat('\n\n')

print('================================')
print('optkeras')
format(optkeras)
cat('\n\n')


print('================================')
print('signals_used')
as.data.frame(signals_list)
cat('\n\n')

sink()

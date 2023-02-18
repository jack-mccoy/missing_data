# 2023 02 trying ml stuff on chen-Zimmermann data

#==============================================================================#
# Setup ----
#==============================================================================#
rm(list = ls())
repulldata = T # use F if debugging and impatient

library(data.table)
library(tidyverse)
library(ggplot2)
library(zoo)
library(janitor)

# ml libraries
library(ranger)
library(tensorflow)
library(keras)
library(lightgbm)

# model independent options
opt = list(
  model = 'keras4' # 'ranger' or 'lm' or 'lightgbm' or 'keras1'-'keras4'
  , data_file = '../output/bcsignals/bcsignals_em.csv'
  , yearm_begin = '1995-06' # insamp_end begin
  , yearm_end = '2020-06'  
  , refit_period = 12 # months  
  , insamp_months_min = 12*10 
  , validate_months = 12*5  
  , window_type = 'expanding' # rolling or expanding
  , cores_frac = 0.5
  , verbose = 1
)

# output names
customname = NULL
outroot = '../output/forecast/'

# settings for specific algos
# max.depth = NULL dramatically slows things down
optranger = list( 
  num.tree = 300, max.depth = c(1,3,6), mtry = c(3,10,30,50)
)

optkeras = list(
  lambda1 = c(1e-5, 1e-3)
  , learning_rate = c(0.001, 0.01)
  , batch_size = 10000, epochs = 100, patience = 5, ensemble_num = 10
)

optlightgbm = list(
  objective = 'regression' # can't seem to get huber to work ok on sims
  , num_tree = c(500, 750, 1000) # allowing num_tree < 500 ends up performing pretty bad
  , learning_rate = c(0.01, 0.10) # default 0.1 (a.k.a. shrinkage)
  , max_depth = c(1,2) # default -1 and <= 0 implies no limit
  , force_col_wise = TRUE # TRUE suppresses warning (something about multi-threading)
)

optlm = list(
  blank = NA_real_
)

# keep this subset for testing, use NULL for keep all
signals_keep = NULL # e.g. c('size','bm','mom12m')


#==============================================================================#
# Data pull ----
#==============================================================================#

if (repulldata) {# Read in data 
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
}


#==============================================================================#
# Auto Setup ----
#==============================================================================#
# uses info from data pull, so is after data pull


if (substr(opt$model,1,5)=='keras'){
  optsuffix = 'keras'
} else {
  optsuffix = opt$model
}

# create optmodel and make list of tuning settings
evalme = paste0('optmodel = opt', optsuffix) 
eval(parse(text = evalme))
tuneset = optmodel %>% expand.grid(stringsAsFactors = F) %>% arrange(across(everything()))


# replace null options with defaults
if (is.null(optmodel$ensemble_num)){
  optmodel$ensemble_num = 1
}

# create folder
if (is.null(customname)) {
  outfolder = paste(
    opt$model
    , str_remove(basename(opt$data_file), '.csv')
    , sep = '-'
  )
} else {
  outfolder = customname
}
dir.create(paste0(outroot), showWarnings = F)
dir.create(paste0(outroot,outfolder), showWarnings = F)

# save settings (make sure you're running the right thing!)
sink(paste0(outroot,outfolder,'/settings.txt'))
print('================================')
print('opt')
format(opt)
cat('\n\n')

print('================================')
print('optmodel')
format(optmodel)
cat('\n\n')


print('================================')
print('signals_used')
as.data.frame(signals_list)
cat('\n\n')

sink()

# print settings (make sure you're running the right thing!)
print('================================')
print('opt')
format(opt)
cat('\n\n')

print('================================')
print('optmodel')
format(optmodel)
cat('\n\n')


print('================================')
print('signals_used')
print(signals_list)



#==============================================================================# 
# fit_and_forecast function ----
#==============================================================================#

# function for choosing forecast
fit_and_forecast = function(modelname, hyperpar, fit_start, fit_end, fcast_start, fcast_end, seed = NULL){
  
  # create subsamples
  insamp = dat[yyyymm >=  fit_start & yyyymm <= fit_end]
  oos  = dat[yyyymm >=  fcast_start & yyyymm <= fcast_end]  
  
  if (modelname == 'lm'){
    form = paste0('bh1m ~ ', paste(signals_list, collapse = '+'))  
    fit = lm(form, data=insamp)
    Ebh1m = predict(fit, oos)
    
  } else if (modelname == 'ranger'){
    
    fit = ranger(x = as.matrix(insamp %>% select(all_of(signals_list)))
                 , y = as.matrix(insamp %>% select(bh1m))
                 , num.trees = hyperpar$num.tree
                 , num.threads = floor(parallel::detectCores()*opt$cores_frac)
                 , max.depth = hyperpar$max.depth
                 , mtry = hyperpar$mtry
    )
    
    pred.ra = predict(fit, oos) 
    Ebh1m = pred.ra$predictions
    
  } else if (modelname == 'lightgbm'){
    
    tempparam = hyperpar[names(hyperpar) != 'ensemble_num']
    tempparam$num_threads = floor(parallel::detectCores()*opt$cores_frac)
    fit = lightgbm(
      data = as.matrix(insamp %>% select(all_of(signals_list)))
      , label = as.matrix(insamp %>% select(bh1m))
      , params = tempparam
      , verbose = opt$verbose
    )
    
    # predict
    Ebh1m = predict(fit, as.matrix(oos %>% select(all_of(signals_list))))
    
    
  } else if (substr(modelname,1,5) == 'keras') {
    # neural network has 5 different types
    
    # build model
    #   see Gu Kelly Xiu page 22 / 2244
    nlayer = substr(modelname,6,6) %>% as.numeric
    if (nlayer == 1){
      keras_model = keras_model_sequential() %>%
        layer_dense(32, activation = 'relu'
                    , kernel_regularizer = regularizer_l1(hyperpar$lambda1)) %>%
        layer_dense(1)
    } else if (nlayer == 2){
      keras_model = keras_model_sequential() %>%
        layer_dense(32, activation = 'relu'
                    , kernel_regularizer = regularizer_l1(hyperpar$lambda1)) %>%
        layer_dense(16, activation = 'relu'
                    , kernel_regularizer = regularizer_l1(hyperpar$lambda1)) 
        layer_dense(1)      
    } else if (nlayer == 3){

      keras_model = keras_model_sequential() %>%
        layer_dense(32, activation = 'relu'
                    , kernel_regularizer = regularizer_l1(hyperpar$lambda1)) %>%
        layer_dense(16, activation = 'relu'
                    , kernel_regularizer = regularizer_l1(hyperpar$lambda1)) %>%
        layer_dense(8, activation = 'relu'
                    , kernel_regularizer = regularizer_l1(hyperpar$lambda1))
        layer_dense(1)      
    } else if (nlayer == 4){
      keras_model = keras_model_sequential() %>%
        layer_dense(32, activation = 'relu'
                    , kernel_regularizer = regularizer_l1(hyperpar$lambda1)) %>%
        layer_dense(16, activation = 'relu'
                    , kernel_regularizer = regularizer_l1(hyperpar$lambda1)) %>%
        layer_dense(8, activation = 'relu'
                    , kernel_regularizer = regularizer_l1(hyperpar$lambda1)) %>%
        layer_dense(4, activation = 'relu'
                    , kernel_regularizer = regularizer_l1(hyperpar$lambda1)) %>%
        layer_dense(1)      
    } else if (nlayer == 5){
      keras_model = keras_model_sequential() %>%
        layer_dense(32, activation = 'relu'
                    , kernel_regularizer = regularizer_l1(hyperpar$lambda1)) %>%
        layer_dense(16, activation = 'relu'
                    , kernel_regularizer = regularizer_l1(hyperpar$lambda1)) %>%
        layer_dense(8, activation = 'relu'
                    , kernel_regularizer = regularizer_l1(hyperpar$lambda1)) %>%
        layer_dense(4, activation = 'relu'
                    , kernel_regularizer = regularizer_l1(hyperpar$lambda1)) %>%
        layer_dense(2, activation = 'relu'
                    , kernel_regularizer = regularizer_l1(hyperpar$lambda1)) %>%        
        layer_dense(1)      
    }
    
    keras_model %>%  compile(
      loss = 'mean_squared_error',
      optimizer = optimizer_adam(
        learning_rate = hyperpar$learning_rate
      )
    )
    
    # fit
    if (!is.null(seed)){
      tensorflow::set_random_seed(seed)
    }
    history <- keras_model %>% fit(
      as.matrix(insamp %>% select(all_of(signals_list))),
      as.matrix(insamp %>% select(bh1m)),
      verbose = opt$verbose,
      epochs = hyperpar$epochs,
      batch_size = hyperpar$batch_size,
      callbacks = callback_early_stopping(
        patience = hyperpar$patience, monitor = 'val_loss'
      ),
      validation_data = list(
        as.matrix(oos %>% select(all_of(signals_list))),
        as.matrix(oos %>% select(bh1m))
      )
    )
    
    Ebh1m = predict(
      keras_model, as.matrix(oos %>% select(all_of(signals_list)))
      , verbose = opt$verbose
    )
    
  }
  
  # output data.table
  oos_small = copy(oos[,.(permno,yyyymm,bh1m)])
  oos_small[, Ebh1m := Ebh1m]
  
  return(oos_small)
}  



#==============================================================================# 
# Loop over sample sets ----
#==============================================================================#

# create list of sample sets
sample_list = data.frame(
  insamp_end = seq(
    as.yearmon(opt$yearm_begin), as.yearmon(opt$yearm_end), opt$refit_period/12
  )
) %>% 
  mutate(
    oos_begin = insamp_end + 1/12
    , oos_end = oos_begin + opt$refit_period/12 - 1/12
    , tune_split = insamp_end - opt$validate_months/12
  ) 
  

# make expanding 
if (opt$window_type == 'rolling'){
  sample_list = sample_list %>% 
    mutate(insamp_begin = insamp_end - opt$insamp_months_min/12)
} else if (opt$window_type == 'expanding'){
  sample_list = sample_list %>% 
    mutate(insamp_begin = min(insamp_end) - opt$insamp_months_min/12)
}
sample_list = sample_list %>% 
  select(insamp_begin, tune_split, everything())


# loop
tic = Sys.time()
outlist = list()
tunelist = list()
file.remove(paste0(outroot,outfolder,'/forecast-loop.log'))
for (sampi in 1:dim(sample_list)[1]){
  sampcur = sample_list[sampi, ]
  
  print(paste0('forecasting in ', sampcur$insamp_end))  
  
  # tune hyperparbest (or not)
  if (dim(tuneset)[1] > 1){
    tunesum = tibble()
    for (seti in 1:dim(tuneset)[1]){
      
      print(paste('tuning set', seti, 'of', dim(tuneset)[1], sep = ' '))
      
      # fit before split, forecast after split
      oostemp = fit_and_forecast(opt$model, tuneset[seti,]
                                 , sampcur$insamp_begin, sampcur$tune_split
                                 , sampcur$tune_split+1/12, sampcur$insamp_end)
      
      sumtemp = tuneset[seti,] %>% cbind(
        oostemp %>% summarize(rmse = sqrt(mean((bh1m-Ebh1m)^2))) 
      )
      
      tunesum = rbind(tunesum, sumtemp)
    }
    
    tunesum = tunesum %>% 
      mutate(best = rmse == min(rmse)) %>% 
      arrange(rmse) 
      
    hyperparbest = tunesum %>% filter(best) %>% 
      arrange(across(everything())) %>% 
      filter(row_number() == 1) # break ties
    
  } else {
    # no tuning requested
    tunesum = tuneset %>% mutate(rmse = NA_real_, best = TRUE)
    hyperparbest = tuneset
    
  }
  
  # fit and forecast full in-samp with hyperparbest + ensembling
  oos_ensemb = list()
  for (ensembi in 1:optmodel$ensemble_num){
    print(paste('ensemble', ensembi, 'of', optmodel$ensemble_num, sep = ' '))
    seed = 316*ensembi
    oos_ensemb[[ensembi]] = fit_and_forecast(opt$model, hyperparbest
                              , sampcur$insamp_begin, sampcur$insamp_end
                              , sampcur$oos_begin, sampcur$oos_end
                              , seed)  
    oos_ensemb[[ensembi]]$seed = seed
  }
  oos_ensemb = rbindlist(oos_ensemb)
  ooscur = oos_ensemb[
    , .(Ebh1m = mean(Ebh1m), bh1m = mean(bh1m)), by = c('permno','yyyymm')
  ]
  
  
  # save
  outlist[[sampi]] = ooscur
  tunelist[[sampi]] = tunesum %>% mutate(
    oos_begin = sampcur$oos_begin
  )
  
  # feedback 
  min_elapsed = round(as.numeric(Sys.time() - tic, units = 'mins'), 1)
  log.text <- paste0(
    "ETA = ", round(min_elapsed / sampi * (dim(sample_list)[1]-sampi+1), 1), " min"
    , " min per fit = ", round(min_elapsed / sampi, 2)    
    , " insamp_end = ", sampcur$insamp_end
    , " min elapsed = ", min_elapsed
    , " datetime = ", Sys.time()
  )
  write.table(log.text, paste0(outroot,outfolder, '/forecast-loop.log')
              , append = TRUE, row.names = FALSE, col.names = FALSE)
  print(tunesum)

}

# bind output
forecast = rbindlist(outlist)
tuneresult = rbindlist(tunelist)


#==============================================================================# 
# Portfolios ----
#==============================================================================#

make_ports = function(sweight = 'ew'){
  
  if (sweight == 'ew'){
    temp = forecast %>% mutate(weight = 1)
  } else if (sweight == 'vw'){
    temp = forecast %>% 
      left_join(dat %>% transmute(permno,yyyymm, weight = me)
                , by = c('permno','yyyymm')) 
  }
  
  port = temp %>% 
    filter(!is.na(weight+bh1m)) %>% 
    group_by(yyyymm) %>% 
    mutate(port = ntile(Ebh1m,10)) %>% 
    mutate(port = sprintf('%02.0f',port)) %>% 
    group_by(port,yyyymm) %>% 
    summarize(nstock = sum(!is.na(bh1m))
              , bh1m = weighted.mean(bh1m,weight)
    )
  
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
  
  port = port %>% rbind(LSport) %>% mutate(sweight = sweight)
  
  port
}

port = rbind(make_ports('ew'), make_ports('vw'))
  


#==============================================================================# 
# Save to disk ----
#==============================================================================#

# save tuning results
tuneresult = tuneresult %>% arrange(-best,oos_begin)
fwrite(tuneresult, paste0(outroot,outfolder,'/tuning-results.csv'))

# save forecasts
fwrite(forecast, paste0(outroot,outfolder,'/permno-month-forecast.csv'))

# save portfolios
fwrite(port, paste0(outroot,outfolder,'/portsort.csv'))


# save settings and summary
sink(paste0(outroot,outfolder,'/summary and settings.txt'))
print('================================')
print('ew port sumstats')
port %>% 
  filter(sweight == 'ew') %>% 
  group_by(port) %>% 
  summarize(
    rbar = mean(bh1m), vol = sd(bh1m), SRann = rbar/vol*sqrt(12)
    , nstock = mean(nstock)
    , nmonth = n()
  ) %>% 
  filter(port != 'NA') 
cat('\n\n')

print('================================')
print('vw port sumstats')
port %>% 
  filter(sweight == 'vw') %>% 
  group_by(port) %>% 
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
print('optmodel')
format(optmodel)
cat('\n\n')


print('================================')
print('signals_used')
as.data.frame(signals_list)
cat('\n\n')

sink()


# one signal strategies.  Created 2022 06 Andrew

# Environment ----
rm(list = ls())
library(data.table)
library(tidyverse)
library(ggplot2)

datebegin = 1985
dateend   = 2020

# fix number of stocks in each leg to "control" for "signal loading" of mvn vs simple strats
nstock_per_leg = 500 

# read in data ----

## read ----
# round yyyymm to make sure they match (needed because the write precision differs?)

# lead returns by 1 month
ret_dat <- fread("../data/crsp_data.csv") %>% select(permno,yyyymm,ret) %>% 
  mutate(yyyymm = round(yyyymm-1/12, 3)) %>% 
  rename(bh1m = ret) 

me_dat = fread("../data/crsp_data.csv") %>% select(permno,yyyymm,me) %>% 
  mutate(yyyymm = round(yyyymm, 3)) 

# signals 
signals_mvn = fread('../output/imp_tmp.csv') %>% 
  mutate(yyyymm = round(yyyymm, 3)) 

signals_simple = fread('../output/bc_tmp.csv') %>% 
  mutate(yyyymm = round(yyyymm, 3))

## auto setup ----

# create a list of signals
signallist = signals_simple %>% select(-c(permno,yyyymm)) %>% colnames() %>% sort()
nsignal = dim(signals_simple)[2] - 2

# reg_dat has ret leaded by 1 month, everything else lagged (like a regression)
reg_dat0 = ret_dat %>% 
  filter(
    yyyymm >= datebegin, yyyymm<= dateend
  ) %>% 
  left_join(
    me_dat, by = c('permno','yyyymm')
  ) %>% 
  filter(
    !is.na(bh1m), !is.na(me)
  )

# Loop over imp_type, signals  ---- 
# takes about 2 minutes
dat_all = data.table()

for (imp_type in c('simple','mvn')){
  print(paste0('imp type is ', imp_type))
  
  # merge on signals of the desired type
  if (imp_type == 'simple'){
    reg_dat = reg_dat0  %>% 
      left_join(
        signals_simple, by = c('permno','yyyymm')
      ) 
  } else {
    reg_dat = reg_dat0  %>% 
      left_join(
        signals_mvn, by = c('permno','yyyymm')
      )     
  } # if imp_type
  
  ## loop over signals
  for (signali in 1:nsignal){
    print(paste0('forming ports for signali = ', signali))
    
    # read in data and assign ports
    dat_one = reg_dat %>% 
      select(permno, yyyymm, bh1m, me, signallist[signali])  %>% 
      rename(signalcurr = signallist[signali]) %>% 
      filter(!is.na(signalcurr)) %>% 
      group_by(yyyymm) %>% 
      arrange(yyyymm,signalcurr) %>% 
      mutate(
        worst_rank = rank(signalcurr)
        , best_rank = rank(-signalcurr)
        , port = case_when(
          worst_rank <= 500 ~ 'short'
          , best_rank <= 500 ~ 'long'
          , TRUE ~ NA_character_
        )
      )
    

    # find port stats
    dat_port = dat_one %>% 
      group_by(
        yyyymm, port
      ) %>% 
      summarize(
        bh1m_ew = mean(bh1m)
        , bh1m_vw = weighted.mean(bh1m, me)
        , nstock = sum(!is.na(signalcurr))
        , .groups = 'drop'
      ) %>% 
      mutate(
        imp_type = imp_type
        , signalname = signallist[signali]
      )
    
    # compile
    dat_all = dat_all %>% 
      rbind(
        dat_port
      )
    
  } # for signali
  
} # for imp_type



# Summarize by imp_type, signal----

ret_var = 'bh1m_ew'

# find bad signal-dates (not enough stocks in the simple imp_type)
bad_signal_dates = dat_all %>% 
  filter(imp_type == 'simple') %>% 
  pivot_wider(
    id_cols = c(yyyymm,signalname), names_from = port, values_from = nstock
  ) %>% 
  mutate(
    bad = if_else(long == 500 & short == 500, 0 , 1)
  ) %>% 
  select(yyyymm,signalname,bad)



# reorganize data
# and create long short strategies, flagging signal-dates if not enough observations 
dat_ls = dat_all %>% 
  pivot_longer(
    cols = starts_with('bh1m'), names_prefix = 'bh1m_'
    , names_to = 'sweight', values_to = 'bh1m'
  )%>% 
  pivot_wider(
    names_from = port, values_from = c('bh1m','nstock')
  ) %>% 
  mutate(
    bh1m_ls = bh1m_long - bh1m_short
  ) %>% 
  filter(
    nstock_long == nstock_per_leg, nstock_short == nstock_per_leg
  ) %>% 
  left_join(
    bad_signal_dates, by = c('yyyymm','signalname')
  ) 

# check how many are bad
mean(dat_ls$bad)


# summarize by sweight, imp_type, signalname
signal_sum = dat_ls %>% 
  filter(bad == 0) %>% 
  group_by(sweight, imp_type, signalname) %>% 
  summarize(
    rbar = mean(bh1m_ls)*12
    , vol = sd(bh1m_ls)*sqrt(12)    
    , ndate = n()/12
    , tstat = rbar / vol * sqrt(ndate)
    , sr = rbar / vol 
  ) 


# turn to wide format for differences
signal_sum_wide = signal_sum %>% 
  select(-c(vol,ndate,tstat)) %>% 
  pivot_wider(
    names_from = imp_type, values_from = c('rbar','sr')
  ) %>% 
  mutate(
    rbar_diff = rbar_mvn - rbar_simple
    , sr_diff = sr_mvn - sr_simple
  )


## tables for output ====

# summary stats
qlist = c(0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95)
imp_sum = signal_sum_wide %>% 
  group_by(sweight) %>% 
  summarize(
    percentile = qlist*100
    , blank = NA*qlist
    , rbar_mvn = quantile(rbar_mvn, qlist, na.rm = TRUE)    
    , rbar_simple = quantile(rbar_simple, qlist, na.rm = TRUE)
    , rbar_mvn_less_simple = quantile(rbar_diff, qlist, na.rm = TRUE)        
    , blank2 = NA*qlist
    , sr_mvn = quantile(sr_mvn, qlist, na.rm = TRUE)    
    , sr_simple = quantile(sr_simple, qlist, na.rm = TRUE)
    , sr_mvn_less_simple = quantile(sr_diff, qlist, na.rm = TRUE)
    , blank3 = NA*qlist    
  ) 

# reorganize
imp_sum2 = imp_sum %>% 
  filter(sweight == 'ew') %>% 
  t() %>% 
  rbind(
    imp_sum %>% 
      filter(sweight == 'vw') %>% 
      t() 
  )

imp_sum2

write.csv(imp_sum2, file = '../output/plots/one_signal.csv' )



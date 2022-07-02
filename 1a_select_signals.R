# created from 1_select_signals.R to make summary stats for missing values
# this analysis uses all continuous predictors (does not limit to 
# the most observed 75 or whatever)
# Andrew 2022 06

#==============================================================================#
# Packages ----
#==============================================================================#

rm(list = ls())

library(data.table)
library(zoo)
library(tidyverse)
library(xtable) # for long table

source("functions.R")


#==============================================================================#
# Read in data ----
#==============================================================================#

# Read in the downloadable signals 
signals <- fread("../data/signed_predictors_dl_wide.csv") 
setnames(signals, colnames(signals), tolower(colnames(signals)))
signals[, yyyymm := as.yearmon(as.Date(as.character(yyyymm*10 + 1), "%Y%m%d"))]

# merge on crsp predictors
crsp_data <- fread("../data/crsp_data.csv")[, .SD, .SDcols = !c("ret", "me")]
crsp_data[, yyyymm := as.yearmon(yyyymm)]
signals <- merge(signals, crsp_data, by = c("permno", "yyyymm"))

# read signal doc
signaldoc = fread('../data/SignalDoc.csv')  %>% 
  mutate(
    signalname = Acronym
  ) %>% 
  select(
    signalname, Authors, Year, starts_with('Cat.')
  ) %>% 
  mutate(signalname = tolower(signalname)) %>% 
  filter(
    Cat.Signal == 'Predictor'
  )

# make indicator matrix of observations (this is all we need here)
obs = signals %>% mutate_at(vars(-c('yyyymm','permno')), function (x) as.numeric(!is.na(x)) )

gc()

#==============================================================================#
# Drop bad signals ----
#==============================================================================#
# Drop discrete and signals w/ gaps
# and here drop anything else we want to in the future

# user
min_obs_ok = 100
date_start = 1985
date_end = 2020

# find continuous
list_cts = signaldoc[Cat.Form == 'continuous'] %>% pull(signalname)

# find totally-missing-months
obs_by_date = obs[
  , by = yyyymm
  , lapply(.SD, sum)
  , .SDcols = !c('permno','yyyymm')
] %>% 
  pivot_longer(cols = -c('yyyymm'), names_to = 'signalname', values_to = 'obs') %>% 
  group_by(signalname)

# find start / end dates by signal (not used right now)
date_start_end = left_join(
  obs_by_date %>% 
    filter(obs > min_obs_ok) %>% 
    arrange(signalname,yyyymm) %>% 
    filter(row_number() == 1) %>% 
    transmute(signalname, date_first_obs = yyyymm)
  , obs_by_date %>% 
    filter(obs > min_obs_ok) %>% 
    arrange(signalname,-yyyymm) %>% 
    filter(row_number() == 1) %>% 
    transmute(signalname, date_last_obs = yyyymm)
)


# find gaps
gaps_all = obs_by_date %>% 
  left_join(
    date_start_end, by = 'signalname'
  )  %>% 
  filter(
    yyyymm >= date_start, yyyymm <= date_end, obs < 2
  ) %>% 
  arrange(signalname, yyyymm) %>% 
  print(n=100)

list_gap = gaps_all %>% distinct(signalname, .keep_all = T) %>% pull(signalname)

# out to console
list_cts
intersect(list_cts,list_gap)

# remove bad guys
obs2 = obs %>% select(permno,yyyymm,all_of(list_cts)) %>% select(-intersect(list_cts,list_gap))

length(names(obs2)) - 2


#==============================================================================#
# Find Stats ----
#==============================================================================#

## Make small dataset ----

# keep only representative years
datelist = as.yearmon(c(
  'Jun 1970'
  , 'Jun 1980', 'Jun 1985'
  , 'Jun 1990', 'Jun 1995'
  , 'Jun 2000', 'Jun 2005'
  , 'Jun 2010', 'Jun 2015'
))

# make small dataset for sum stats
small = obs2 %>% filter(yyyymm %in% datelist) %>% 
  select(-permno) %>% filter(streversal == 1)


## One signal summaries ----

# pctiles
ptile = c(5, 10, 25, 50, 75, 90, 95)


# one-signal sum
nobs_by_signal = small %>% 
  group_by(yyyymm) %>% 
  summarize_all(sum) %>% 
  pivot_longer(
    cols = -c(yyyymm,streversal), names_to = 'signalname', values_to = 'nobs'
  ) %>% 
  mutate(
    pct_obs = nobs/streversal*100
  )

# table for output
pct_obs_one_signal_ptile = nobs_by_signal %>% 
  group_by(yyyymm) %>% 
  summarize(
    ptile = ptile
    , pct_obs = quantile(pct_obs, ptile/100)
  ) %>% 
  pivot_wider(
    names_from = ptile, values_from = pct_obs
  ) %>% t()

pct_obs_one_signal_ptile 

write.csv(pct_obs_one_signal_ptile, '../output/pct_obs_one_signal_ptile.csv')




## Many signal summaries ----
datelist2 = as.yearmon(c('Jun 1990', 'Jun 2000', 'Jun 2010'))


# focus on sets of predictors determined by obs in 1985
nobs_by_signal_overall = small %>% 
  filter(yyyymm == 'Jun 1985') %>% 
  select(-yyyymm) %>% 
  summarize_all(sum) %>% 
  pivot_longer(
    cols = -c(streversal), names_to = 'signalname', values_to = 'nobs'
  ) %>% 
  mutate(
    pct_obs = nobs/streversal*100
  ) %>% 
  arrange(
    -pct_obs
  ) %>% 
  mutate(
    rank = row_number()
  )

# multiple signal counting
dat = tibble()

for (datecurr in datelist){
  
  obs_curr = small %>% 
    filter(yyyymm == datecurr) %>% 
    select(-c(yyyymm))
  
  for (J in c(5,25,50,75,100,125, 150)){
    
    # J is desired number of signals
    # We keep the J most observed signals 
    signallist = nobs_by_signal_overall %>% filter(rank <= J) %>% pull(signalname) %>% 
      unique()
    
    # then count the multiple signal "situation" in the given date
    obs_select = obs_curr %>% select(signallist) 
    
    pct_cc = sum(apply(obs_select, 1, prod)) / dim(obs_select)[1] * 100
    pct_obs = sum(obs_select) / prod(dim(obs_select)) * 100
    
    # save
    dat_mini = tibble(datecurr, J, pct_cc, pct_obs)
    dat = rbind(dat, dat_mini)
    
  } # for J
  
} # for datecurr


# reshape for output
dat_long = dat %>% 
  pivot_longer(
    cols = starts_with('pct'), names_to = 'type', values_to = 'pct', names_prefix = 'pct_'
  )


pct_obs_by_nsignal = dat_long %>% 
  mutate(
    datecurr = paste0('y', floor(datecurr))
  ) %>% 
  pivot_wider(
    names_from = datecurr, values_from = pct
  ) %>% 
  arrange(J,desc(type)) %>% 
  pivot_wider(
    names_from = type, values_from = starts_with('y')
  )


write.csv(pct_obs_by_nsignal, '../output/pct_obs_by_nsignal.csv')

pct_obs_by_nsignal

#==============================================================================#
# Output final lists of signals ----
#==============================================================================#

## find lists ----
list_1985 = nobs_by_signal %>% 
  filter(year(yyyymm) == 1985) %>% 
  group_by(signalname) %>% 
  summarize(
    mean_pct_obs = mean(pct_obs)
  ) %>% 
  arrange(-mean_pct_obs) %>% 
  mutate(
    rank = row_number()
    , best100 = if_else(row_number() <= 100, 1, 0)
  ) %>% 
  select(signalname, mean_pct_obs, rank, best100) 

list_full = nobs_by_signal %>% 
  filter(yyyymm >= date_start, yyyymm <= date_end) %>% 
  group_by(signalname) %>% 
  summarize(
    mean_pct_obs = mean(pct_obs)
  ) %>% 
  arrange(-mean_pct_obs) %>% 
  mutate(
    rank = row_number()
    , best100 = if_else(row_number() <= 100, 1, 0)
  ) %>% 
  select(signalname, mean_pct_obs, rank, best100) 


# merge to documentation for checking and tables (add more info later)
signaldoc2 = signaldoc %>% 
  left_join(list_1985 %>% transmute(signalname, rank1985 = rank) ) %>% 
  left_join(list_full %>% transmute(signalname, rankfull = rank) )
  

## output ----
writeLines(list_1985$signalname, '../data/signals_best100_1985.txt')

writeLines(list_full$signalname, '../data/signals_best100_full.txt')

write.csv(signaldoc2, '../data/signals_best_doc.csv', row.names = F)




#==============================================================================#
# Check (to console, for now) ----
#==============================================================================#

doc3 = fread('../data/signals_best_doc.csv')

signaldoc %>% filter(Cat.Data == 'Analyst') %>% 
  left_join(
    list_1985
  ) %>% 
  select(signalname, mean_pct_obs, rank, best100)



datasum = doc3 %>% filter(best100_1985 == 1) %>%group_by(Cat.Data) %>% summarize(nsignal_1985 = n()) %>% 
  left_join(
    doc3 %>% filter(best100_full == 1) %>%group_by(Cat.Data) %>% summarize(nsignal_full = n()) 
  )   %>% 
  left_join(
    doc3  %>% group_by(Cat.Data) %>% summarize(nsignal_all = n()) 
  ) %>% 
  print(n=Inf)

econsum = doc3 %>% filter(best100_1985 == 1) %>%group_by(Cat.Economic) %>% summarize(nsignal_1985 = n()) %>% 
  left_join(
    doc3 %>% filter(best100_full == 1) %>%group_by(Cat.Economic) %>% summarize(nsignal_full = n()) 
  )   %>% 
  left_join(
    doc3  %>% group_by(Cat.Economic) %>% summarize(nsignal_all = n()) 
  ) %>% 
  print(n=100)


keysignals = c(
  'accruals','bmdec','mom12m','assetgrowth'
  ,'earningssurprise', 'analystrevision', 'feps'
  , 'gp','roaq'
  , 'smileslope', 'skew1'
)

keysum = doc3 %>% filter(signalname %in% keysignals) %>% 
  select(signalname, starts_with('best')) %>% 
  print(n=100)


#==============================================================================#
# Tables for paper ----
#==============================================================================#


# Create Latex output table 1: Clear Predictors
outputtable1 <- xtable(final_list %>%
                         filter(rank <= 75) )


print(outputtable1,
      include.rownames = FALSE,
      include.colnames = FALSE,
      hline.after = NULL,
      only.contents = TRUE,
      file = '../output/final_list.tex')

# Check on stuff manually ----

## list signal names for examples  ----
nobs_by_signal %>% 
  filter(yyyymm == 'Jun 2000') %>% 
  arrange(pct_obs)


nobs_by_signal %>% 
  filter(yyyymm == 'Jun 2000') %>% 
  arrange(-pct_obs)

nobs_by_signal %>% 
  filter(yyyymm == 'Jun 2000') %>% 
  filter(pct_obs <= 65, pct_obs >= 60)

## list all in the June 1980 reference month ----
final_list %>% as.data.frame()

## find when particular predictors appear in the data ----

signals %>% 
  select(permno,yyyymm,earningsstreak,sfe,forecastdispersion,feps) %>% 
  mutate(
    yyyymm = as.numeric(yyyymm)
  ) %>% 
  filter(
    yyyymm >= 1970, yyyymm%%1 == 0
  ) %>% 
  group_by(yyyymm) %>% 
  summarize_all(
    function (x) sum(!is.na(x))
  ) %>% 
  mutate_at(
    vars(-c('yyyymm','permno')), funs(./permno*100)
  ) %>% 
  print(n=50)

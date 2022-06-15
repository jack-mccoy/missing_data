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

# keep only representative years
datelist = as.yearmon(c('Jun 1970', 'Jun 1980', 'Jun 1990', 'Jun 2000', 'Jun 2010'))


#==============================================================================#
# Read in data ----
#==============================================================================#

# Read in the signals we want 
signals <- fread("../data/signed_predictors_dl_wide.csv") 
setnames(signals, colnames(signals), tolower(colnames(signals)))

# Ease of use
signals[, yyyymm := as.yearmon(as.Date(as.character(yyyymm*10 + 1), "%Y%m%d"))]

crsp_data <- fread("../data/crsp_data.csv")[, .SD, .SDcols = !c("ret", "me")]
crsp_data[, yyyymm := as.yearmon(yyyymm)]

# merge crsp signals onto downloadable signals
signals <- merge(signals, crsp_data, by = c("permno", "yyyymm"))

# read signal doc
signaldoc = fread('../data/SignalDoc.csv')  %>% 
  transmute(
    signalname = Acronym, Cat.Form, Authors, Year, Cat.Signal
  ) %>% 
  mutate(signalname = tolower(signalname))
  

#==============================================================================#
# Process ----
#==============================================================================#

# keep only continuous predictor signals
signallist = signaldoc %>% filter(Cat.Signal == 'Predictor', Cat.Form =='continuous') %>% 
  pull(signalname)

small = signals %>% 
  filter(yyyymm %in% datelist) %>% 
  select(c('permno', 'yyyymm', all_of(signallist)))

# make indicator matrix of observations
obs = small %>% 
  mutate_at(vars(-c('yyyymm')), function (x) as.numeric(!is.na(x)) )

#==============================================================================#
# Find Stats and Output ----
#==============================================================================#

## One signal summaries ----

# pctiles
ptile = c(5, 10, 25, 50, 75, 90, 95)


# one-signal sum
nobs_by_signal = obs %>% 
  group_by(yyyymm) %>% 
  summarize_all(sum) %>% 
  pivot_longer(
    cols = -c(yyyymm,permno), names_to = 'signalname', values_to = 'nobs'
  ) %>% 
  mutate(
    pct_obs = nobs/permno*100
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


# focus on sets of predictors determined by obs in 1980
nobs_by_signal_overall = obs %>% 
  filter(yyyymm == 'Jun 1980') %>% 
  select(-yyyymm) %>% 
  summarize_all(sum) %>% 
  pivot_longer(
    cols = -c(permno), names_to = 'signalname', values_to = 'nobs'
  ) %>% 
  mutate(
    pct_obs = nobs/permno*100
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

  obs_curr = obs %>% 
    filter(yyyymm == datecurr) %>% 
    select(-c(permno,yyyymm)) %>% 
    filter(streversal==1)
  
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

## Output final list of signals ----


final_list = nobs_by_signal_overall %>% 
  left_join(
    signaldoc, by = 'signalname'
  ) %>% 
  mutate(
    pct_obs_1980 = pct_obs
    , AuthorYear = paste0(Authors, ' ', Year)
  ) %>% 
  select(
    signalname, AuthorYear, pct_obs_1980, rank
  )
  

final_list

write.csv(final_list, '../output/final_list.csv', row.names = F)

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

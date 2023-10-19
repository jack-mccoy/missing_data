# created from 1_select_signals.R to make summary stats for missing values
# this analysis uses all continuous predictors (does not limit to 
# the most observed 75 or whatever)
# Andrew 2022 06

# Update 2023 04: now we use the same 125 signals throughout for ease of reading.

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
# Hardcodes ----
#==============================================================================#

# user
date_start <- 1985
date_end <- 2022

#==============================================================================#
# Read in, create obs table, and select 125 signals  ----
#==============================================================================#

# Get the file paths
getFilePaths(signal_list_exists = FALSE) # signal list not yet made

## read ----

# Read in the downloadable signals 
signals <- fread(paste0(FILEPATHS$data_path, 
    "raw/signed_predictors_dl_wide_filtered.csv")) 
setnames(signals, colnames(signals), tolower(colnames(signals)))
signals[, yyyymm := as.yearmon(yyyymm)] # NOTE: data is pre-formatted now in 0a_download_data.R

# merge on crsp predictors
crsp_data <- fread(paste0(FILEPATHS$data_path, "raw/crsp_data.csv"))[, .SD, .SDcols = !c("ret", "me")]
crsp_data[, yyyymm := as.yearmon(yyyymm)]
signals <- merge(signals, crsp_data, by = c("permno", "yyyymm"))

# read signal doc
signaldoc = fread(paste0(FILEPATHS$data_path, 'raw/SignalDoc.csv'))  %>% 
  mutate(
    signalname = Acronym
  ) %>% 
  select(
    signalname, Authors, Year, starts_with('Cat.'), starts_with('Sample')
  ) %>% 
  mutate(signalname = tolower(signalname)) %>% 
  filter(
    Cat.Signal == 'Predictor'
  ) %>% 
  # fix typo for earningssurprise
  mutate(
    Cat.Data = if_else(signalname == 'earningssurprise', 'Accounting', Cat.Data)
  )

# make indicator matrix of observations (this is all we need here)
obs = signals %>% 
  mutate_at(vars(-c('yyyymm')), function (x) as.numeric(!is.na(x)) ) %>% 
  # keep only if streversal is observed
  filter(streversal == 1)

## drop discrete ----
list_cts = signaldoc[Cat.Form == 'continuous'] %>% pull(signalname)
obs = obs %>% select(yyyymm,permno,all_of(list_cts))

## drop signals with gaps ----

# find totally-missing-months
obs_by_date = obs[
  , by = yyyymm
  , lapply(.SD, sum)
  , .SDcols = !c('permno','yyyymm')
] %>% 
  pivot_longer(cols = -c('yyyymm'), names_to = 'signalname', values_to = 'obs') %>% 
  group_by(signalname)

# find gaps
gaps_all = obs_by_date %>% 
  filter(
    yyyymm >= date_start, yyyymm <= date_end, obs < 2
  ) %>% 
  arrange(signalname, yyyymm) 

list_gap = gaps_all %>% distinct(signalname, .keep_all = T) %>% pull(signalname)

obs = obs %>% select(-c(list_gap)) 

gc()

## keep only most obs in 1985 ----

# we do this double averaging for historical reasons
# even though a single average is probably just fine.

 
mostobs1985 = obs[floor(yyyymm) == 1985]  %>% 
  group_by(yyyymm) %>% 
  summarize_all(sum) %>% 
  pivot_longer(
    cols = -c(yyyymm,permno), names_to = 'signalname', values_to = 'nobs'
  ) %>% 
  mutate(
    pct_obs = nobs/permno*100
  )   %>% 
  group_by(signalname) %>% 
  summarize(
    mean_pct_obs = mean(pct_obs)
  ) %>% 
  arrange(-mean_pct_obs) %>% 
  mutate(rank = row_number())  %>% 
  #filter(rank <= 125) %>% 
  setDT()
  
obs = obs %>% select(permno,yyyymm,all_of(mostobs1985$signalname))


#==============================================================================#
# Final lists of signals ----
#==============================================================================#

# merge to documentation for checking and tables (add more info later)
signaldoc2 = signaldoc %>% 
  inner_join(mostobs1985 %>% transmute(signalname, pct_obs_1985 = mean_pct_obs, rank1985 = rank) ) %>% 
  arrange(rank1985)
  
writeLines(
  signaldoc2 %>% filter(rank1985 <= 125) %>% pull(signalname),
  paste0(FILEPATHS$data_path, 'signals_best125_1985.txt')
)
  
writeLines(
  signaldoc2 %>% filter(rank1985 <= 150) %>% pull(signalname), 
  paste0(FILEPATHS$data_path, 'signals_best150_1985.txt')
)

writeLines(
  signaldoc2 %>% pull(signalname),
  paste0(FILEPATHS$data_path, 'signals_all_1985.txt')
)

#==============================================================================#
# Output lists of predictors ----
#==============================================================================#

dir.create(paste0(FILEPATHS$out_path, "tables/"))

final_list = signaldoc2 %>% 
  mutate(
    AuthorYear = paste0(Authors, ' ', Year)
  ) %>% 
  select(
    signalname, AuthorYear, pct_obs_1985, rank1985
  ) %>% 
  arrange(rank1985) 

# Create Latex output table 1: Clear Predictors
outputtable1 <- xtable(
  final_list 
  , digits = c(0,0,0,1,0)
)

print(outputtable1,
      include.rownames = FALSE,
      include.colnames = FALSE,
      hline.after = NULL,
      only.contents = TRUE,
      file = paste0(FILEPATHS$out_path, 'tables/long_list.tex'))

#==============================================================================#
# Check (to console, for now) ----
#==============================================================================#

## data coverage ----
signaldoc2  %>% group_by(Cat.Data) %>% summarize(nsignal_1985 = n()) 

## economic coverage ----
signaldoc2  %>% group_by(Cat.Economic) %>% summarize(nsignal_1985 = n()) %>% 
  print(n=50)

## select coverage ----

selectsignals = c(
  'accruals','bmdec','mom12m','assetgrowth'
  ,'earningssurprise', 'analystrevision', 'feps'
  , 'gp','roaq'
  , 'smileslope', 'skew1'
)

keysum = signaldoc2 %>% filter(signalname %in% selectsignals) %>% 
  select(signalname, starts_with('pct'), starts_with('rank')) %>% 
  arrange(rank1985) %>% 
  print(n=50) 



## analyst signals ----

signaldoc2 %>% filter(Cat.Data == 'Analyst') %>% 
  select(signalname, pct_obs_1985, rank1985) %>% 
  arrange(rank1985)


#==============================================================================#
# Tables: sum stats ----
#==============================================================================#

## Make small dataset ----

# keep only representative years
datelist = as.yearmon(c(
  'Jun 1970', 'Jun 1975'
  , 'Jun 1980', 'Jun 1985'
  , 'Jun 1990', 'Jun 1995'
  , 'Jun 2000', 'Jun 2005'
  , 'Jun 2010', 'Jun 2015'
))

# make small dataset for sum stats
small = obs %>% filter(yyyymm %in% datelist) 


## One signal summaries ----
dateselect = as.yearmon(c('Jun 1970','Jun 1975','Jun 1980','Jun 1985','Jun 1990', 'Jun 1995', 'Jun 2000'))
ptile = c(5, 25, 50, 75, 95)

# one-signal sum
nobs_by_signal = small %>% 
  group_by(yyyymm) %>% 
  summarize_all(sum) %>% 
  pivot_longer(
    cols = -c(yyyymm,permno), names_to = 'signalname', values_to = 'nobs'
  ) %>% 
  mutate(
    pct_obs = nobs/permno*100
  )

pct_obs_one_signal_ptile = nobs_by_signal %>% 
  filter(yyyymm %in% dateselect) %>% 
  group_by(yyyymm) %>% 
  summarize(
    ptile = ptile
    , pct_obs = quantile(pct_obs, ptile/100)
  ) %>%
  pivot_wider(
    names_from = ptile, values_from = pct_obs
  ) %>% t()

pct_obs_one_signal_ptile 

write.csv(pct_obs_one_signal_ptile, 
    paste0(FILEPATHS$out_path, 'tables/pct_obs_one_signal_ptile.csv'))

## Many signal summaries ----
dateselect = c(1975, 1985, 1995)
Jselect = seq(25,150,25)

# focus on sets of predictors determined by obs in 1985
nobs_by_signal_overall = small %>% 
  filter(yyyymm == 'Jun 1985') %>% 
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
nobs_by_J_signals = dat %>% 
  pivot_longer(
    cols = starts_with('pct'), names_to = 'type', values_to = 'pct', names_prefix = 'pct_'
  )


pct_obs_by_nsignal = nobs_by_J_signals %>% 
  filter(floor(datecurr) %in% dateselect, J %in% Jselect) %>% 
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


write.csv(pct_obs_by_nsignal, 
    paste0(FILEPATHS$out_path, 'tables/pct_obs_by_nsignal.csv'))

pct_obs_by_nsignal

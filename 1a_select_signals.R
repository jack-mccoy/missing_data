# outputs signals_best100_1985.txt and signals_best100_full.txt

#==============================================================================#
# Packages ----
#==============================================================================#

library(data.table)
library(optparse)
library(zoo)

#==============================================================================#
# Functions ----
#==============================================================================#

source("functions.R")

#==============================================================================#
# Count stocks ----
#==============================================================================#

## read in data ====
# doc
doc = fread('../data/SignalDoc.csv') %>% 
  rename(signalname = Acronym) %>% 
  filter(
    Cat.Signal == 'Predictor'
  ) %>% 
  select(signalname, starts_with('Cat'), -Cat.Signal) %>% 
  mutate(
    signalname = tolower(signalname)
  ) 

# portfolio dat (long only)
portret = fread('../data/PredictorPortsFull.csv') %>% 
  filter(port != 'LS') %>% 
  mutate(signalname = tolower(signalname))

# count stocks
nstock = portret %>% 
  group_by(signalname, date) %>% 
  summarize(nstock = sum(Nlong)) %>% 
  mutate(year = year(date))

nstocksum = nstock %>% 
  group_by(signalname,year) %>% 
  summarize(nstock = mean(nstock)) 


nstockselect = nstocksum %>% 
  filter(
    year %in% c(1980, 1985, 1990)
  ) %>% 
  pivot_wider(
    names_from = year, names_prefix = 'nstock_', values_from = nstock
  ) %>% 
  left_join(
    nstocksum %>% filter(year >= 1985, year <= 2020) %>% summarize(nstock_1985_2020 = mean(nstock))
  )


yearok = nstocksum %>% 
  arrange(signalname, year) %>% 
  filter(nstock >= 1000) %>% 
  group_by(signalname) %>%   
  filter(row_number() == 1) %>% 
  transmute(signalname, year_n_gt_1000 = year)

# merge
doc2 = doc %>% 
  left_join(yearok) %>% 
  left_join(nstockselect) %>% 
  mutate_at(vars(starts_with('nstock')), ~replace_na(.,0))

#==============================================================================#
# Select signals ----
#==============================================================================#


# first remove discrete and coskewacx 
#  coskewacx has some gaps due to some mysterious bug
signallist0 = doc2 %>% 
  filter(Cat.Form == 'continuous') %>% select(-Cat.Form) %>% 
  filter(signalname != 'coskewacx') 

# screen for most valid signals
list_1985 = signallist0 %>% arrange(-nstock_1985) %>% filter(row_number() <= 100) %>% 
  select(signalname) %>% mutate(best100_1985 = 1) 

list_full = signallist0 %>% arrange(-nstock_1985_2020) %>% filter(row_number() <= 100) %>% 
  select(signalname) %>% mutate(best100_full = 1) 

# merge to documentation for checking and tables
doc3 = doc2 %>% 
  left_join(list_1985) %>% 
  left_join(list_full) %>% 
  mutate_at(vars(starts_with('best')), ~replace_na(.,0))





#==============================================================================#
# Output ----
#==============================================================================#


writeLines(list_1985$signalname, '../data/signals_best100_1985.txt')

writeLines(list_1985$signalname, '../data/signals_best100_full.txt')

write.csv(doc3, '../data/signals_best_doc.csv', row.names = F)


#==============================================================================#
# Check (to console, for now) ----
#==============================================================================#

doc3 = fread('../data/signals_best_doc.csv')


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


keysignals = c('accruals','bmdec','mom12m','assetgrowth','earningssurprise', 'analystrevision', 'feps', 'gp','roaq')

keysum = doc3 %>% filter(signalname %in% keysignals) %>% 
  select(signalname, starts_with('best')) %>% 
  print(n=100)
# loop over imputation-model combinations
# check ../output/forecast/ to track progress


# Setup -------------------------------------------------------------------

library(dplyr)

# models to run
#   ranger is super slow
modlist = data.frame(
  model = c('lm','pcr','spcr','keras1','keras2','keras3','keras4','keras5','lightgbm','ranger')
) %>% 
  mutate(modid = row_number())

# imputations to use
bcroot = '../output/bcsignals/'
signal_file_list = c(
  'bcsignals_none.csv'
  ,'bcsignals_em.csv'
  , 'bcsignals_availcase.csv'
)

# first and last in-sample ends
yearm_begin = '1995-06'
yearm_end = '2020-06'


setlist = expand.grid(
  model = modlist$model
  , signal_file = paste0(bcroot, signal_file_list)
  , yearm_begin = yearm_begin
  , yearm_end = yearm_end
  , stringsAsFactors = F
) %>% 
  mutate(
    output_folder = 'auto'
  ) %>% 
  # this stuff is just to put things in order
  left_join(modlist, by = 'model') %>% 
  arrange(modid) %>% 
  select(modid, everything())

# output for sanity
print(setlist, right = F)
sink('../output/forecast/setlist.txt')
print(setlist, right = F)
sink()


# Loop over settings ------------------------------------------------------

for (seti in 1:dim(setlist)[1]){
  
  # run one setting
  # outputs to ../output/forecast/[model-bcname]/
  runme = paste0(' Rscript 5a_one_forecast.R'
                 , ' --model=', setlist[seti, 'model']
                 , ' --signal_file=', setlist[seti, 'signal_file']
                 , ' --output_folder=', setlist[seti, 'output_folder']
                 , ' --yearm_begin=', setlist[seti, 'yearm_begin']
                 , ' --yearm_end=', setlist[seti, 'yearm_end']
  )   
  print(runme)
  
  shell(runme)
  
}

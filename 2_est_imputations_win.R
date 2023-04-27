# Run imputations on windows
#   I guess a bat or powershell script would be a bit nicer but I know how to do this - ac
# not super happy with the unlink (force delete temporary files permanently) but oh well


# tbc: ppca 


# run bllp  ---------------------------------------------------------------

# this should only take like one hour

runme = paste0(
  'Rscript 2d_bllp_imp.R'
  , ' --num_PCs=', 6
  , ' --out_path=', '../output/'
  , ' --out_name=', 'bcsignals_bllp.csv'
  , ' --impute_vec=', '../output/signals.txt'
)

shell(runme)

# run ar1 em --------------------------------------------------------------


unlink("../output/em_intermediate/", recursive = T, force = T) # clean up again
for (year_cur in 2014:2019){
  
  runme = paste0(
    'Rscript 2a_ar1_em_est.R'
    , ' --em_type=', 'ar1'
    , ' --output_ts_pred=', TRUE # output the ts prediction for analysis
    , ' --impute_yr=', year_cur
    , ' --impute_vec=', '../output/signals.txt'
    , ' --ar1_sample_length=5'
    , ' --maxiter=10000'
    , ' --tol=1e-4'
    , ' --cores_frac=0.5'
  )
  
  shell(runme)
  
}

# and bind
runme = paste0('Rscript 2b_bind_em_years.R '
    , '--big_file_name=bcsignals_emar1.csv '
    , '--ts_file_name=ts_prediction_ar1.csv')
shell(runme)



# run regular em ----------------------------------------------------------


unlink("../output/em_intermediate/", recursive = T, force = T) # first clean up
for (year_cur in 1985:2020){
  
  runme = paste0(
    'Rscript 2a_ar1_em_est.R'
    , ' --em_type=', 'regular'
    , ' --output_ts_pred=', FALSE    
    , ' --impute_yr=', year_cur
    , ' --impute_vec=', '../output/signals.txt'
    , ' --ar1_sample_length=1'
    , ' --maxiter=10000'
    , ' --tol=1e-4'
    , ' --cores_frac=0.5'
  )
  
  shell(runme)
  
}

# and bind
shell(
  'Rscript 2b_bind_em_years.R --big_file_name=bcsignals_em.csv'
)




# run select to clean -----------------------------------------------------
library(tidyverse)

# •	Jan 2014
# •	March 2015, May 2015
# •	October 2016
# •	May 2017
# •	July 2018
# •	May 2019


yearm_list = list(
  2014, 1
  , 2015, '3,5'
  , 2016, 10
  , 2017, 05
  , 2018, 07
  , 2019, 05
) %>% 
  matrix(ncol = 2, byrow = T)


unlink("../output/em_intermediate/", recursive = T, force = T) # first clean up

for (yearid in 1:nrow(yearm_list)){
  
  year_cur = yearm_list[yearid,1]
  months_cur = yearm_list[yearid,2]
  
  runme = paste0(
    'Rscript 2a_ar1_em_est.R'
    , ' --em_type=', 'regular'
    , ' --output_ts_pred=', FALSE    
    , ' --impute_yr=', year_cur
    , ' --impute_months=', months_cur
    , ' --impute_vec=', '../output/signals.txt'
    , ' --ar1_sample_length=1'
    , ' --maxiter=10000'
    , ' --tol=1e-4'
    , ' --cores_frac=0.5'
  )
  
  shell(runme)
  
}
# x --------------------------------------------------
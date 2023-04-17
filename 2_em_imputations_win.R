# Run imputations on windows
#   I guess a bat or powershell script would be a bit nicer but I know how to do this - ac
# not super happy with the unlink (force delete temporary files permanently) but oh well



# run ar1 em --------------------------------------------------------------


unlink("../output/em_intermediate/", recursive = T, force = T) # clean up again
for (year_cur in 1985:2020){
  
  runme = paste0(
    'Rscript 2a_ar1_em_est.R'
    , ' --em_type=', 'ar1'
    , ' --output_ts_pred=', TRUE # output the ts prediction for analysis
    , ' --impute_yr=', year_cur
    , ' --impute_vec=', '../output/signals.txt'
    , ' --ar1_sample_length=60'
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




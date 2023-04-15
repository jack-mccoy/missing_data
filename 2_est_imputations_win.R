# Run imputations on windows
#   I guess a bat or powershell script would be a bit nicer but I know how to do this - ac

# The actual imputations are run on each month separately
#   The loop here passes a given year into 2a_em_est.R
#   2a_em_est.R estimates all 12 months in the given year in parallel

for (year_cur in 1985:2020){
  
  runme = paste0(
    'Rscript 2a_em_est.R'
    , ' --impute_yr=', year_cur
    , ' --impute_vec=', '../output/signals.txt'
    , ' --maxiter=10000'
    , ' --out_path=../output/impute_ests/'
    , ' --tol=1e-4'
    , ' --boxcox'
    , ' --cores_frac=0.5'
  )
  
  shell(runme)
  
}
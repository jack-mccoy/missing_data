# Run imputations on windows
#   I guess a bat or powershell script would be a bit nicer but I know how to do this - ac

# 1985 causes errors?

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
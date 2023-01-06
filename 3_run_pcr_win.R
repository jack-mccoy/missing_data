# script for running 3*.R scripts on Windows
# based on 3_run_pcr.sh


# Settings ----------------------------------------------------------------

start_yr=1995
end_yr=2020
main_path="../output/"
tmp_file_imp=paste0(main_path, "imp_tmp.csv")
tmp_file_bc=paste0(main_path, "bc_tmp.csv")

 
# Parameters for principal component regressions
signals_file=paste0(main_path, "signals.txt") # file with list of signals to use
n_pcs=50 # number of PCs in maximal regression
n_years=10 # number of years for principal components/predictive regs
quantile_prob=0.1 # quantile to form long/short portfolios (0.1 means deciles)

# The sample start year is however many years before start_yr (to make PCs)
sample_start_yr=((start_yr-n_years))

# Paths
params_path=paste0(main_path, "impute_ests/")
output_path=paste0(main_path, "pcr_returns/")


# Data prep ---------------------------------------------------------------

# Submit the data prepping job
runme = paste0(
  'Rscript 3a_prep_data_em.R'
  , ' --impute_vec=', signals_file
  , ' --sample_start_year=', sample_start_yr
  , ' --sample_end_year=', end_yr
  , ' --params_path=', params_path
  , ' --tmp_file_bc=', tmp_file_bc
  , ' --tmp_file_imp=', tmp_file_imp
)

shell(runme)


# Run PCRs each year ------------------------------------------------------

for (yr in start_yr:end_yr){
  
  print(paste0('==== Starting PCR for year ', yr))
  tic = Sys.time()
  
  for (mon in 1:12){
    
    print(paste0('  month ', mon))
    
    # simple mean (mn) imputations
    runme = paste0(' Rscript 3b_pcr.R'
    , ' --signals_keep=', signals_file
    , ' --data_file=', tmp_file_bc
    , ' --prefix=', "pcr_mn_"
    , ' --out_path=', output_path
    , ' --iter_year=', yr
    , ' --iter_month=', mon
    , ' --n_yrs=', n_years
    , ' --quantile_prob=', quantile_prob
    , ' --n_pcs=', n_pcs
    )
    
    shell(runme)

    # EM imputations
    runme = paste0(' Rscript 3b_pcr.R'
                   , ' --signals_keep=', signals_file
                   , ' --data_file=', tmp_file_imp
                   , ' --prefix=', "pcr_em_"
                   , ' --out_path=', output_path
                   , ' --iter_year=', yr
                   , ' --iter_month=', mon
                   , ' --n_yrs=', n_years
                   , ' --quantile_prob=', quantile_prob
                   , ' --n_pcs=', n_pcs
    )
    
    shell(runme)
  }
  
  toc = Sys.time()
  
  print(paste0('==== Finished PCR for year ', yr))
  print(paste0('====  minutes required: ', round(as.numeric(toc-tic, units = 'mins'))))
  
}



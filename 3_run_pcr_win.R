# script for running 3*.R scripts on Windows
# based on 3_run_pcr.sh


# Settings ----------------------------------------------------------------

start_yr=1995
end_yr=2020
main_path="../output/"

# Parameters for principal component regressions
signals_file=paste0(main_path, "signals.txt") # file with list of signals to use
n_pcs=50 # number of PCs in maximal regression
n_years=10 # number of years for principal components/predictive regs
quantile_prob=0.1 # quantile to form long/short portfolios (0.1 means deciles)

# The sample start year is however many years before start_yr (to make PCs)
sample_start_yr=((start_yr-n_years))

# Paths
params_path=paste0(main_path, "impute_ests/")
bcsignals_out_path=paste0(main_path, "bcsignals/")
pcr_out_path=paste0(main_path, "pcr_returns/")

dir.create(bcsignals_out_path, showWarnings = F)
dir.create(pcr_out_path, showWarnings = F)

# Big datasets bcsignals filenames
tmp_file_bc=paste0(bcsignals_out_path, "bcsignals_none.csv")
tmp_file_em=paste0(bcsignals_out_path, "bcsignals_em.csv")
tmp_file_ac =paste0(bcsignals_out_path, "bcsignals_availcase.csv")

# fraction of cores to use
cores_frac = 0.5

# switch to using Huang et al scaled pca
scaled_pca = T

# Data prep ---------------------------------------------------------------

# prep BC only data
runme = paste0(
  'Rscript 3a_prep_big_data.R'
  , ' --impute_type=none'
  , ' --impute_vec=', signals_file
  , ' --sample_start_year=', sample_start_yr
  , ' --sample_end_year=', end_yr
  , ' --params_path=', params_path
  , ' --bcsignals_filename=', tmp_file_bc
  , ' --cores_frac=', cores_frac
)

shell(runme)

# prep EM data
runme = paste0(
  'Rscript 3a_prep_big_data.R'
  , ' --impute_type=em'
  , ' --impute_vec=', signals_file
  , ' --sample_start_year=', sample_start_yr
  , ' --sample_end_year=', end_yr
  , ' --params_path=', params_path
  , ' --bcsignals_filename=', tmp_file_em
  , ' --cores_frac=', cores_frac  
)

shell(runme)

# prep availcase data
runme = paste0(
  'Rscript 3a_prep_big_data.R'
  , ' --impute_type=availcase'
  , ' --impute_vec=', signals_file
  , ' --sample_start_year=', sample_start_yr
  , ' --sample_end_year=', end_yr
  , ' --params_path=', params_path
  , ' --bcsignals_filename=', tmp_file_ac
  , ' --cores_frac=', cores_frac  
)

shell(runme)


# Run PCRs each year ------------------------------------------------------

for (yr in start_yr:end_yr){
  
  print(paste0('==== Starting PCR for year ', yr))
  tic = Sys.time()
  
  for (mon in 1:12){
    
    print(paste0('  month ', mon))
    
    # Available case imputations
    runme = paste0(' Rscript 3b_pcr.R'
                   , ' --signals_keep=', signals_file
                   , ' --data_file=', tmp_file_ac
                   , ' --prefix=', "pcr_ac_"
                   , ' --out_path=', pcr_out_path
                   , ' --iter_year=', yr
                   , ' --iter_month=', mon
                   , ' --n_yrs=', n_years
                   , ' --quantile_prob=', quantile_prob
                   , ' --n_pcs=', n_pcs
                   , ' --cores_frac=', cores_frac
                   , ' --scaled_pca=', scaled_pca
    )
    shell(runme)    
    
    # simple mean imputations
    runme = paste0(' Rscript 3b_pcr.R'
                   , ' --signals_keep=', signals_file
                   , ' --data_file=', tmp_file_bc
                   , ' --prefix=', "pcr_mn_"
                   , ' --out_path=', pcr_out_path
                   , ' --iter_year=', yr
                   , ' --iter_month=', mon
                   , ' --n_yrs=', n_years
                   , ' --quantile_prob=', quantile_prob
                   , ' --n_pcs=', n_pcs
                   , ' --cores_frac=', cores_frac    
                   , ' --scaled_pca=', scaled_pca                   
    )
    shell(runme)
    
    # EM imputations
    runme = paste0(' Rscript 3b_pcr.R'
                   , ' --signals_keep=', signals_file
                   , ' --data_file=', tmp_file_em
                   , ' --prefix=', "pcr_em_"
                   , ' --out_path=', pcr_out_path
                   , ' --iter_year=', yr
                   , ' --iter_month=', mon
                   , ' --n_yrs=', n_years
                   , ' --quantile_prob=', quantile_prob
                   , ' --n_pcs=', n_pcs
                   , ' --cores_frac=', cores_frac         
                   , ' --scaled_pca=', scaled_pca                   
    )
    shell(runme)
        
  }
  
  toc = Sys.time()
  
  print(paste0('==== Finished PCR for year ', yr))
  print(paste0('====  minutes required: ', round(as.numeric(toc-tic, units = 'mins'))))
  
}



# script for running 3*.R scripts on Windows
# based on 3_run_pcr.sh


# Settings ----------------------------------------------------------------

start_yr=1995
end_yr=2020
main_path="../output/"

# Parameters for principal component regressions
signals_file=paste0(main_path, "signals_20.txt") # file with list of signals to use
n_pcs=50 # number of PCs in maximal regression
n_years=10 # number of years for principal components/predictive regs
quantile_prob=0.1 # quantile to form long/short portfolios (0.1 means deciles)

# The sample start year is however many years before start_yr (to make PCs)
sample_start_yr=((start_yr-n_years))

# Paths
params_path=paste0(main_path, "impute_ests/")
bcsignals_out_path=paste0(main_path, "bcsignals/")
pc_ret_path = paste0(main_path, 'pc_returns/')

dir.create(bcsignals_out_path, showWarnings = F)
dir.create(pc_ret_path, showWarnings = F)

# Big datasets bcsignals filenames
tmp_file_bc=paste0(bcsignals_out_path, "bcsignals_none.csv")
tmp_file_em=paste0(bcsignals_out_path, "bcsignals_em.csv")
tmp_file_ac =paste0(bcsignals_out_path, "bcsignals_availcase.csv")

# fraction of cores to use
cores_frac = 0.5

# data prep or skip
make_bcsignals = TRUE

# Data prep ---------------------------------------------------------------

if (make_bcsignals){
  
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
  
}

# Make Portfolios ---------------------------------------------------------

## Specification list ------------------------------------------------------

spec_list = expand.grid(
  forecast = c('spca1','spca2', 'pca')
  , imp = c('em','ac','mn')
  , stringsAsFactors = FALSE
)

# spec_list = expand.grid(
#   forecast = c('spca1')
#   , imp = c('em')
#   , stringsAsFactors = FALSE
# ) 

tic_loop = Sys.time()
file.remove(paste0(main_path, 'pcr.log'))
write.table('pcr spec list', paste0(main_path, 'pcr.log')
            , append = TRUE, row.names = FALSE, col.names = FALSE)
write.table(spec_list, paste0(main_path, 'pcr.log')
            , append = TRUE, row.names = FALSE, col.names = FALSE)

## Forecast and form ports each spec-month ------------------------------------------------------
for (yr in start_yr:end_yr){
  
  print(paste0('==== Starting PCR for year ', yr))
  tic = Sys.time()
  
  for (mon in 1:12){
    
    # feedback for sanity
    min_elapsed = round(as.numeric(Sys.time() - tic_loop, units = 'mins'), 1)
    mon_done = (yr*12 + mon) - start_yr*12 -1
    mon_remain = 12*(end_yr - start_yr + 1) - mon_done
    print(paste0('  month ', mon))
    log.text <- paste0(
      Sys.time()
      , " year = ", yr
      , " mon = ", mon
      , " min elapsed = ", min_elapsed
      , " min per mon = ", round(min_elapsed / mon_done, 2)
      , " min remaining = ", round(min_elapsed / mon_done * mon_remain, 1)
    )
    write.table(log.text, paste0(main_path, 'pcr.log')
                , append = TRUE, row.names = FALSE, col.names = FALSE)
    
    for (speci in 1:dim(spec_list)[1]){
      
      cur_spec = spec_list[speci, ]
      
      print(cur_spec)
      
      # pc returns output path
      out_path = paste0(pc_ret_path, cur_spec$forecast, '_', cur_spec$imp, '/')
      dir.create(out_path, showWarnings = FALSE)
      
      # forecast settings
      if (cur_spec$forecast == 'pca'){
        scaled_pca = 'FALSE'
        scaled_pca_weight = 'ew'
      } else if (cur_spec$forecast == 'spca1'){
        scaled_pca = 'TRUE'
        scaled_pca_weight = 'ew'
      } else if (cur_spec$forecast == 'spca2'){
        scaled_pca = 'TRUE'
        scaled_pca_weight = 'vw'
      }
      
      # imputation settings 
      if (cur_spec$imp == 'em'){
        data_file = paste0(main_path, 'bcsignals/', 'bcsignals_em.csv')
      } else if (cur_spec$imp == 'ac'){
        data_file = paste0(main_path, 'bcsignals/', 'bcsignals_availcase.csv')
      } else if (cur_spec$imp == 'mn'){
        data_file = paste0(main_path, 'bcsignals/', 'bcsignals_none.csv')
      }
      
      # forecast and get returns
      runme = paste0(' Rscript 3b_pcr.R'
                     , ' --out_path=', out_path
                     , ' --scaled_pca=', scaled_pca
                     , ' --scaled_pca_weight=', scaled_pca_weight
                     , ' --signals_keep=', signals_file
                     , ' --data_file=', data_file
                     , ' --iter_year=', yr
                     , ' --iter_month=', mon
                     , ' --n_yrs=', n_years
                     , ' --quantile_prob=', quantile_prob
                     , ' --n_pcs=', n_pcs
                     , ' --cores_frac=', cores_frac    
      )
      shell(runme)
      
    }
    
    
  }
  
  toc = Sys.time()
  
  print(paste0('==== Finished PCR for year ', yr))

  
}



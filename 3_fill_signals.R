# make filled permno-yyyymm-[signals] data

# Settings ----------------------------------------------------------------

# ac: perhaps start and end should be automated
sample_start_yr =1985 
end_yr=2020
main_path="../output/"

# Parameters for principal component regressions
signals_file=paste0(main_path, "signals.txt") # file with list of signals to use

# Paths
params_path=paste0(main_path, "impute_ests/")
bcsignals_out_path=paste0(main_path, "bcsignals/")

dir.create(bcsignals_out_path, showWarnings = F)

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


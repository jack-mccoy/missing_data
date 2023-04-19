# reads in many bcsignals_emar*.csv files made by *ar1_em_est.R and 
# makes one big csv and also cleans up.



# options parsing ---------------------------------------------------------
library(optparse) # Python-like option parsing

option_list <- list(
    optparse::make_option(c("--out_path"),
        type = "character", default = "../output/",
        help = "output directory"),
    optparse::make_option(c("--big_file_name"),
        type = "character", default = "bcsignals_em.csv",
        help = "bcsignals_em.csv or bcsignals_emar1.csv"),
    optparse::make_option(c("--ts_file_name"),
        type = "character", default = "ts_prediction_em.csv",
        help = "optional, name of ts prediction file name")  
)


opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)


# bind emar1 annual csvs and save -------------------------------------------------------------------------
library(data.table)

# read in many csv files
csv.list = list.files(
    paste0(opt$out_path, 'em_intermediate/'),
    pattern = 'bcsignals_emar1',
    full.names = T
)

lst <- lapply(csv.list, fread)

# bind together and write
bcsignals_emar1 = rbindlist(lst)

fwrite(bcsignals_emar1, paste0(opt$out_path, 'bcsignals/', opt$big_file_name))

# bind ts prediction annual csvs and save -------------------------------------------------------------------------
# (only if they exist)

# read in many csv files
ts.list = list.files( 
    paste0(opt$out_path, 'em_intermediate/'),
    pattern = 'ts_prediction_',
    full.names = T
)

if (length(ts.list)>0){
  
  lst <- lapply(ts.list, fread)
  
  # bind together and write
  ts_prediction = rbindlist(lst)
  
  fwrite(ts_prediction, paste0(opt$out_path, 'bcsignals/', opt$ts_file_name))
  
} # if (length(ts.lst)>0)


# check ------------------------------------------------------------------
library(tidyverse)

names(bcsignals_emar1) 

signallist = names(bcsignals_emar1)
signallist = signallist[3:length(signallist)]


bcsignals_emar1[
  , lapply(.SD, mean)
  , .SDcols = signallist
]



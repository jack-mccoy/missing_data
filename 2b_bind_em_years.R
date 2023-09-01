# reads in many [--pattern]*.csv files made by *ar1_em_est.R and 
# makes one big csv and also cleans up.

#===============================================================================
# Libraries
#===============================================================================

library(data.table)
library(optparse) # Python-like option parsing
source("functions.R") # for file path function

#===============================================================================
# options parsing 
#===============================================================================

# Set up options
option_list <- list(
    #optparse::make_option(c("--data_path"),
    #    type = "character", default = "../output/",
    #    help = "data directory, assumed to have subdir `em_intermediate/`"),
    optparse::make_option(c("--em_type"),
        type = "character", default = "regular",
        help = "type of EM imputation to bind, one of (regular, ar1)"),
    optparse::make_option(c("--ts_file_name"),
        type = "character", default = "ts_prediction_em.csv",
        help = "optional, name of ts prediction file name")  
)

# Unpack options
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Input checking plus creating file names and grep patterns
if (opt$em_type == "regular") {
    big_file_name <- "bcsignals_em.csv"
    prefix <- "bcsignals_em_"
} else if (opt$em_type == "ar1") {
    big_file_name <- "bcsignals_emar1.csv"
    prefix <- "bcsignals_emar1_"
} else {
    stop("Option `--em_type` must be one of (regular, ar1)\n")
}

# File paths
getFilePaths()

#===============================================================================
# bind emar1 annual csvs and save 
#===============================================================================

# read in many csv files
csv_list <- list.files(
    paste0(FILEPATHS$data_path, 'em_intermediate/'),
    pattern = prefix,
    full.names = TRUE
)

# Read all CSVs to a list
lst <- lapply(csv_list, fread)

# bind together and write
bcsignals_imp <- rbindlist(lst)
fwrite(bcsignals_imp, paste0(FILEPATHS$data_path, 'bcsignals/', big_file_name))

#===============================================================================
# bind ts prediction annual csvs and save (only if they exist)
#===============================================================================

# read in many csv files
ts.list <- list.files( 
    paste0(FILEPATHS$out_path, 'em_intermediate/'),
    pattern = 'ts_prediction_',
    full.names = T
)

if (length(ts.list)>0){
  
  lst <- lapply(ts.list, fread)
  
  # bind together and write
  ts_prediction = rbindlist(lst)
  fwrite(ts_prediction, paste0(FILEPATHS$out_path, 'bcsignals/', opt$ts_file_name))
  
}

#===============================================================================
# check to console
#===============================================================================

names(bcsignals_imp) 

signallist <- names(bcsignals_imp)[3:ncol(bcsignals_imp)]
signallist <- signallist[3:length(signallist)]

bcsignals_imp[ , lapply(.SD, mean), .SDcols = signallist]



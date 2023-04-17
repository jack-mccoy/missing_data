# reads in many bcsignals_emar*.csv files made by *ar1_em_est.R and 
# makes one big csv and also cleans up.



# options parsing ---------------------------------------------------------
library(optparse) # Python-like option parsing

option_list <- list(
  optparse::make_option(c("--big_file_name"),
                        type = "character", default = "bcsignals_em.csv",
                        help = "bcsignals_em.csv or bcsignals_emar1.csv")
)


opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)


# do stuff -------------------------------------------------------------------------

library(data.table)


# read in many csv files
csv.list = list.files(
  '../output/em_intermediate/', pattern = 'bcsignals_emar1'
  , full.names = T
)

lst <- lapply(csv.list, fread)

# bind together and write
bcsignals_emar1 = rbindlist(lst)

fwrite(bcsignals_emar1
       , paste0('../output/bcsignals/', opt$big_file_name)
       )


# check ------------------------------------------------------------------
library(tidyverse)

names(bcsignals_emar1) 

signallist = names(bcsignals_emar1)
signallist = signallist[3:length(signallist)]


bcsignals_emar1[
  , lapply(.SD, mean)
  , .SDcols = signallist
]


# -------------------------------------------------------------------------


temp = signallist[6]

temp

bcsignals_emar1 %>% 
  select(permno,yyyymm,all_of(temp)) %>% 
  rename(x := !!temp) %>% 
  filter(is.na(x)) %>% 
  head()

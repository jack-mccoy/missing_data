# (also downloads Chen-Zimmermann signals)
#==============================================================================#
# Packages and Setup ====
#==============================================================================#

library(data.table)
library(getPass)
library(RPostgres)
library(zoo)
library(googledrive)
library(tidyverse)

dir.create('../data/', showWarnings = F)
dir.create('../output/', showWarnings = F)


#==============================================================================#
# Download CRSP ====
#==============================================================================#

## Download size and returns ----

wrds_user <- getPass::getPass("WRDS username: ")
wrds_pass <- getPass::getPass("WRDS pass: ")

wrds <- DBI::dbConnect(RPostgres::Postgres(),
    host = "wrds-pgdata.wharton.upenn.edu",
    db = "wrds",
    port = 9737,
    user = wrds_user, 
    pass = wrds_pass)

crspm <- as.data.table(DBI::dbGetQuery(wrds, "
    SELECT a.permno, a.date, a.ret, a.shrout, a.prc, 
        b.exchcd,
        c.dlstcd, c.dlret   -- from delistings table
    FROM crsp.msf AS a
    LEFT JOIN crsp.msenames AS b
        ON a.permno=b.permno
        AND b.namedt<=a.date
        AND a.date<=b.nameendt
    LEFT JOIN crsp.msedelist AS c
        ON a.permno=c.permno
        AND date_trunc('month', a.date) = date_trunc('month', c.dlstdt)
;"))

DBI::dbDisconnect(wrds)

## Delisting returns ----

crspm[, # Adapted from Andrew's github code
      dlret := dplyr::case_when(
        is.na(dlret) 
        & (dlstcd == 500 | (520 <= dlstcd & dlstcd <= 584))
        & (exchcd == 1 | exchcd == 2) 
        ~ -0.25,
        is.na(dlret)
        & (dlstcd == 500 | (dlstcd >=520 & dlstcd <=584))
        & (exchcd == 3)
        ~ -0.55,
        dlret < -1 & !is.na(dlret) ~ -1,
        TRUE ~ 0
      )
][ # return is the return plus delisting return
  , ret := ret + dlret
][ # in cases where no ret but non-zero dlret
  is.na(ret) & dlret != 0,
  ret := dlret
]

## Lastly, the signals that aren't in Andrew's master file ----
# Signals are SIGNED (higher signal ==> higher return)

crsp_final <- crspm[, .(
  permno,
  yyyymm = as.yearmon(as.Date(date)),
  ret = 100 * ret,
  me = abs(prc) * shrout,
  price = -1 * log(abs(prc))
)][, ":="(
  streversal = -1 * ifelse(is.na(ret), 0, ret),
  size = ifelse(me > 0, -1 * log(me), NA)
)]

# Memory
rm(crspm)


## Output ----
fwrite(crsp_final, "../data/crsp_data.csv")


#==============================================================================#
# Download Chen-Zimmermann Data ====
#==============================================================================#

# root of March 2022 release on Gdrive
pathRelease = 'https://drive.google.com/drive/folders/1O18scg9iBTiBaDiQFhoGxdn4FdsbMqGo'


## dl signal documentation ====
target_dribble = pathRelease %>% drive_ls() %>% 
  filter(name=='SignalDoc.csv')

drive_download(target_dribble, path = '../data/SignalDoc.csv', overwrite = T)

## dl signals (except the crsp ones) ====

# 2 gig dl, can take a few minutes
# download
target_dribble = pathRelease %>% drive_ls() %>% 
  filter(name == 'Firm Level Characteristics') %>% drive_ls() %>% 
  filter(name == 'Full Sets') %>% drive_ls() %>% 
  filter(name == 'signed_predictors_dl_wide.zip') 
dl = drive_download(target_dribble, path = '../data/deleteme.zip', overwrite = T)

# unzip, clean up
unzip('../data/deleteme.zip', exdir = '../data')
file.remove('../data/deleteme.zip')


## download ports ====
# this is just for counting and making signal lists
target_dribble = pathRelease %>% drive_ls() %>% 
  filter(name == 'Portfolios') %>% drive_ls() %>% 
  filter(name == 'Full Sets OP') %>% drive_ls() %>% 
  filter(name == 'PredictorPortsFull.csv')

drive_download(target_dribble, path = '../data/PredictorPortsFull.csv', overwrite = T)


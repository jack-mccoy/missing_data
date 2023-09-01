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
source("functions.R")

# Gets the file paths and creates directories as needed
getFilePaths(signal_list_exists = FALSE) # signal list not yet made

# Create directory for downloaded data that has not been transformed or imputed
dir.create(paste0(FILEPATHS$data_path, "raw/"), showWarnings = FALSE)

# log in ----
wrds_user <- getPass::getPass("WRDS username: ")
wrds_pass <- getPass::getPass("WRDS pass: ")

wrds <- DBI::dbConnect(RPostgres::Postgres(),
    host = "wrds-pgdata.wharton.upenn.edu",
    db = "wrds",
    port = 9737,
    user = wrds_user, 
    pass = wrds_pass)

#==============================================================================#
# Download CRSP ====
#==============================================================================#

## Download size and returns ----

crspm <- as.data.table(DBI::dbGetQuery(wrds, "
    SELECT a.permno, a.date, a.ret, a.shrout, a.prc, a.hsiccd, 
        b.exchcd, b.shrcd,
        c.dlstcd, c.dlret   -- from delistings table
    FROM crsp.msf AS a
    LEFT JOIN crsp.msenames AS b
        ON a.permno=b.permno
        AND b.namedt<=a.date
        AND a.date<=b.nameendt
    LEFT JOIN crsp.msedelist AS c
        ON a.permno=c.permno
        AND date_trunc('month', a.date) = date_trunc('month', c.dlstdt)
    WHERE /* Standard share and exchange code filter */
        b.shrcd IN (10, 11, 12) AND
        b.exchcd IN (1, 2, 3)
;"))

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

# Lastly, the signals that aren't in Andrew's master file ----
# Signals are SIGNED (higher signal ==> higher return)

crsp_final <- crspm[, .(
    permno,
    hsiccd,
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

# Output ----
fwrite(crsp_final, paste0(FILEPATHS$data_path, "raw/crsp_data.csv"))

#==============================================================================#
# Fama-French Factors ====
#==============================================================================#

# Fama-French factors
ff5_mom <- as.data.table(dbGetQuery(wrds, "
    SELECT date,
        mktrf * 100 as mktrf, /* The PCR returns are in pct out of 100 */
        smb * 100 as smb,
        hml * 100 as hml,
        rmw * 100 as rmw,
        cma * 100 as cma,
        umd * 100 as umd
    FROM ff_all.fivefactors_monthly
    ORDER BY date
;"))

DBI::dbDisconnect(wrds)

fwrite(ff5_mom, paste0(FILEPATHS$data_path, "raw/ff5_factors.csv"))

rm(ff5_mom)

#==============================================================================#
# Download Chen-Zimmermann Data ====
#==============================================================================#

# root of Aug 2023 release on Gdrive
path_release <- "https://drive.google.com/drive/folders/1EP6oEabyZRamveGNyzYU0u6qJ-N43Qfq"

# Allow non-interactive download
drive_deauth()

# dl signal documentation ====

target_dribble = path_release %>% drive_ls() %>% 
  filter(name=='SignalDoc.csv')

drive_download(target_dribble, 
    path = paste0(FILEPATHS$data_path, 'raw/SignalDoc.csv'),
    overwrite = TRUE)

# dl signals (except the crsp ones) ====

# 2 gig dl, can take a few minutes
# download
target_dribble <- path_release %>% drive_ls() %>% 
    filter(name == 'Firm Level Characteristics') %>% drive_ls() %>% 
    filter(name == 'Full Sets') %>% drive_ls() %>% 
    filter(name == 'signed_predictors_dl_wide.zip') 
dl <- drive_download(target_dribble, 
    path = paste0(FILEPATHS$data_path, 'raw/deleteme.zip'),
    overwrite = TRUE)

# unzip, clean up
unzip(paste0(FILEPATHS$data_path, 'raw/deleteme.zip'), 
    exdir = paste0(FILEPATHS$data_path, 'raw/'))
file.remove(paste0(FILEPATHS$data_path, 'raw/deleteme.zip'))

# download ports ====
# this is just for counting and making signal lists
target_dribble = path_release %>% drive_ls() %>% 
    filter(name == 'Portfolios') %>% drive_ls() %>% 
    filter(name == 'Full Sets OP') %>% drive_ls() %>% 
    filter(name == 'PredictorPortsFull.csv')

drive_download(target_dribble, 
    path = paste0(FILEPATHS$data_path, 'raw/PredictorPortsFull.csv'),
    overwrite = TRUE)

# memory
rm(target_dribble, dl)

#==============================================================================#
# Filtering ====
#==============================================================================#

# Maximize memory
gc()

signals <- fread(paste0(FILEPATHS$data_path, "raw/signed_predictors_dl_wide.csv"))

# For merging
if (class(signals$yyyymm) != "yearmon") {
    signals[, yyyymm := as.yearmon(as.character(yyyymm * 10 + 1), format = "%Y%m%d")]
}

# merge in (the merge is the filter because it's an inner join by default)
cat("Before share code filter, there are", nrow(signals), "observations",
    "in the signals data set")

signals_filt <- merge(signals, crsp_final[, .(permno, yyyymm)], 
    by = c("permno", "yyyymm")) 

cat("After share code filter, there are", nrow(signals_filt), "observations",
    "in the signals data set")

fwrite(signals_filt, 
    paste0(FILEPATHS$data_path, "raw/signed_predictors_dl_wide_filtered.csv"))


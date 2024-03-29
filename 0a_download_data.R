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

# log in if .pgpass not set up (better to have it set up) ----
if (file.exists("~/.pgpass")) {
    wrds_user <- NULL
    wrds_pass <- NULL
} else {
    wrds_user <- getPass::getPass("WRDS username: ")
    wrds_pass <- getPass::getPass("WRDS pass: ")
}

# Set up connection
wrds <- DBI::dbConnect(RPostgres::Postgres(),
    host = "wrds-pgdata.wharton.upenn.edu",
    db = "wrds",
    port = 9737,
    user = wrds_user, 
    pass = wrds_pass)

#==============================================================================#
# Size breakpoints from Ken French
#==============================================================================#

# Website path
size_path <- "https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/ME_Breakpoints_CSV.zip"

# Download and unzip
download.file(size_path, paste0(FILEPATHS$data_path, 'raw/deleteme.zip'))
decompress_file(directory = paste0(FILEPATHS$data_path, "raw/"),
    file = "deleteme.zip")

# Read in data and clean ====

size_quantiles <- fread(paste0(FILEPATHS$data_path, 'raw/ME_Breakpoints.CSV'),
    skip = 1,
    col.names = c("yyyymm", "n", paste0("p", seq(5, 100, by = 5))))

# Match the other yyyymm formats
size_quantiles[, yyyymm := as.yearmon(as.character(yyyymm*10+1), format="%Y%m%d")]

# Output just the ones we need
fwrite(size_quantiles[, .SD, .SDcols = c("yyyymm", "p20", "p50")],
    paste0(FILEPATHS$data_path, "raw/ME_Breakpoints_20_50_clean.csv"))

# memory
rm(size_quantiles)

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
        b.shrcd IN (10, 11, 12) AND /* 10, 11, 12: ordinary common shares */
        b.exchcd IN (1, 2, 3, 31, 32, 33) /* 1: NYSE, 2: AMEX, 3: Nasdaq; 3X are 'when-issued' equivalents */ 
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
decompress_file(directory = paste0(FILEPATHS$data_path, "raw/"),
    file = "deleteme.zip")

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
    "in the signals data set\n")

signals_filt <- merge(signals, crsp_final[, .(permno, yyyymm)], 
    by = c("permno", "yyyymm")) 

cat("After share code filter, there are", nrow(signals_filt), "observations",
    "in the signals data set\n")

fwrite(signals_filt, 
    paste0(FILEPATHS$data_path, "raw/signed_predictors_dl_wide_filtered.csv"))


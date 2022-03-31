
#===============================================================================
# Packages
#===============================================================================

library(data.table)
library(getPass)
library(RPostgres)
library(zoo)

#===============================================================================
# Download
#===============================================================================

# Download size and returns ----

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

# Delisting returns ----

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

#===============================================================================
# Output
#===============================================================================

fwrite(crsp_final, "crsp_data.csv")


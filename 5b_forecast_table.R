
library(data.table)
library(kableExtra)
library(stringr)
library(zoo)

source("functions.R")

#===============================================================================
# Paths
#===============================================================================

# Load the projects file paths
getFilePaths()

# Set up the file list
fcast_dir <- paste0(FILEPATHS$out_path, "forecast/")
fcasts <- system(paste("ls", fcast_dir), intern=TRUE) 
fcasts_files <- paste0(fcast_dir, fcasts, "/permno-month-forecast.csv")

#===============================================================================
# Read in data
#===============================================================================

# Read in data
crsp_data <- fread(paste0(FILEPATHS$data_path, "raw/crsp_data.csv"))[, c("permno", "yyyymm", "me")]
fcast_data <- rbindlist(lapply(fcasts, function(f) { 
    file <- paste0(fcast_dir, f, "/permno-month-forecast.csv")
    method <- stringr::str_split(f, "-")[[1]][1]
    imp <- stringr::str_split(f, "_")[[1]][2]
    dt <- tryCatch(fread(file)[, ":="(method = method, imp = imp)], 
        error = function(e) return(NULL))
    return(dt)
}))

# Merge and make dates nice 
fcast_data <- merge(fcast_data, crsp_data, by = c("permno", "yyyymm"))

#===============================================================================
# Create long-short portfolios
#===============================================================================

# Define deciles
fcast_data[,
    ":="(
        pctile10 = quantile(Ebh1m, 0.1, na.rm = TRUE),
        pctile90 = quantile(Ebh1m, 0.9, na.rm = TRUE)
    ),
    by = .(yyyymm, method, imp)
][,
    decile := ifelse(Ebh1m <= pctile10, 1, ifelse(Ebh1m >= pctile90, 10, NA)) 
]

# Average portfolios by prediction decile for month, method, and imp
# make negative (short) if lowest decile
ports <- fcast_data[
    decile %in% c(1, 10), 
    .(
        ew = mean(ifelse(decile == 1, -1, 1) * bh1m, na.rm = TRUE),
        vw = weighted.mean(ifelse(decile == 1, -1, 1) * bh1m, w = me, na.rm = T)
    ),
    by = .(yyyymm, method, imp, decile)
]

# Mean portfolios for each leg
ports_mean <- ports[, .(ew = mean(ew), vw = mean(vw)), by = .(method, imp, decile)]

# Long-short portfolios (annualized)
ls_ports_mean <- ports_mean[, 
    .(ew = 12*sum(ew), vw = 12*sum(vw)),
    by = .(method, imp)
]

# Convert to nicer table with imputation rows and fcast method cols
table_ew <- dcast(ls_ports_mean, "imp ~ method", value.var = "ew")
table_vw <- dcast(ls_ports_mean, "imp ~ method", value.var = "vw")

#===============================================================================
# Output
#===============================================================================

# In case it doesn't already exist
dir.create(paste0(FILEPATHS$out_path, "tables"))

fwrite(table_ew, paste0(FILEPATHS$out_path, 'tables/fcast_table_ew.csv'))
fwrite(table_vw, paste0(FILEPATHS$out_path, 'tables/fcast_table_vw.csv'))



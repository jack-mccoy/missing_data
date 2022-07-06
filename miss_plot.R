
#===============================================================================
# Packages and functions
#===============================================================================

library(data.table)
library(ggplot2)
library(optparse)
library(zoo)

source("functions.R")

#===============================================================================
# Hardcodes
#===============================================================================

plot_path <- "../output_best100_1985/plots/"

dates <- as.yearmon(c("Jun 1990", "Jun 2000", "Jun 2010"))

#===============================================================================
# Signals
#===============================================================================

signals_keep <- scan("../data/signals_best100_1985.txt", character()) 

#===============================================================================
# Read in the data
#===============================================================================

dat <- fread("../output_best100_1985/bc_tmp.csv", select = c("permno", "yyyymm", signals_keep))[
    as.yearmon(yyyymm) %in% dates # Keep only the months we need
]

dat[, yyyymm := as.yearmon(yyyymm)]

#===============================================================================
# Output plots
#===============================================================================

setwd(plot_path)

for (dd in 1:length(dates)) {
    date_file <- gsub(" ", "", dates[dd]) 
    tmp <- missPlot(dat[yyyymm == dates[dd]], 
        rhs_vars = signals_keep,
        title = "")
    ggsave(plot = tmp, filename = paste0("missplot_", date_file, ".pdf"), 
        width = 7, height = 5, unit = "in")
    rm(tmp)
}


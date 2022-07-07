
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

plot_path <- "../output/plots/"

dates <- as.yearmon(c("Jun 1990", "Jun 2000", "Jun 2010"))


#===============================================================================
# Read in the data
#===============================================================================

dat <- fread("../output/bc_tmp.csv")[
  as.yearmon(yyyymm) %in% dates # Keep only the months we need
]


dat[, yyyymm := as.yearmon(yyyymm)]

#===============================================================================
# Output plots
#===============================================================================

# setwd(plot_path)

signals_keep = names(dat)
signals_keep = signals_keep[!signals_keep %in% c('permno','yyyymm')]

for (dd in 1:length(dates)) {
    date_file <- gsub(" ", "", dates[dd]) 
    tmp <- missPlot(dat[yyyymm == dates[dd]], 
        rhs_vars = signals_keep,
        title = "")
    ggsave(plot = tmp, filename = paste0(plot_path, "missplot_", date_file, ".pdf"), 
        width = 10, height = 7, unit = "in")
    rm(tmp)
}


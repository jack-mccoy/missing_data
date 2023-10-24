#===============================================================================#
# Setup ----
#===============================================================================#

library(data.table)
library(ggplot2)
library(optparse)
library(zoo)
library(tidyverse)
source("functions.R")

# set file paths and choose yearmons
getFilePaths()
plot_path <- paste0(FILEPATHS$out_path, "plots/")
dates <- as.yearmon(c("Jun 1985","Jun 1990", "Jun 2000", "Jun 2010"))

dir.create(plot_path, showWarnings = F)

#===============================================================================#
# Read in the data ----
#===============================================================================#

dat <- fread(paste0(FILEPATHS$data_path, "bcsignals/bcsignals_none.csv"))[
  as.yearmon(yyyymm) %in% dates # Keep only the months we need
]
dat[, yyyymm := as.yearmon(yyyymm)]

#===============================================================================#
# Output plots ----
#===============================================================================#

# Missingness plot function
missPlot <- function(data, rhs_vars, xlab = "Stock i.d.", ylab = "Predictor i.d.",
                     title = "Missingness Map"
) {
  data2 <- copy(as.data.table(data))
  
  nas <- is.na(data2[,..rhs_vars])
  
  sort_by_na <- do.call("order", as.data.table(is.na(data2)))
  
  nas <- as.data.table(nas)[sort_by_na]
  
  col_miss <- colSums(nas)
  col_ind <- order(-col_miss) 
  nas <- nas[, ..col_ind]
  signal_id <- data.table(variable = colnames(nas), 
                          id = factor(ncol(nas):1, levels = ncol(nas):1)) 
  nas[, firm_id := 1:nrow(nas)]
  
  ratio <- nrow(nas)*1/(ncol(nas) - 1)
  nas_long <- merge(melt(nas, id = "firm_id"), signal_id, by = "variable") 
  
  miss_plot <- ggplot(data = nas_long, aes(x = firm_id, y = id)) + 
    geom_raster(aes(fill = value)) +
    labs(x = xlab, y = ylab, color = "") + 
    scale_fill_manual(
      name = "",
      labels = c("Observed","Missing"),
      values = c('lightgrey', MATRED)
    ) + 
    theme_bw() +
    theme(
      text = element_text(size = 22),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      legend.position = "top"
    ) +
    scale_x_continuous(expand = c(0,0), breaks = c(1,seq(1000,12000,1000))) +
    scale_y_discrete(breaks = c(1, seq(0, ncol(nas), by = 25))) +
    ggtitle(title)# + 
  #coord_fixed(ratio)
  
  return(miss_plot)
}

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

#===============================================================================#
# Testing plot by data category ----
#===============================================================================#

# import get data categories
signaldoc = fread(paste0(FILEPATHS$data_path, "raw/SignalDoc.csv "))
signaldoc[ , signalname := tolower(Acronym)]
signaldoc = signaldoc[ , .(signalname, Cat.Data, Cat.Economic)]

# group together minor categories and make factor
data_list = c('Price','Accounting','Analyst','Trading','Other')
signaldoc[ !Cat.Data %in% data_list, Cat.Data := 'Other']
signaldoc[ , Cat.Data 
  := factor(Cat.Data, levels = c('Price','Accounting','Analyst','Trading','Other'))]

# function inputs
signals_keep = names(dat)
signals_keep = signals_keep[!signals_keep %in% c('permno','yyyymm')]

data = dat[yyyymm == dates[1]]
rhs_vars = signals_keep
title = ''
xlab = 'Stock i.d.'
ylab = 'Predictor i.d.'

# function body sketch
data2 <- copy(as.data.table(data))
nas <- is.na(data2[,..rhs_vars])

# sort stocks by number of na's
sort_by_na <- do.call("order", as.data.table(is.na(data2)))
nas <- as.data.table(nas)[sort_by_na]

# sort predictors by data category and then by number of na's
signalsort = merge(
  data.table(signalname = colnames(nas), nas = colSums(nas)), 
  signaldoc, by = 'signalname', all.x = T)
signalsort[ , id_old := .I]
setorder(signalsort, -Cat.Data, -nas)
signalsort[ , id := factor(dim(signalsort)[1]:1, 
  levels = dim(signalsort)[1]:1)]

nas = setcolorder(nas, signalsort$signalname)
# signal_id <- data.table(variable = colnames(nas), 
#                         id = factor(ncol(nas):1, levels = ncol(nas):1)) 
nas[, firm_id := 1:nrow(nas)]

# prep for plotting
nas_long <- merge(melt(nas, id = "firm_id"), signalsort[, .(signalname,Cat.Data,id)], 
  by.x = "variable", by.y = "signalname") 
y_data_change = nas_long %>% group_by(Cat.Data) %>% 
  summarize(nsignal = n_distinct(id)) %>% 
  mutate(y_pos = length(signals_keep) - cumsum(nsignal)+0.5) %>% 
  pull(y_pos)

# plot
color_datacat = 'black'
miss_plot <- ggplot(data = nas_long, aes(x = firm_id, y = id)) + 
  geom_raster(aes(fill = value)) +
  labs(x = xlab, y = ylab, color = "") + 
  scale_fill_manual(
    name = "",
    labels = c("Observed","Missing"),
    values = c('lightgrey', MATRED)
  ) + 
  theme_bw() +
  theme(
    text = element_text(size = 22),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.position = "top"
  ) +
  scale_x_continuous(expand = c(0,0), breaks = c(1,seq(1000,12000,1000))) +
  scale_y_discrete(breaks = c(1, seq(0, ncol(nas), by = 25))) +
  ggtitle(title) +
  # horizontal lines for y_data_change
  geom_hline(yintercept = y_data_change, color = color_datacat) +
  annotate("text", x = 100, y = c(160, y_data_change[1:4])-5.5, 
           label = data_list[seq(1,5,1)], 
           color = color_datacat, size = 5, hjust = 0, vjust = 0) 

ggsave(miss_plot, 
    filename = paste0(FILEPATHS$out_path, "plots/missplotcat_", data2$yyyymm[1], ".pdf"),
    width = 10, height = 7, unit = "in")

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
dates <- as.yearmon(c("Jun 1985","Jun 1990", "Jun 2000", "Jun 2010", "Jun 2020"))

dir.create(plot_path, showWarnings = F)

#===============================================================================#
# Read in the data ----
#===============================================================================#

dat <- fread(paste0(FILEPATHS$data_path, "bcsignals/bcsignals_none.csv"))[
  as.yearmon(yyyymm) %in% dates # Keep only the months we need
]
dat[, yyyymm := as.yearmon(yyyymm)]

# import get data categories
signaldoc = fread(paste0(FILEPATHS$data_path, "raw/SignalDoc.csv "))
signaldoc[ , signalname := tolower(Acronym)]
signaldoc = signaldoc[ , .(signalname, Cat.Data, Cat.Economic)]

# group together minor categories and make factor
data_list = c('Price','Accounting','Analyst','Trading','Other')
signaldoc[ !Cat.Data %in% data_list, Cat.Data := 'Other']
signaldoc[ , Cat.Data 
  := factor(Cat.Data, levels = c('Price','Accounting','Analyst','Trading','Other'))]

# other setup
signals_keep = names(dat)
signals_keep = signals_keep[!signals_keep %in% c('permno','yyyymm')]
nsignal = length(signals_keep)

#==============================================================================#
# Plotting Function ----
#===============================================================================#

# function for outputting missplot pdf
make_missplot_pdf = function(yearm_select){    

  # create (permno,signalname,miss) long data
  misslong = dat[yyyymm == as.yearmon(yearm_select)] %>% select(-yyyymm)  %>% 
    pivot_longer(cols = -permno, names_to = 'signalname', values_to = 'value') %>% 
    mutate(miss = is.na(value)) %>% 
    left_join(signaldoc, by = 'signalname') %>% 
    select(permno, signalname, miss, Cat.Data)     

  # sort signals by data category and then by number of na's
  signalsum = misslong %>% 
    group_by(Cat.Data, signalname) %>% summarize(nmiss_signal = sum(miss)) %>% ungroup() %>% 
    arrange(desc(Cat.Data),-nmiss_signal) %>% 
    # mutate(signal_id = as.factor(row_number(), levels = nsignal:1)) 
    mutate(signal_id = factor(nsignal-row_number()+1, levels = nsignal:1)) 

  # sort stocks by data category and then by number of na's
  stocksum = misslong %>%
    group_by(permno) %>% summarize(nmiss_stock = sum(miss)) %>% ungroup() %>%
    arrange(nmiss_stock) %>%
    mutate(stock_id = row_number())
      
  # merge on sorted ids
  misslong2 = misslong %>% 
    left_join(signalsum %>% select(signalname,signal_id), by = 'signalname') %>%
    left_join(stocksum %>% select(permno,stock_id), by = 'permno') 

  # find signal ids for when data category changes
  y_data_change = misslong2 %>% 
    mutate(nsignal = n_distinct(signal_id)) %>% 
    group_by(Cat.Data) %>% 
    summarize(nsignalcat = n_distinct(signal_id)) %>% 
    mutate(y_pos = nsignal - cumsum(nsignalcat)+0.2)  %>% 
    pull(y_pos)

  # plot
  color_datacat = 'black'
  miss_plot <- ggplot(data = misslong2, aes(x = stock_id, y = signal_id)) + 
    geom_raster(aes(fill = miss)) +
    labs(x = 'Stock i.d.', y = 'Predictor i.d.', color = "") + 
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
      legend.position = "top",
      plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")
    ) +
    scale_x_continuous(expand = c(0,0), breaks = c(1,seq(1000,12000,1000))) +
    scale_y_discrete(breaks = c(1, seq(0, nsignal, by = 25))) +
    # horizontal lines for y_data_change
    geom_hline(yintercept = y_data_change, color = color_datacat) +
    annotate("text", x = 100, y = c(160, y_data_change[1:4])-5.5, 
            label = data_list[seq(1,5,1)], 
            color = color_datacat, size = 5, hjust = 0, vjust = 0) 

  ggsave(miss_plot, 
      filename = paste0(FILEPATHS$out_path, "plots/missplotcat_", yearm_select, ".pdf"),
      width = 10, height = 7, unit = "in")

} # end function make_missplot_pdf

#==============================================================================#
# Run Function Several Times ----
#===============================================================================#

# make plot for each date
for (dd in 1:length(dates)) {    
    make_missplot_pdf(dates[dd])
}

print('Finished making missplots by data category')
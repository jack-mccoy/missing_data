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
dates <- as.yearmon(paste0('Jun ', seq(1985,2020,5)))

dir.create(plot_path, showWarnings = F)

#===============================================================================#
# Read in the data ----
#===============================================================================#

dat <- fread(paste0(FILEPATHS$data_path, "bcsignals/bcsignals_none.csv"))[
  , yyyymm := as.yearmon(yyyymm)
]

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

  # subset data
  permno_select = dat[yyyymm == as.yearmon(yearm_select)]$permno
  datsmall = dat[yyyymm <= as.yearmon(yearm_select) & permno %in% permno_select, ]  %>% 
    pivot_longer(cols = -c(permno, yyyymm), names_to = 'signalname', values_to = 'value') %>% 
    mutate(obs = !is.na(value)) %>% select(permno,signalname,yyyymm,obs) %>% setDT()  
  
  # find most recent observation
  setorder(datsmall, permno, signalname, -yyyymm)
  lastobs = datsmall[obs==TRUE, .SD[1], by = .(permno, signalname)]
  setnames(lastobs, old = 'yyyymm', new = 'lastobs')
  lastobs[ , obs := NULL]

  # merge current slice of data with lastobs
  datselect = merge(datsmall[yyyymm == as.yearmon(yearm_select)], lastobs,
                    by = c('permno','signalname'), all.x = T) %>% 
                    mutate(miss_months = as.numeric(as.yearmon(yearm_select) - lastobs)*12)

  # merge on data categories
  datselect = merge(datselect, signaldoc[ , .(signalname, Cat.Data)], by = 'signalname')

  # sort signals by data category and then by number of na's
  signalsum = datselect %>% 
    group_by(Cat.Data, signalname) %>% summarize(nmiss_signal = sum(!obs)) %>% ungroup() %>% 
    arrange(desc(Cat.Data),-nmiss_signal) %>% 
    # mutate(signal_id = as.factor(row_number(), levels = nsignal:1)) 
    mutate(signal_id = factor(nsignal-row_number()+1, levels = nsignal:1)) 

  # sort stocks by data category and then by number of na's
  stocksum = datselect %>%
    group_by(permno) %>% summarize(nmiss_stock = sum(!obs)) %>% ungroup() %>%
    arrange(nmiss_stock) %>%
    mutate(stock_id = row_number())
      
  # merge on sorted ids
  datselect2 = datselect %>% 
    left_join(signalsum %>% select(signalname,signal_id), by = 'signalname') %>%
    left_join(stocksum %>% select(permno,stock_id), by = 'permno') 

  # find signal ids for when data category changes
  y_data_change = datselect2 %>% 
    mutate(nsignal = n_distinct(signal_id)) %>% 
    group_by(Cat.Data) %>% 
    summarize(nsignalcat = n_distinct(signal_id)) %>% 
    mutate(y_pos = nsignal - cumsum(nsignalcat)+0.2)  %>% 
    pull(y_pos) 

  # define miss months groups
  datselect2 = datselect2 %>% mutate(
    miss_group = case_when(
      miss_months == 0 ~ 'Obs Now',
      miss_months <= 12 ~ '1-12',
      miss_months <= 24 ~ '13-24',
      miss_months > 24 ~ '>24',
      is.na(miss_months) ~ 'Never'
    )
  )

  # plot
  # color reference: http://sape.inf.usi.ch/quick-reference/ggplot2/colour
  color_datacat = MATBLUE # annotation color
  miss_plot = ggplot(
    datselect2 %>% 
    mutate(miss_group = factor(
      miss_group, levels = c('Obs Now','1-12','13-24','>24','Never')
    ))
    , aes(x=stock_id, y=signal_id)) + 
    geom_raster(aes(fill=miss_group))+
    labs(x = 'Stock i.d.', y = 'Predictor i.d.', color = "") +
    scale_fill_manual(
      name = "         Months Since Last Observation",
      labels = c("Obs Now","1-12","13-24",">24","Never Obs"),
      values = c('gray90', 'lightgreen','orange','darkred','gray25')
    ) + 
    guides(fill = guide_legend(title.position = "top", nrow = 1)) +
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
      filename = paste0(FILEPATHS$out_path, "plots/missplotcat_",
                        str_replace(yearm_select,' ','_'), ".pdf"),
      width = 10, height = 8, unit = "in")

  # check share that are never observed 
  # temp = datselect2 %>% group_by(miss_group) %>% summarize(n = n()) 
  # temp %>% mutate(frac = n/sum(temp$n[1:4]))

} # end function make_missplot_pdf

#==============================================================================#
# Run Function Several Times ----
#===============================================================================#

# make plot for each date
for (dd in 1:length(dates)) {    
    make_missplot_pdf(dates[dd])
}

print('Finished making missplots by data category')
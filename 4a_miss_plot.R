#===============================================================================#
# Setup ----
#===============================================================================#

library(data.table)
library(ggplot2)
library(optparse)
library(zoo)

source("functions.R")

# Hardcodes 
plot_path <- "../output/plots/"
dates <- as.yearmon(c("Jun 1985","Jun 1990", "Jun 2000", "Jun 2010"))


#===============================================================================#
# Read in the data ----
#===============================================================================#

dat <- fread("../output/bc_tmp.csv")[
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



#===============================================================================#
# Packages ----
#===============================================================================#

library(data.table)
library(ggplot2)
library(gridExtra)
library(zoo)
source('functions.R')

#===============================================================================#
# Hardcodes ----
#===============================================================================#

data_dir <- "../output/pcr_returns/"

yrmons <- gsub(
  "[[:space:]]", "",
  as.character(seq(as.yearmon("Jan 1995"), as.yearmon("Dec 2020"), by = 1/12))
)

# check 
unique(fread('../output/impute_ests/bcn_scale_Apr2000.csv')$variable)

#===============================================================================#
# Pull in and format data ----
#===============================================================================#

# Each CSV is stored separately for a given month from different run of 2b_pcr.R
pcr_em <- list()
for (i in yrmons) {
  ym <- as.yearmon(paste0(substr(i, 1, 3), " ", substr(i, 4, 7)))
  tryCatch(
    {pcr_em[[i]] <- fread(paste0(data_dir, "pcr_em_", i, ".csv"))[, yyyymm := ym]},
    error = function(e) ""
  )
}

pcr_mn <- list()
for (i in yrmons) {
  ym <- as.yearmon(paste0(substr(i, 1, 3), " ", substr(i, 4, 7)))
  tryCatch(
    {pcr_mn[[i]] <- fread(paste0(data_dir, "pcr_mn_", i, ".csv"))[, yyyymm := ym]},
    error = function(e) ""
  )
}

pcr_all <- rbind(
    melt(rbindlist(pcr_em), 
        id.vars = c("pc", "n_signals", "yyyymm"), 
        variable.name = "weighting", value.name = "ls_ret")[, type := "EM"],
    melt(rbindlist(pcr_mn), 
        id.vars = c("pc", "n_signals", "yyyymm"), 
        variable.name = "weighting", value.name = "ls_ret")[, type := "Mean"]
)

# Add in the cumulative returns
pcr_all[
  order(type, weighting, pc, yyyymm),
  cumret := log(cumprod(1 + ifelse(is.na(ls_ret/100), 0, ls_ret/100))),
  by = .(type, weighting, pc)
]


#===============================================================================#
# Plots ----
#===============================================================================#

## Aggregate plots ----
scale_gg = 0.55
Npc_max = 80 

# Combine all the data to get average returns and Sharpe ratio over time
agg_data <- pcr_all[
  pc <= Npc_max,
  .( # Get as annualized mean return and std. dev.
    ls_mn = mean(ls_ret, na.rm = T) * 12,
    ls_sd = sd(ls_ret, na.rm = T) * sqrt(12)
  ),
  by = .(type, pc, weighting)
] %>% 
  mutate(
    ls_sharpe = ls_mn / ls_sd
    , weighting = dplyr::case_when(
        weighting == "vw_ls" ~ "Value",
        weighting == "ew_ls" ~ "Equal"
    )
    , type = dplyr::case_when(
      type == 'Mean' ~ 'Naive'
      , type == 'EM' ~ 'EM'
    )
  ) 
   

# All the line plots will have same basic look
plot_base <- ggplot(agg_data, aes(x = pc, colour = weighting, linetype = type)) + 
  theme_bw() + 
  theme(
    legend.position = c(23,85)/100,
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.spacing.y = unit(0.01, 'cm'),
    legend.spacing.x = unit(0.2, 'cm'),
    legend.box = 'horizontal',
    legend.key.size = unit(0.4, 'cm'),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9)
  ) + 
  labs(
    colour = "Stock Weights",
    linetype = "Imputation",
    x = "Number of PCs"
  ) +
  guides(
    colour = guide_legend(order = 1), linetype = guide_legend(order = 2)
  ) +
  scale_color_manual(values = c(MATRED, MATBLUE))

mn <- plot_base + geom_line(aes(y = ls_mn)) + 
  ylab("Annualized Mean Return (%)")
sd <- plot_base + geom_line(aes(y = ls_sd)) +
  ylab("Annualized Std. Dev. (%)") 
sharpe <- plot_base + geom_line(aes(y = ls_sharpe)) + 
  ylab("Annualized Sharpe Ratio") 

out_grid <- marrangeGrob(
  grobs = list(mn, sd, sharpe),
  ncol = 1, nrow = 3,
  top = "LS returns (using deciles) from principal component regressions",
  vp = grid::viewport(width = unit(5.5, "in"), height = unit(10, "in"))
)


ggsave(plot = mn, 
       filename = "../output/plots/pcr_expected_rets.pdf",
       width = 8, height = 5, unit = "in", scale = scale_gg)

ggsave(plot = sharpe, 
       filename = "../output/plots/pcr_sharpes.pdf",
       width = 8, height = 5, unit = "in", scale = scale_gg)



## Cumulative returns over time ----

cumret_plot <- ggplot(
    pcr_all[pc %in% c(2, 10, 25, 37)],
    aes(x = yyyymm, y = cumret)
  ) +
  geom_line(aes(colour = factor(pc), linetype = paste0(weighting, ", ", type))) + 
  theme_bw() + 
  labs(
    x = "Month",
    y = "ln(1 + cumulative mean return)",
    colour = "Number of PCs",
    linetype = "Weighting",
    title = "Hedge portfolio returns over time"
  ) + 
  theme(
    legend.position = c(0.2, 0.8),
    legend.background = element_blank(),
    legend.key = element_blank()
  )





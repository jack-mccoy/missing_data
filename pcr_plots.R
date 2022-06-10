
#===============================================================================
# Packages
#===============================================================================

library(data.table)
library(ggplot2)
library(gridExtra)
library(zoo)

#===============================================================================
# Hardcodes
#===============================================================================

data_dir <- "./output/pcr_returns/"

yrmons <- gsub(
  "[[:space:]]", "",
  as.character(seq(as.yearmon("Jan 1990"), as.yearmon("Dec 2018"), by = 1/12))
)

#===============================================================================
# Pull in and format data
#===============================================================================

# Each CSV is stored separately for a given month from different run of 2b_pcr.R
pcr_em <- list()
for (i in yrmons) {
  ym <- as.yearmon(paste0(substr(i, 1, 3), " ", substr(i, 4, 7)))
  tryCatch(
    {pcr_em[[i]] <- fread(paste0(data_dir, "pcr_", i, ".csv"))[, yyyymm := ym]},
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

# Combine all the data to get average returns and Sharpe ratio over time
agg_data <- pcr_all[
  pc <= 100, 
  .( # Get as annualized mean return and std. dev.
    ls_mn = mean(ls_ret, na.rm = T) * 12,
    ls_sd = sd(ls_ret, na.rm = T) * sqrt(12)
  ),
  by = .(type, pc, weighting)
][, ":="( # Sharpe ratios
  ls_sharpe = ls_mn / ls_sd
)][,
  weighting := dplyr::case_when(
    weighting == "vw_ls" ~ "Value weighted",
    weighting == "ew_ls" ~ "Equal weighted"
  )
]

#===============================================================================
# Plots
#===============================================================================

# Aggregate plots ----

# All the line plots will have same basic look
plot_base <- ggplot(agg_data, aes(x = pc, colour = weighting, linetype = type)) + 
  theme_bw() + 
  theme(
    #legend.position = c(0.8, 0.6),
    legend.background = element_blank(),
    legend.key = element_blank()
  ) + 
  labs(
    colour = "Weighting",
    linetype = "Imputation Type",
    x = "Number of PCs"
  )

mn <- plot_base + geom_line(aes(y = ls_mn)) + 
  ylab("Annualized mean return") +
  ggtitle("Mean returns")
sd <- plot_base + geom_line(aes(y = ls_sd)) +
  ylab("Annualized std. dev.") +
  ggtitle("Standard deviations")
sharpe <- plot_base + geom_line(aes(y = ls_sharpe)) + 
  ylab("Annualized Sharpe ratio") +
  ggtitle("Sharpe ratios from PCR")

# Cumulative returns over time ----

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

#===============================================================================
# Save
#===============================================================================

out_grid <- marrangeGrob(
  grobs = list(mn, sd, sharpe),
  ncol = 1, nrow = 3,
  top = "LS returns (using deciles) from principal component regressions",
  vp = grid::viewport(width = unit(5.5, "in"), height = unit(10, "in"))
)

ggsave(plot = out_grid, 
  filename = "output/plots/pcr_sharpes.pdf",
  width = 6.5, height = 11, unit = "in")

ggsave(plot = sharpe, filename = "output/plots/pcr_sharpes_standalone.pdf",
    width = 7, height = 6, unit = "in")


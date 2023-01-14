
#===============================================================================#
# Packages ----
#===============================================================================#

library(data.table)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(RPostgres)
library(zoo)
source('functions.R')

#===============================================================================#
# Hardcodes ----
#===============================================================================#

out_path <- "../output/"
data_dir <- paste0(out_path, "pcr_returns/")

yrmons <- gsub(
  "[[:space:]]", "",
  as.character(seq(as.yearmon("Jan 1995"), as.yearmon("Dec 2020"), by = 1/12))
)

# check 
unique(fread(paste0(out_path, 'impute_ests/bcn_scale_Apr2000.csv'))$variable)

# create plots folder
dir.create(paste0(out_path, 'plots/'))

#===============================================================================#
# Pull in and format data ----
#===============================================================================#

# Fama-French factors
library(getPass)
user = getPass('wrds username: ')
pass = getPass('wrds password: ')

wrds_con <- dbConnect(Postgres(),
                      host = 'wrds-pgdata.wharton.upenn.edu',
                      port = 9737,
                      dbname = 'wrds',
                      user = user,
                      pass = pass)

ff5_mom <- as.data.table(dbGetQuery(wrds_con, "
    SELECT date,
        mktrf * 100 as mktrf, /* The PCR returns are in pct out of 100 */
        smb * 100 as smb,
        hml * 100 as hml,
        rmw * 100 as rmw,
        cma * 100 as cma,
        umd * 100 as umd
    FROM ff_all.fivefactors_monthly
    ORDER BY date
;"))


# function for importing one (imputation, pcr) spec
import_ret_csv = function(ret_csv_path){
  
  list_csv_files <- list.files(path =ret_csv_path, full.names = TRUE)
  ret = lapply(
    list_csv_files
    , function(fname){
      fname = list_csv_files[1]
      ym = as.yearmon(substr(fname, nchar(fname)-10, nchar(fname)-4))
      ret = fread(fname)
      ret$yyyymm = ym
      ret$date = as.Date(ret$yyyymm)
      ret
    }
  )
  ret = do.call(rbind, ret)
  ret = melt(ret
             , id.vars = c('pc','n_signals','yyyymm','date')
             , variable.name = 'weighting'
             , value.name = 'ls_ret'
             )
}

pcr_em = import_ret_csv(paste0(data_dir, 'em/'))
pcr_ac = import_ret_csv(paste0(data_dir, 'ac/'))
pcr_mn = import_ret_csv(paste0(data_dir, 'mn/'))

pcr_all = rbind(
  pcr_em[, type := "EM Algo"]
  , pcr_ac[, type := "Available Case"]
  , pcr_mn[, type := "Simple Mean"]  
)


# Add in the cumulative returns
pcr_all[
  order(type, weighting, pc, yyyymm),
  cumret := log(cumprod(1 + ifelse(is.na(ls_ret/100), 0, ls_ret/100))),
  by = .(type, weighting, pc)
]


#===============================================================================#
# Regressions for alphas ----
#===============================================================================#

# Merge data to align PC returns with FF5 factors
reg_data <- merge(pcr_all, ff5_mom, by = 'date')

# Set the sequence of PCs to iterate through
pcs <- seq(min(pcr_all$pc), max(pcr_all$pc))

# FF5 + Mom regresions ---

# Alphas from EM-based strategy
alphas_em_ff5 <- melt(
    rbindlist(lapply(pcs, function(x) {
        dat = reg_data[pc == x & type == 'EM Algo']
        alpha_ew <- lm('ls_ret ~ mktrf + smb + hml + rmw + cma + umd',
            dat[weighting == 'ew_ls'])$coefficients['(Intercept)']
        alpha_vw <- lm('ls_ret ~ mktrf + smb + hml + rmw + cma + umd', 
            dat[weighting == 'vw_ls'])$coefficients['(Intercept)']
        return(data.table(pc = x, Equal = alpha_ew, Value = alpha_vw))
    })),
    id = 'pc', value.name = 'alpha_ff5', variable.name = 'weighting'
)[, ":="(type = 'EM Algo')]

# Alphas from mean-based strategy
alphas_mn_ff5 <- melt(
    rbindlist(lapply(pcs, function(x) {
        dat = reg_data[pc == x & type == 'Simple Mean']
        alpha_ew <- lm('ls_ret ~ mktrf + smb + hml + rmw + cma + umd',
            dat[weighting == 'ew_ls'])$coefficients['(Intercept)']
        alpha_vw <- lm('ls_ret ~ mktrf + smb + hml + rmw + cma + umd', 
            dat[weighting == 'vw_ls'])$coefficients['(Intercept)']
        return(data.table(pc = x, Equal = alpha_ew, Value = alpha_vw))
    })), 
    id = 'pc', value.name = 'alpha_ff5', variable.name = 'weighting'
)[, ":="(type = 'Simple Mean')]

# Alphas from available case strategy
alphas_avail_ff5 <- melt(
    rbindlist(lapply(pcs, function(x) {
        dat = reg_data[pc == x & type == 'Available Case']
        alpha_ew <- lm('ls_ret ~ mktrf + smb + hml + rmw + cma + umd',
            dat[weighting == 'ew_ls'])$coefficients['(Intercept)']
        alpha_vw <- lm('ls_ret ~ mktrf + smb + hml + rmw + cma + umd', 
            dat[weighting == 'vw_ls'])$coefficients['(Intercept)']
        return(data.table(pc = x, Equal = alpha_ew, Value = alpha_vw))
    })),
    id = 'pc', value.name = 'alpha_ff5', variable.name = 'weighting'
)[, ":="(type = 'Available Case')]

# CAPM regressions ---

# Alphas from EM-based strategy
alphas_em_capm <- melt(
    rbindlist(lapply(pcs, function(x) {
        dat = reg_data[pc == x & type == 'EM Algo']
        alpha_ew <- lm('ls_ret ~ mktrf',
            dat[weighting == 'ew_ls'])$coefficients['(Intercept)']
        alpha_vw <- lm('ls_ret ~ mktrf', 
            dat[weighting == 'vw_ls'])$coefficients['(Intercept)']
        return(data.table(pc = x, Equal = alpha_ew, Value = alpha_vw))
    })), 
    id = 'pc', value.name = 'alpha_capm', variable.name = 'weighting'
)[, ":="(type = 'EM Algo')]

# Alphas from mean-based strategy
alphas_mn_capm <- melt(
    rbindlist(lapply(pcs, function(x) {
        dat = reg_data[pc == x & type == 'Simple Mean']
        alpha_ew <- lm('ls_ret ~ mktrf',
            dat[weighting == 'ew_ls'])$coefficients['(Intercept)']
        alpha_vw <- lm('ls_ret ~ mktrf', 
            dat[weighting == 'vw_ls'])$coefficients['(Intercept)']
        return(data.table(pc = x, Equal = alpha_ew, Value = alpha_vw))
    })), 
    id = 'pc', value.name = 'alpha_capm', variable.name = 'weighting'
)[, ":="(type = 'Simple Mean')]

# Alphas from available case strategy
alphas_avail_capm <- melt(
    rbindlist(lapply(pcs, function(x) {
        dat = reg_data[pc == x & type == 'Available Case']
        alpha_ew <- lm('ls_ret ~ mktrf',
            dat[weighting == 'ew_ls'])$coefficients['(Intercept)']
        alpha_vw <- lm('ls_ret ~ mktrf', 
            dat[weighting == 'vw_ls'])$coefficients['(Intercept)']
        return(data.table(pc = x, Equal = alpha_ew, Value = alpha_vw))
    })), 
    id = 'pc', value.name = 'alpha_capm', variable.name = 'weighting'
)[, ":="(type = 'Available Case')]

#===============================================================================#
# Plots
#===============================================================================#

# Aggregate plots ----
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
][, ":="(
    ls_sharpe = ls_mn / ls_sd,
    weighting = dplyr::case_when(
        weighting == "vw_ls" ~ "Value",
        weighting == "ew_ls" ~ "Equal"
    )
)][ # Add in alphas
    rbind(alphas_mn_capm, alphas_em_capm, alphas_avail_capm),
    on = c('type', 'pc', 'weighting')
][
    rbind(alphas_mn_ff5, alphas_em_ff5, alphas_avail_ff5),
    on = c('type', 'pc', 'weighting')
][, ":="( # annualize alphas and order type as would be good in legend
    alpha_capm = alpha_capm * 12,
    alpha_ff5 = alpha_ff5 * 12,
    type = factor(type, levels = c('EM Algo', 'Simple Mean', 'Available Case'))
)]

# All the line plots will have same basic look
plot_base <- ggplot(agg_data, aes(x = pc, colour = weighting, linetype = type)) + 
  theme_bw() + 
  theme(
    legend.position = c(27, 85)/100,
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
  scale_linetype_manual(values = c('solid', 'longdash', 'dotted')) +
  scale_size_manual(values = c(0.8, 0.8, 0.5)) +
  scale_color_manual(values = c(MATRED, MATBLUE))

# Specific plots
mn <- plot_base + geom_line(aes(y = ls_mn)) + 
    ylab("Annualized Mean Return (%)")
sd <- plot_base + geom_line(aes(y = ls_sd)) +
    ylab("Annualized Std. Dev. (%)") 
sharpe <- plot_base + geom_line(aes(y = ls_sharpe)) + 
    ylab("Annualized Sharpe Ratio") 
alpha_capm <- plot_base + geom_line(aes(y = alpha_capm)) +
    ylab('Annualized CAPM Alpha (%)')
alpha_ff5 <- plot_base + geom_line(aes(y = alpha_ff5)) +
    ylab('Annualized FF5 + Mom Alpha (%)')

out_grid <- marrangeGrob(
  grobs = list(mn, sd, sharpe),
  ncol = 1, nrow = 3,
  top = "LS returns (using deciles) from principal component regressions",
  vp = grid::viewport(width = unit(5.5, "in"), height = unit(10, "in"))
)

ggsave(plot = mn,
    filename = paste0(out_path, "plots/pcr_expected_rets.pdf"),
    width = 8, height = 5, unit = "in", scale = scale_gg)

ggsave(plot = sharpe, 
    filename = paste0(out_path, "plots/pcr_sharpes.pdf"),
    width = 8, height = 5, unit = "in", scale = scale_gg)

ggsave(plot = alpha_capm,
    filename = paste0(out_path, "plots/pcr_alpha_capm.pdf"),
    width = 8, height = 5, unit = "in", scale = scale_gg)

ggsave(plot = alpha_ff5, 
    filename = paste0(out_path, "plots/pcr_alpha_ff5_mom.pdf"),
    width = 8, height = 5, unit = "in", scale = scale_gg)

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





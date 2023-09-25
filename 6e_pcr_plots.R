
#===============================================================================#
# Packages ----
#===============================================================================#

library(data.table)
library(ggplot2)
library(gridExtra)
library(RPostgres)
library(zoo)
library(dplyr) 
library(stringr)
library(tidyr)

source('functions.R')

#===============================================================================#
# Hardcodes ----
#===============================================================================#

# Project-level file paths
getFilePaths()

pc_ret_path <- paste0(FILEPATHS$out_path, "pca_returns/")
plot_path <- paste0(FILEPATHS$out_path, 'plots/')

dir.create(plot_path, showWarnings = FALSE)

imp_names <- setNames(
    c("EM Algo","Simple Mean","BLLP loc B-XS"),
    c("em","none","bllp6")
)

#===============================================================================#
# Pull in and format data ----
#===============================================================================#

# ff factors
ff5_mom <- fread(paste0(FILEPATHS$data_path, 'raw/ff5_factors.csv'))

# load up specs in pc_ret_path
spec_dat <- data.table(dir = list.files(pc_ret_path))[, ":="(
    forecast = str_split_fixed(basename(dir), '_', 2)[,1],
    imp = str_split_fixed(basename(dir), '_', 2)[,2]
)]
 
# function for importing one spec
import_ret_csv <- function(cur_spec) {
  
    list_csv_files <- list.files(path = paste0(pc_ret_path,cur_spec$dir), 
        full.names = TRUE)
    ret <- lapply(list_csv_files,
        function(fname){
            fname = as.vector(fname)
            ym = str_replace(substr(fname, nchar(fname)-10, nchar(fname)-4), 
                '_', '-')
            ym = as.yearmon(ym)
            ret = fread(fname)
            ret$yyyymm = ym
            ret$date = as.Date(ret$yyyymm)
            ret
        }
    )
    ret <- do.call(rbind, ret)
    ret <- melt(ret,
        id.vars = c('pc','n_signals','yyyymm','date'),
        variable.name = 'weighting',
        value.name = 'ls_ret'
    )
    ret[, weighting := gsub("ew_ls", "Equal", weighting)]
    ret[, weighting := gsub("vw_ls", "Value", weighting)]
    ret$forecast <- cur_spec$forecast
    if (cur_spec$imp %in% names(imp_names)) { 
        ret$imp <- imp_names[cur_spec$imp]
    } else {
        ret$imp <- cur_spec$imp
    }
    return(ret)
}

# import all specs
ret <- rbindlist(lapply(1:dim(spec_dat)[1], 
    function(i) import_ret_csv(spec_dat[i, ]))) 

# cumulative returns
ret[
    order(forecast, imp, weighting, pc, yyyymm),
    cumret := log(cumprod(1 + ifelse(is.na(ls_ret/100), 0, ls_ret/100))),
    by = .(forecast, imp, weighting, pc)
]

#===============================================================================#
# Regressions for alphas ----
#===============================================================================#

# Merge data to align PC returns with FF5 factors
# check timing -ac
reg_data <- merge(ret, ff5_mom, by = 'date') %>% 
    filter(!is.na(ls_ret), !is.na(mktrf))

sum_data <- reg_data[,
    list(
        rbar = summary(lm(ls_ret ~ 1))$coefficients['(Intercept)', 'Estimate']*12,
        alpha_capm = summary(lm(ls_ret ~ mktrf))$coefficients['(Intercept)', 'Estimate']*12,
        alpha_ff5 = summary(
          lm(ls_ret ~ mktrf + smb + hml + rmw + cma + umd)
        )$coefficients['(Intercept)', 'Estimate']*12,
        sharpe = mean(ls_ret)/sd(ls_ret)*sqrt(12),
        vol = sd(ls_ret)*sqrt(12)
    ),
    by = c('forecast', 'imp', 'pc', 'weighting')
][
    imp %in% imp_names
]

sum_data[, imp := factor(imp, levels = imp_names)]

#===============================================================================#
# Plots
#===============================================================================#

# Aggregate plots ----
scale_gg <- 0.55
Npc_max <- 80 

fore_list <- unique(sum_data$forecast)
imp_list <- unique(sum_data$imp)

for (cur_fore in fore_list) {

    if (cur_fore == "pca") {
        pos <- c(25,85)/100
    } else {
        pos <- c(70, 59)/100
    }

    # Elements common to all plots
    common_theme <- theme_bw() + 
        theme(
            legend.position = pos,
            legend.background = element_blank(),
            legend.key = element_blank(),
            legend.spacing.y = unit(0.005, 'cm'),
            legend.spacing.x = unit(0.1, 'cm'),
            legend.box = 'horizontal',
            legend.key.width = unit(0.55, 'cm'),
            legend.key.height = unit(0.38, 'cm'),
            legend.text = element_text(size = 7),
            legend.title = element_text(size = 8)
        )
  
    # All the line plots will have same basic look
    plot_base_main <- ggplot(
            sum_data[forecast == cur_fore & imp != "BLLP loc B-XS"], 
            aes(x = pc, colour = weighting, linetype = imp)
        ) + 
        common_theme  +
        guides(
          colour = guide_legend(order = 1),
          linetype = guide_legend(order = 2)
        ) +
        scale_linetype_manual(values = c('solid', 'longdash', 'dotted')) +
        scale_size_manual(values = c(0.8, 0.8, 0.5)) +
        scale_color_manual(values = c(MATRED, MATBLUE)) + 
        labs(
          colour = "Stock Weights",
          linetype = "Imputation",
          x = "Number of PCs"
        )
    plot_base_appendix <- ggplot(sum_data[forecast == cur_fore], 
            aes(x = pc, colour = weighting, linetype = imp)) + 
        common_theme  +
        guides(
          colour = guide_legend(order = 1),
          linetype = guide_legend(order = 2)
        ) +
        scale_linetype_manual(values = c('solid', 'longdash', 'dotted')) +
        scale_size_manual(values = c(0.8, 0.8, 0.5)) +
        scale_color_manual(values = c(MATRED, MATBLUE)) + 
        labs(
          colour = "Stock Weights",
          linetype = "Imputation",
          x = "Number of PCs"
        )

    # Main figure plots ====
  
    mn_main <- plot_base_main + geom_line(aes(y = rbar)) + 
        scale_y_continuous(breaks = seq(10, 50, 10), limits=c(0,55)) +
        ylab("Annualized Mean Return (%)")
    stdev_main <- plot_base_main + geom_line(aes(y = vol)) +
        ylab("Annualized Std. Dev. (%)") 
    sharpe_main <- plot_base_main + geom_line(aes(y = sharpe)) + 
        scale_y_continuous(breaks = seq(0.5, 3, 0.5), limits=c(0,3.25)) +
        ylab("Annualized Sharpe Ratio") 
    alpha_capm_main <- plot_base_main + geom_line(aes(y = alpha_capm)) +
        ylab('Annualized CAPM Alpha (%)')
    alpha_ff5_main <- plot_base_main + geom_line(aes(y = alpha_ff5)) +
        ylab('Annualized FF5 + Mom Alpha (%)')
  
    ggsave(plot = mn_main,
        filename = paste0(plot_path, cur_fore, "_expected_rets_main.pdf"),
        width = 8, height = 5, unit = "in", scale = scale_gg)
  
    ggsave(plot = sharpe_main, 
        filename = paste0(plot_path, cur_fore, "_sharpes_main.pdf"),
        width = 8, height = 5, unit = "in", scale = scale_gg)
    
    ggsave(plot = alpha_capm_main,
        filename = paste0(plot_path, cur_fore, "_alpha_capm_main.pdf"),
        width = 8, height = 5, unit = "in", scale = scale_gg)
    
    ggsave(plot = alpha_ff5_main, 
        filename = paste0(plot_path, cur_fore, "_alpha_ff5_mom_main.pdf"),
        width = 8, height = 5, unit = "in", scale = scale_gg)

    # Appendix ====
  
    mn_appendix <- plot_base_appendix + geom_line(aes(y = rbar)) + 
        scale_y_continuous(breaks = seq(10, 50, 10), limits=c(0,55)) +
        ylab("Annualized Mean Return (%)")
    stdev_appendix <- plot_base_appendix + geom_line(aes(y = vol)) +
        ylab("Annualized Std. Dev. (%)") 
    sharpe_appendix <- plot_base_appendix + geom_line(aes(y = sharpe)) + 
        scale_y_continuous(breaks = seq(0.5, 3, 0.5), limits=c(0,3.25)) +
        ylab("Annualized Sharpe Ratio") 
    alpha_capm_appendix <- plot_base_appendix + geom_line(aes(y = alpha_capm)) +
        ylab('Annualized CAPM Alpha (%)')
    alpha_ff5_appendix <- plot_base_appendix + geom_line(aes(y = alpha_ff5)) +
        ylab('Annualized FF5 + Mom Alpha (%)')
  
    ggsave(plot = mn_appendix,
        filename = paste0(plot_path, cur_fore, "_expected_rets_appendix.pdf"),
        width = 8, height = 5, unit = "in", scale = scale_gg)
  
    ggsave(plot = sharpe_appendix, 
        filename = paste0(plot_path, cur_fore, "_sharpes_appendix.pdf"),
        width = 8, height = 5, unit = "in", scale = scale_gg)
    
    ggsave(plot = alpha_capm_appendix,
        filename = paste0(plot_path, cur_fore, "_alpha_capm_appendix.pdf"),
        width = 8, height = 5, unit = "in", scale = scale_gg)
    
    ggsave(plot = alpha_ff5_appendix, 
        filename = paste0(plot_path, cur_fore, "_alpha_ff5_mom_appendix.pdf"),
        width = 8, height = 5, unit = "in", scale = scale_gg)
    
}

# Cumulative returns over time ----

cumret_plot <- ggplot(
    ret[pc %in% c(3, 5, max(ret$pc)) & forecast == 'pca'],
    aes(x = yyyymm, y = cumret)
  ) +
  geom_line(aes(colour = factor(pc), linetype = paste0(weighting, ", ", imp))) + 
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

# check
ret %>% 
  group_by(weighting, forecast, imp, pc) %>% 
  summarize(
    rbar = mean(ls_ret)
  ) %>% 
  pivot_wider(names_from = 'forecast', values_from = 'rbar') %>% 
  print(n=100)

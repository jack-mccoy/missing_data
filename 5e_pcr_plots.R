
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

pc_ret_path <- paste0(FILEPATHS$out_path, "pca_forecasts/")
plot_path <- paste0(FILEPATHS$out_path, 'plots/')

dir.create(plot_path, showWarnings = FALSE)

imp_names <- setNames(
    c("EM Algo","Simple Mean","BLLP loc B-XS"),
    c("em","none","bllp6")
)

# The "zoomed-in" PC range at the left
vec_zoomed_in <- 1:25

# Number of PCs to skip for "zoomed out" range
zoom_out_skip <- 10

#===============================================================================#
# Pull in and format data ----
#===============================================================================#

# ff factors
ff5_mom <- fread(paste0(FILEPATHS$data_path, 'raw/ff5_factors.csv'))
setnames(ff5_mom, "date", "yyyymm")
ff5_mom$yyyymm <- as.yearmon(ff5_mom$yyyymm)

# load up specs in pc_ret_path
spec_dat <- data.table(dir = list.files(pc_ret_path))[, ":="(
    forecast = str_split_fixed(basename(dir), '_', 3)[,1],
    imp = str_split_fixed(basename(dir), '_', 3)[,2],
    firmset = str_split_fixed(basename(dir), '_', 3)[,3]
)]
 
# function for importing one spec
import_ret_csv <- function(cur_spec) {
  
    list_csv_files <- list.files(path = paste0(pc_ret_path,cur_spec$dir), 
        full.names = TRUE)
    ret <- rbindlist(lapply(list_csv_files,
        function(fname){
            fname <- as.vector(fname)
            tmp <- fread(fname)
            tmp[, ":="(
                yyyymm = as.yearmon(yyyymm),
                date = as.yearmon(yyyymm)
            )]
            return(tmp)
        }
    ))
    ret$forecast <- cur_spec$forecast
    ret$firmset <- cur_spec$firmset
    if (cur_spec$imp %in% names(imp_names)) { 
        ret$imp <- imp_names[cur_spec$imp]
    } else {
        ret$imp <- cur_spec$imp
    }
    return(ret)
}

# import all firm-month level forecasts
fcasts <- rbindlist(lapply(1:dim(spec_dat)[1], 
    function(i) import_ret_csv(spec_dat[i, ])))

# Filter to make sure that we have all the necessary groups for separate
# then group all the "separate" together
fcasts[ # Count number of groups for separately estimated forecasts
    firmset != "all", 
    n_groups := length(unique(firmset)),
    by = .(yyyymm, forecast, imp, pc)
]
fcasts <- fcasts[ # separately estd forecasts should have 3 groups each
    (n_groups == 3 & firmset != "all") | 
    (is.na(n_groups) & firmset == "all")
][, # Now group all the "separate" estimates together for when forming portfolios
    firmset := dplyr::case_when(
        firmset %in% c("micro", "small", "big") ~ "separate",
        TRUE ~ "all"
    )
]

#===============================================================================
# Create long-short portfolios
#===============================================================================

# Define deciles
fcasts[,
    ":="(
        pctile10 = quantile(Ebh1m, 0.1, na.rm = TRUE),
        pctile90 = quantile(Ebh1m, 0.9, na.rm = TRUE)
    ),
    by = .(yyyymm, forecast, imp, firmset, pc)
][,
    decile := ifelse(Ebh1m <= pctile10, 1, ifelse(Ebh1m >= pctile90, 10, NA)) 
]

# Average portfolios by prediction decile for month, method, and imp
# make negative (short) if lowest decile
ls_ports <- fcasts[
    decile %in% c(1, 10), 
    .(
        ew = mean(ifelse(decile == 1, -1, 1) * bh1m, na.rm = TRUE),
        vw = weighted.mean(ifelse(decile == 1, -1, 1) * bh1m, w = me, na.rm = T)
    ),
    by = .(yyyymm, forecast, imp, firmset, decile, pc)
][,
    .(ew = sum(ew), vw = sum(vw)),
    by = .(yyyymm, forecast, imp, firmset, pc)
]

# Melt to have weighting as a categorical variable (instead of separate cols)
ret <- melt(ls_ports,
    id.vars = c("forecast", "imp", "firmset", "pc", "yyyymm"),
    variable.name = "weighting",
    value.name = "ls_ret")

# cumulative returns
ret[
    order(forecast, imp, weighting, pc, firmset, yyyymm),
    cumret := log(cumprod(1 + ifelse(is.na(ls_ret/100), 0, ls_ret/100))),
    by = .(forecast, imp, weighting, pc, firmset)
]

#===============================================================================#
# Regressions for alphas ----
#===============================================================================#

# Merge data to align PC returns with FF5 factors
# check timing -ac
reg_data <- merge(ret, ff5_mom, by = "yyyymm") %>% 
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
    by = c('forecast', 'imp', 'pc', 'weighting', 'firmset')
][
    imp %in% imp_names
]

# Better labeling and factor ordering for plots
sum_data[, imp := factor(imp, levels = imp_names)]
sum_data[, weighting := dplyr::case_when(
    weighting == "ew" ~ "Equal",
    weighting == "vw" ~ "Value"
)]

#===============================================================================#
# Plots
#===============================================================================#

scale_gg <- 0.55

# Loop through different forecasts (pca and spca) and firm sets (all vs. separate)
fore_list <- unique(sum_data$forecast)
firmset_list <- unique(sum_data$firmset)

# the "zoomed out" PC range
zooms <- list(
    'out' = unique(c(
        1, # one-PC model
        seq(zoom_out_skip, max(sum_data$pc), zoom_out_skip), # Every X PCs
        max(sum_data$pc)) # max PC model
    ),
    'in' = 1:25
)

for (cur_firmset in firmset_list) {
    for (cur_fore in fore_list) {
        for (zoom in names(zooms)) {
    
            if (zoom == "out") {
                pos <- c(75, 62)/100
            } else {
                pos <- c(25, 85)/100
            }
    
            # Elements common to all plots
            common_theme <- theme_bw() + 
                theme(
                    legend.position = pos,
                    legend.background = element_blank(),
                    legend.key = element_blank(),
                    legend.spacing.y = unit(0.005, 'cm'),
                    legend.spacing.x = unit(0.06, 'cm'),
                    legend.box = 'horizontal',
                    legend.key.width = unit(0.55, 'cm'),
                    legend.key.height = unit(0.38, 'cm'),
                    legend.text = element_text(size = 6),
                    legend.title = element_text(size = 7)
                )
      
            # All the line plots will have same basic look
            plot_base_main <- ggplot(
                    sum_data[
                        pc %in% zooms[[zoom]] &
                        forecast == cur_fore &
                        firmset == cur_firmset & 
                        imp != "BLLP loc B-XS"
                    ], 
                    aes(x = pc, colour = weighting, shape = weighting, linetype = imp)
                ) + 
                common_theme  +
                guides(
                    colour = guide_legend(order = 1),
                    shape = guide_legend(order = 1),
                    linetype = guide_legend(order = 2)
                ) +
                scale_linetype_manual(values = c('solid', 'twodash', 'dotted')) +
                scale_size_manual(values = c(0.8, 0.8, 0.5)) +
                scale_color_manual(values = c(MATRED, MATBLUE)) + 
                labs(
                    colour = "Stock Weights",
                    linetype = "Imputation",
                    shape = "Stock Weights",
                    x = "Number of PCs"
                )
            plot_base_appendix <- ggplot(
                    sum_data[
                        pc %in% zooms[[zoom]] &
                        forecast == cur_fore &
                        firmset == cur_firmset
                    ], 
                    aes(x = pc, colour = weighting, shape = weighting, linetype = imp)
                ) + 
                common_theme  +
                guides(
                    colour = guide_legend(order = 1),
                    shape = guide_legend(order = 1),
                    linetype = guide_legend(order = 2)
                ) +
                scale_linetype_manual(values = c('solid', 'twodash', 'dotted')) +
                scale_size_manual(values = c(0.8, 0.8, 0.5)) +
                scale_color_manual(values = c(MATRED, MATBLUE)) + 
                labs(
                    colour = "Stock Weights",
                    linetype = "Imputation",
                    shape = "Stock Weights",
                    x = "Number of PCs"
                )
    
            # Main figure plots ====

            mn_pt <- geom_point(aes(y = rbar))
            sharpe_pt <- geom_point(aes(y = sharpe))
            mn_y <- scale_y_continuous(breaks = seq(0, 50, 10), limits=c(-5,55))
            if (zoom == "out") {
                sharpe_y <- scale_y_continuous(breaks = seq(0, 3, 0.5), limits=c(-0.25,3.25))
            } else {
                mn_pt <- geom_point(aes(y = rbar))
                sharpe_pt <- geom_point(aes(y = sharpe))
                if (cur_fore == "pca") {
                    sharpe_y <- scale_y_continuous(breaks = seq(0, 2, 0.4), limits=c(-0.1,2.1))
                } else {
                    sharpe_y <- scale_y_continuous(breaks = seq(0, 3, 0.5), limits=c(-0.25,3.25))
                }
            }
      
            mn_main <- plot_base_main + 
                geom_line(aes(y = rbar)) + 
                mn_pt + 
                mn_y +
                ylab("Annualized Mean Return (%)")
            sharpe_main <- plot_base_main + geom_line(aes(y = sharpe)) + 
                sharpe_pt + 
                sharpe_y +
                ylab("Annualized Sharpe Ratio") 
            alpha_capm_main <- plot_base_main + geom_line(aes(y = alpha_capm)) +
                geom_point(aes(y = alpha_capm)) + 
                ylab('Annualized CAPM Alpha (%)')
            alpha_ff5_main <- plot_base_main + geom_line(aes(y = alpha_ff5)) +
                geom_point(aes(y = alpha_ff5)) + 
                ylab('Annualized FF5 + Mom Alpha (%)')
      
            ggsave(plot = mn_main,
                filename = paste0(plot_path, 
                    cur_fore, "_", cur_firmset, "_zoom", zoom, 
                    "_expected_rets_main.pdf"),
                width = 8, height = 5, unit = "in", scale = scale_gg)
      
            ggsave(plot = sharpe_main, 
                filename = paste0(plot_path, 
                    cur_fore, "_", cur_firmset, "_zoom", zoom, 
                    "_sharpes_main.pdf"),
                width = 8, height = 5, unit = "in", scale = scale_gg)

            ggsave(plot = alpha_capm_main,
                filename = paste0(plot_path, 
                    cur_fore, "_", cur_firmset, "_zoom", zoom, 
                    "_alpha_capm_main.pdf"),
                width = 8, height = 5, unit = "in", scale = scale_gg)
      
            ggsave(plot = alpha_ff5_main, 
                filename = paste0(plot_path, 
                    cur_fore, "_", cur_firmset, "_zoom", zoom, 
                    "_alpha_ff5_main.pdf"),
                width = 8, height = 5, unit = "in", scale = scale_gg)
    
            # Appendix ====
      
            mn_appendix <- plot_base_appendix + geom_line(aes(y = rbar)) + 
                mn_pt +
                mn_y +
                ylab("Annualized Mean Return (%)")
            sharpe_appendix <- plot_base_appendix + geom_line(aes(y = sharpe)) + 
                sharpe_pt +
                sharpe_y +
                ylab("Annualized Sharpe Ratio") 
            alpha_capm_appendix <- plot_base_appendix + geom_line(aes(y = alpha_capm)) +
                geom_point(aes(y = alpha_capm)) + 
                ylab('Annualized CAPM Alpha (%)')
            alpha_ff5_appendix <- plot_base_appendix + geom_line(aes(y = alpha_ff5)) +
                geom_point(aes(y = alpha_ff5)) + 
                ylab('Annualized FF5 + Mom Alpha (%)')
      
            ggsave(plot = mn_appendix,
                filename = paste0(plot_path, 
                    cur_fore, "_", cur_firmset, "_zoom", zoom,
                    "_expected_rets_appendix.pdf"),
                width = 8, height = 5, unit = "in", scale = scale_gg)
      
            ggsave(plot = sharpe_appendix, 
                filename = paste0(plot_path, 
                    cur_fore, "_", cur_firmset, "_zoom", zoom,
                    "_sharpes_appendix.pdf"),
                width = 8, height = 5, unit = "in", scale = scale_gg)

            ggsave(plot = alpha_capm_appendix,
                filename = paste0(plot_path, 
                    cur_fore, "_", cur_firmset, "_zoom", zoom, 
                    "_alpha_capm_appendix.pdf"),
                width = 8, height = 5, unit = "in", scale = scale_gg)
      
            ggsave(plot = alpha_ff5_appendix, 
                filename = paste0(plot_path, 
                    cur_fore, "_", cur_firmset, "_zoom", zoom, 
                    "_alpha_ff5_appendix.pdf"),
                width = 8, height = 5, unit = "in", scale = scale_gg)
            
        }
    }
}


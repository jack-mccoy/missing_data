# Distribution of correlations, before and after imputation.  Created 2022 06 Andrew

# Environment ----
rm(list = ls())
library(data.table)
library(tidyverse)
library(ggplot2)


datelist = c(1990.5, 2000.5, 2010,5)

## Figures Setup ====
library(gridExtra)
library(latex2exp)
library(extrafont)


MATBLUE = rgb(0,0.4470,0.7410)
MATRED = rgb(0.8500, 0.3250, 0.0980)
MATYELLOW = rgb(0.9290, 0.6940, 0.1250)

NICEBLUE = "#619CFF"
NICEGREEN = "#00BA38"
NICERED = "#F8766D"


chen_theme =   theme_minimal() +
  theme(
    text = element_text(family = "Palatino Linotype")
    , panel.border = element_rect(colour = "black", fill=NA, size=1)    
    
    # Font sizes
    , axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26),
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    legend.text = element_text(size = 18),
    
    # Tweaking legend
    legend.position = c(0.7, 0.8),
    legend.text.align = 0,
    legend.background = element_rect(fill = "white", color = "black"),
    legend.margin = margin(t = 5, r = 20, b = 5, l = 5), 
    legend.key.size = unit(1.5, "cm")
    , legend.title = element_blank()    
  ) 


# read in data ----

## read ----
# round yyyymm to make sure they match (needed because the write precision differs?)

# signals 
signals_mvn = fread('../data/imp_tmp.csv') %>% 
  mutate(yyyymm = round(yyyymm, 3)) 

signals_obs = fread('../data/bc_tmp.csv') %>% 
  mutate(yyyymm = round(yyyymm, 3))

## select dates and make small dataset ----

small = signals_mvn[yyyymm %in% datelist] %>% 
  mutate(imp_type = 'mvn') %>% 
  rbind(
    signals_obs[yyyymm %in% datelist] %>% 
      mutate(imp_type = 'none')
  )
  
keylist = small %>% 
  distinct(
    imp_type, yyyymm
  )
  
# Calculate ----
  
small %>% names()
  
# find correlations  
cor_dat = data.table()
eig_dat = data.table()
for (keyi in 1:dim(keylist)[1]){
  
  # select data
  imp_curr = keylist$imp_type[keyi]
  date_curr = keylist$yyyymm[keyi]
  
  signalmat = small %>% 
    filter(imp_type == imp_curr, yyyymm == date_curr) %>% 
    select(
      -c(permno,yyyymm,imp_type)
    ) %>% 
    as.matrix() 
  
  # calc correlation
  cmat = cor(signalmat, use = 'pairwise.complete')
  clist = cmat[lower.tri(cmat)]
  
  # calc eigendecomp
  eig = eigen(cmat)$values
  pct_var = cumsum(eig)/length(eig)*100
  
  # organize
  dt_curr = data.table(
    imp_type = imp_curr
    , date = date_curr 
    , cor = clist
  )
  
  dt2_curr = data.table(
    imp_type = imp_curr
    , date = date_curr 
    , n_PC = 1:length(eig)
    , pct_var = pct_var
  )
  
  cor_dat = rbind(cor_dat, dt_curr)
  eig_dat = rbind(eig_dat, dt2_curr)
  
} # for keyi



# plot ----

# histogram data
edge = seq(-1,1,0.1)
histdat = cor_dat %>% 
  group_by(imp_type, date) %>% 
  summarize(
    mids = hist(cor, edge)$mids
    , density = hist(cor, edge)$density
  ) %>% 
  mutate(
    imp_type = factor(
      imp_type 
      , levels = c('mvn','none')
      , labels = c('MVN Imputed','Obs Avail Case')
    )
  )

eig_dat2 = eig_dat %>% 
  mutate(
    imp_type = factor(
      imp_type 
      , levels = c('mvn','none')
      , labels = c('MVN Imputed','Obs Avail Case')
    )
  )


for (date_curr in datelist){
  
  # corr dist 
  p = ggplot(histdat %>% filter(date == date_curr), aes(x=mids,y=density)) +
    geom_line(aes(group = imp_type, color = imp_type, linetype = imp_type), size = 2) +
    chen_theme +
    xlab('Pairwise Correlation') +
    ylab('Density') +
    theme(
      legend.position = c(8,8)/10
    ) +
    scale_color_manual(
      values=c(NICEBLUE, 'gray')
    ) +
    scale_linetype_manual(values = c('solid','31'))
  
  p
  
  ggsave(
    filename = paste0('../output/cor_dist_', floor(date_curr), '.pdf')
    , width = 5, height = 4, scale = 1.5, device = cairo_pdf
  )  
  
  # PCA 
  p2 = ggplot(
    eig_dat2 %>% filter(date == date_curr, n_PC <= 10), aes(x=n_PC, y=pct_var)
  ) +
    geom_line(
      aes(group = imp_type, color = imp_type, linetype = imp_type)
      , size = 2
    ) +
    chen_theme +
    xlab('Number of PCs') +
    ylab('Frequency') +
    theme(
      legend.position = c(2.5,8)/10
    ) +
    scale_color_manual(
      values=c(NICEBLUE, 'gray')
    ) +
    scale_linetype_manual(values = c('solid','31'))  +
    labs(x = 'Number of PCs', y = 'Pct Var Explained')  +
    scale_x_continuous(breaks = seq(1, 9 ,2)) +
    coord_cartesian(ylim = c(0,100))
  
  ggsave(
    filename = paste0('../output/pca_', floor(date_curr), '.pdf')
    , width = 5, height = 4, scale = 1.5, device = cairo_pdf
  )   

  
} # for date_curr


# sanity check ----

# check expect share of missing values
small %>% filter(imp_type == 'none') %>% summarize_all(function (x) mean(is.na(x)))
small %>% filter(imp_type == 'mvn') %>% summarize_all(function (x) mean(is.na(x)))

# Distribution of correlations, before and after imputation.  Created 2022 06 Andrew

# Environment ----
rm(list = ls())
library(data.table)
library(tidyverse)
library(ggplot2)
library(zoo)
source('functions.R')

dateliststr = c('Jun1990', 'Jun2000', 'Jun2010')
datelist = as.yearmon(dateliststr)

outpath = '../output/'



## Figures Setup ====
library(gridExtra)
library(latex2exp)
library(extrafont)


# read in data ----

## read in observed ----

# observed signals (transformed)
signals_obs = fread(paste0(outpath, '/bcsignals/bcsignals_none.csv'))
signals_obs = signals_obs[ , yyyymm := as.yearmon(yyyymm)][yyyymm %in% datelist]

# convert to covariance matricies
cor_obs = data.table()
n_obs = data.table()
for (datecurr in datelist){
  
  # select data
  signalmat = signals_obs %>% 
    filter(yyyymm == datecurr) %>% 
    select(-c(permno,yyyymm)) %>% 
    as.matrix() 
  
  # calc correlation
  cmat = cor(signalmat, use = 'pairwise.complete') %>% 
    as.data.table() %>% mutate(yyyymm = as.yearmon(datecurr)) 
  
  # count number of pairwise obs
  obsmat = 1*!is.na(signalmat) %>% as.matrix
  nmat = t(obsmat) %*% obsmat %>% 
    as.data.table() %>% mutate(yyyymm = as.yearmon(datecurr)) 
  
  # bind
  cor_obs = rbind(cor_obs, cmat)
  n_obs = rbind(n_obs, nmat)
  
} # for keyi


## read in imputed ----
cor_imp = data.table()
for (datecurr in dateliststr){
  fname = paste0(outpath, 'impute_ests/estR_', datecurr, '.csv')
  temp = read.csv(fname, row.names = 1) %>% as.matrix() %>% cov2cor() %>% 
    as.data.table() %>% mutate(yyyymm = as.yearmon(datecurr))
  cor_imp = rbind(cor_imp, temp)
} 

## merge ----
cor_dat = rbind(
  cor_obs %>% mutate(imp_type = 'none')
  , cor_imp %>% mutate(imp_type = 'mvn')
) 


# Find eigenvalues and unroll correlations  ---------------------------------------------------------------
  
keylist = cor_dat %>% distinct(imp_type,yyyymm)

clist_dat = data.table()
eig_dat = data.table()
for (keyi in 1:dim(keylist)[1]){
  
  # select data
  imp_curr = keylist$imp_type[keyi]
  date_curr = keylist$yyyymm[keyi]
  cmat = cor_dat %>% 
    filter(imp_type == imp_curr, yyyymm == date_curr) %>% 
    select(-c(imp_type,yyyymm)) %>% 
    as.matrix()
  
  clist = cmat[lower.tri(cmat)]
  
  # calc eigendecomp
  eig = eigen(cmat)$values
  pct_var = cumsum(eig)/length(eig)*100
  
  # organize
  dt_curr = data.table(
    imp_type = imp_curr
    , date = as.yearmon(date_curr)
    , cor = clist
  )
  
  dt2_curr = data.table(
    imp_type = imp_curr
    , date = as.yearmon(date_curr)
    , n_PC = 1:length(eig)
    , pct_var = pct_var
  )
  
  clist_dat = rbind(clist_dat, dt_curr)
  eig_dat = rbind(eig_dat, dt2_curr)

  
} # for keyi


# histogram data
edge = seq(-1,1,0.05)
histdat = clist_dat %>% 
  group_by(imp_type, date) %>% 
  summarize(
    mids = hist(cor, edge)$mids
    , density = hist(cor, edge)$density
  ) %>% 
  mutate(
    imp_type = factor(
      imp_type 
      , levels = c('mvn','none')
      , labels = c('EM Algo','Observed')
    )
  )

eig_dat2 = eig_dat %>% 
  mutate(
    imp_type = factor(
      imp_type 
      , levels = c('mvn','none')
      , labels = c('EM Algo','Observed')
    )
  )

# Plot cor and eig --------------------------------------------------------------------


line_size = 2.8

## corr dist ----
for (date_curr in datelist){
  
  p = ggplot(histdat %>% filter(date == as.yearmon(date_curr)), aes(x=mids,y=density)) +
    geom_line(aes(group = imp_type, color = imp_type, linetype = imp_type), size = line_size) +
    chen_theme +
    xlab('Pairwise Correlation') + ylab('Density') +
    theme(legend.position = c(2.5,8)/10) +
    scale_color_manual(values=c(NICEBLUE, 'gray')) +
    scale_linetype_manual(values = c('solid','31'))+ 
    guides(color=guide_legend(
      keywidth=1.8,
      keyheight=1.2,
      default.unit="cm")
    )
    
  
  ggsave(
    filename = paste0(outpath, '/plots/cor_dist_', floor(date_curr), '.pdf')
    , width = 5, height = 4, scale = 1.5, device = cairo_pdf
  )  
  
} # for datecur


## eig ----
for (date_curr in datelist){
  
  p2 = ggplot(
    eig_dat2 %>% filter(date == as.yearmon(date_curr), n_PC <= 10), aes(x=n_PC, y=pct_var)
  ) +
    geom_line(aes(group = imp_type, color = imp_type, linetype = imp_type), size = line_size) +
    chen_theme +
    xlab('Number of PCs') +
    ylab('Frequency') +
    theme(legend.position = c(2.9,8)/10) +
    scale_color_manual(values=c(NICEBLUE, 'gray')) +
    scale_linetype_manual(values = c('solid','31'))  +
    labs(x = 'Number of PCs', y = 'Pct Var Explained')  +
    scale_x_continuous(breaks = seq(1, 9 ,2)) +
    coord_cartesian(ylim = c(0,100))+ 
    guides(color=guide_legend(
      keywidth=1.8,
      keyheight=1.2,
      default.unit="cm")
    )
  
  ggsave(
    filename = paste0(outpath, 'plots/pca_', floor(date_curr), '.pdf')
    , width = 5, height = 4, scale = 1.5, device = cairo_pdf
  )   
  
  
} # for date_curr



# Plot deltas --------------------------------------------------------------------


## level difference ----
edge = seq(-1,1,0.05/10)

for (date_curr in datelist){

  cimp = cor_dat %>% filter(yyyymm == date_curr, imp_type == 'mvn') %>% select(-c(yyyymm, imp_type)) %>% as.matrix()
  cobs = cor_dat %>% filter(yyyymm == date_curr, imp_type == 'none') %>% select(-c(yyyymm, imp_type)) %>% as.matrix()
  
  dc = cimp-cobs
  dclist = tibble(dc = dc[lower.tri(dc)])
  
  dpct = dc/abs(cobs)*100

  
  plotme = dclist %>% 
    summarize(
      mids = hist(dc, edge)$mids
      , density = hist(dc, edge)$density
    )
  
  p_dc = ggplot(plotme, aes(x=mids, y=density)) +
    geom_line(size = 2, color = 'gray') +
    coord_cartesian(xlim = c(-0.2,+0.2))+
    chen_theme +
    xlab('Difference in Correlation (EM vs Observed)')
  
  ggsave(
    filename = paste0(outpath, '/plots/dcor_dist_', floor(date_curr), '.pdf')
    , width = 7, height = 4, scale = 1.5, device = cairo_pdf
  )  
  
} # end for date_curr


## pct difference ----
edge = seq(-500,500,10)

for (date_curr in datelist){
  
  cimp = cor_dat %>% filter(yyyymm == date_curr, imp_type == 'mvn') %>% select(-c(yyyymm, imp_type)) %>% as.matrix()
  cobs = cor_dat %>% filter(yyyymm == date_curr, imp_type == 'none') %>% select(-c(yyyymm, imp_type)) %>% as.matrix()
  
  dc = cimp-cobs
  dpct = dc/abs(cobs)*100
  
  plotme = tibble(dpct = dpct[lower.tri(dpct)])  %>% 
    filter(dpct >= min(edge), dpct <= max(edge)) %>% 
    summarize(
      mids = hist(dpct, edge)$mids
      , density = hist(dpct, edge)$density
    )
  
  p_dpct = ggplot(plotme, aes(x=mids, y=density)) +
    geom_line(size = 2, color = 'gray') +
    coord_cartesian(xlim = 200*c(-1,+1))+
    chen_theme +
    xlab('% Difference in Corr (EM vs Observed)')
  
  ggsave(
    filename = paste0(outpath, '/plots/dpct_dist_', floor(date_curr), '.pdf')
    , width = 7, height = 4, scale = 1.5, device = cairo_pdf
  )  
  
} # end for date_curr


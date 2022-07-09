# Simulating correlations to illustrate link between dimensionality and imputation slopes
# This does not depend on any of the other scripts
# Created 2022 05


# Setup ====
library(tidyverse)

# random correlations based on Lewandowski et al 
# https://stats.stackexchange.com/questions/124538/how-to-generate-a-large-full-rank-random-correlation-matrix-with-some-strong-cor
random_corr = function(d, betaparam){
  P = matrix(0, nrow = d, ncol = d);           # storing partial correlations
  S = diag(d);
  B = matrix(rbeta(d*d, betaparam, betaparam), nrow = d)
  
  for (k in seq(1,d-1)){
    for (i in seq(k+1,d)){
      P[k,i] = B[k,i]; # sampling from beta
      P[k,i] = (P[k,i]-0.5)*2;     # linearly shifting to [-1, 1]
      p = P[k,i];
      
      if ((k-1) >= 1){
        for (l in seq((k-1),1,-1)){ # converting partial correlation to raw correlation
          p = p * sqrt((1-P[l,i]^2)*(1-P[l,k]^2)) + P[l,i]*P[l,k];
        }
      }
      
      S[k,i] = p;
      S[i,k] = p;
    }
  }
  
  return = S
}

# function for simulating correlations and calculating pca and imputations
simulate_slopes = function(
    d = 75, betaparam = 2, nobs = 38, seed = 111, nsim = 20
){
  # nobs = number of observable signals used in each stock's imputation
  
  set.seed(seed)
  
  # simulate nsim times
  for (simi in 1:nsim){
    
    # make random correlation matrix
    corrmat = random_corr(d,betaparam)  
    corrvec = corrmat[lower.tri(corrmat)] 
    
    # pca decomp
    eig = eigen(corrmat)$values
    pct_var = cumsum(eig)/d*100
    
    # slope distributions
    iall = 1:d
    for (i in 1:d){
      iobs = iall[iall != i]
      iobs = sample(iobs, nobs)
      Cov = corrmat[i,iobs]
      Sinv = solve(corrmat[iobs,iobs])
      
      if (i==1){
        slopes = Cov %*% Sinv
      } else{
        slopes = c(slopes, Cov %*% Sinv)
      }
    }
    
    # store
    if (simi == 1){
      corrall = t(corrvec)
      pct_all = pct_var
      slopes_all = slopes
      
    } else {
      
      corrall = rbind(corrall, t(corrvec))
      pct_all = rbind(pct_all, pct_var)
      slopes_all = rbind(slopes_all, slopes)
      
    } # if simi == 1
    
    
  } # for simi
  
  dat = list(
    corr = corrall, pct_var = pct_all, slopes = slopes_all
  )
  
  return = dat
  
} # function simulate_slopes

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



# Sim + Plots ====

## settings ====
d = 125
p_width = 4
p_height = 1.3
p_scale = 3

## Low Dim Simulation ====

# simulate data
dat = simulate_slopes(betaparam = 1.2, d = d)


# correlation histogram
edge = seq(-1,1,0.1)
temp = hist(dat$corr, edge)
hdat = tibble(
  cmid = temp$mids, freq = temp$counts / sum(temp$counts)
)
p_corr = ggplot(hdat, aes(x=cmid, y=freq)) +
  geom_bar(stat = 'identity', fill = 'gray') +
  labs(x = 'Correlation', y = 'Frequency') +
  chen_theme

# PC decomp
PCmax = 10
pca_dat = tibble(
  num_PCs = 1:d, pct_explained = apply(dat$pct_var, 2, mean) 
)

p_pca = ggplot(
  pca_dat %>% filter(num_PCs <= PCmax), aes(x=num_PCs, y = pct_explained)
) +
  geom_line() +
  chen_theme +
  labs(x = 'Number of PCs', y = 'Pct Var Explained')  +
  scale_x_continuous(breaks = seq(1, 9 ,2)) +
  coord_cartesian(ylim = c(0,100))

# Imputation Slopes
edge = c(-1e4, seq(-2,2,0.1), 1e4)
temp = hist(dat$slopes, edge)
hdat = tibble(
  mids = temp$mids, freq = temp$counts / sum(temp$counts)
) %>% 
  filter(
    abs(mids) <= 10
  )

p_slopes = ggplot(hdat, aes(x=mids, y=freq)) +
  geom_bar(stat = 'identity', fill = 'gray') +
  labs(x = 'Imputation Slope', y = 'Frequency') +
  chen_theme


p_all = grid.arrange(p_corr, p_pca, p_slopes, nrow = 1)


ggsave(
  plot = p_all, filename = '../output/plots/sim impute a.pdf'
  , width = p_width, height = p_height, scale = p_scale, device = cairo_pdf
)



## Med Dim Simulation ====

# simulate data
dat = simulate_slopes(betaparam = 4, d = d)


# correlation histogram
edge = seq(-1,1,0.1)
temp = hist(dat$corr, edge)
hdat = tibble(
  cmid = temp$mids, freq = temp$counts / sum(temp$counts)
)
p_corr = ggplot(hdat, aes(x=cmid, y=freq)) +
  geom_bar(stat = 'identity', fill = 'gray') +
  labs(x = 'Correlation', y = 'Frequency') +
  chen_theme

# PC decomp
PCmax = 10
pca_dat = tibble(
  num_PCs = 1:d, pct_explained = apply(dat$pct_var, 2, mean) 
)

p_pca = ggplot(
  pca_dat %>% filter(num_PCs <= PCmax), aes(x=num_PCs, y = pct_explained)
) +
  geom_line() +
  chen_theme +
  labs(x = 'Number of PCs', y = 'Pct Var Explained')  +
  scale_x_continuous(breaks = seq(1, 9 ,2)) +
  coord_cartesian(ylim = c(0,100))

# Imputation Slopes
edge = c(-1e4, seq(-2,2,0.1), 1e4)
temp = hist(dat$slopes, edge)
hdat = tibble(
  mids = temp$mids, freq = temp$counts / sum(temp$counts)
) %>% 
  filter(
    abs(mids) <= 10
  )

p_slopes = ggplot(hdat, aes(x=mids, y=freq)) +
  geom_bar(stat = 'identity', fill = 'gray') +
  labs(x = 'Imputation Slope', y = 'Frequency') +
  chen_theme


p_all = grid.arrange(p_corr, p_pca, p_slopes, nrow = 1)


ggsave(
  plot = p_all, filename = '../output/plots/sim impute b.pdf'
  , width = p_width, height = p_height, scale = p_scale, device = cairo_pdf
)



## High Dim Simulation ====

# simulate data
dat = simulate_slopes(betaparam = 15, d = d)


# correlation histogram
edge = seq(-1,1,0.1)
temp = hist(dat$corr, edge)
hdat = tibble(
  cmid = temp$mids, freq = temp$counts / sum(temp$counts)
)
p_corr = ggplot(hdat, aes(x=cmid, y=freq)) +
  geom_bar(stat = 'identity', fill = 'gray') +
  labs(x = 'Correlation', y = 'Frequency') +
  chen_theme

# PC decomp
PCmax = 10
pca_dat = tibble(
  num_PCs = 1:d, pct_explained = apply(dat$pct_var, 2, mean) 
)

p_pca = ggplot(
  pca_dat %>% filter(num_PCs <= PCmax), aes(x=num_PCs, y = pct_explained)
) +
  geom_line() +
  chen_theme +
  labs(x = 'Number of PCs', y = 'Pct Var Explained')  +
  scale_x_continuous(breaks = seq(1, 9 ,2)) +
  coord_cartesian(ylim = c(0,100))

# Imputation Slopes
edge = c(-1e4, seq(-2,2,0.1), 1e4)
temp = hist(dat$slopes, edge)
hdat = tibble(
  mids = temp$mids, freq = temp$counts / sum(temp$counts)
) %>% 
  filter(
    abs(mids) <= 10
  )

p_slopes = ggplot(hdat, aes(x=mids, y=freq)) +
  geom_bar(stat = 'identity', fill = 'gray') +
  labs(x = 'Imputation Slope', y = 'Frequency') +
  chen_theme


p_all = grid.arrange(p_corr, p_pca, p_slopes, nrow = 1)

ggsave(
  plot = p_all, filename = '../output/plots/sim impute c.pdf'
  , width = p_width, height = p_height, scale = p_scale, device = cairo_pdf
)



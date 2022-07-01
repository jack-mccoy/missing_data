# Andrew 2022 06: simulate and estimate EM to check algo

# setup ====

library(data.table)
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

# parameter choices ----
d = 5
betaparam = 1.2
n = 1000
pct_miss = 50


## simulate ----
set.seed(116)

# generate true Cov matrix C0
C0 = random_corr(d,betaparam)

# generate x ~ MVN(0, C0)
z = matrix( rnorm(n*d) , nrow = d, ncol = n)
x = t(t(z) %*%chol(C0))

# generate missing, drop stock if all are missing
miss = matrix( runif(n*d) < pct_miss/100 , nrow = d, ncol = n)
drop = apply(miss, 2, prod) 

miss = miss[ , !drop]
x    = x[, !drop]
n = dim(x)[2]

# define obs and missing data matricies
xobs = x
xobs[miss] = NA
xmiss = 0*x
xmiss[!miss] = NA
xhat = xobs
xhat[miss] = xmiss[miss]

## em algo ----

# initialize
Cavail = cov(t(xobs), use = 'pairwise.complete.obs')
Cavail[is.na(Cavail)] = 0
C = Cavail



for (iter in (1:1000)){
  
  # E-step ---
  # set up reduction
  m = matrix(0, nrow = d, ncol = 1)
  s = matrix(0, nrow = d, ncol = d)
  for (i in 1:n){
    
    # estimate missing
    slope = C[miss[,i], !miss[,i] ] %*% solve(C[!miss[,i], !miss[,i]])
    mimiss = slope %*% xobs[!miss[,i], i]
    Vimiss = C[miss[,i], miss[,i]] - slope %*% C[ !miss[,i], miss[,i]]
    # simiss = Vimiss + mimiss %*% t(mimiss)
    
    # unpack into full d x d
    mi = xobs[,i] # plug in obs
    mi[miss[,i]] = mimiss  # fill in missing
    
    vol_correction = matrix(0,d,d)
    vol_correction[miss[,i],miss[,i]] = Vimiss
    si = mi %*% t(mi) + vol_correction

    # "reduce" into sufficient stats
    m = m + mi/n
    s = s + si/n    
      
  } # for i

  # M-step ---
  munew = m
  Cnew = s - m %*% t(m)
  
  
  dist = max(abs(Cnew - C))
  err  = max(abs(Cnew - C0))
  print(c(dist,err))
  
  # update 
  C = Cnew
  mu = munew
  
  if (dist < 1e-8){
    break
  }
  
  
} # for iter



round(C- C0,2)
round(Cavail- C0,2)
round(abs(C-C0)-abs(Cavail- C0),2)



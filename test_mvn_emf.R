# 2023 04 Andrew
# code for testing mvn_emf.f90
# gonna try to make it a bit more flexible


# setup -------------------------------------------------------------------------
source('compile_mvn_emf_win.R')
source('functions.R')
library(dplyr)
set.seed(3)


J = 10 # number of signals
N = 500 # number of stocks
obs_frac = 0.5 # fraction observed

# draw true parameters
Sig0 = matrix(rnorm(J*J), nrow = J)
Sig0 = Sig0 %*% t(Sig0)
Sig0 = round(Sig0*100)

Sig0 = as.matrix(Matrix::nearPD(Sig0)$mat)


mu0 = 1*round(matrix(rnorm(J), nrow = J)*100)

# simulate
L = chol(Sig0)
shock = matrix(rnorm(J*N), nrow = J)
xfull = mu0 %*% matrix(1,1,N) +  t(L) %*% shock
obsID = matrix(runif(J*N), J) < obs_frac
xobs = xfull
xobs[!obsID] = NA

xobs

# check
x = xfull
mu = matrix(rowMeans(x))
xalt = x - mu %*% matrix(1,1,N)
Sig = xalt %*% t(xalt) / N

mu
mu0

Sig
Sig0



# run function ------------------------------------------------------------


y = t(xobs)
# mu = (matrix(0,J,1))
mu = mu0
R0 = diag(1,J,J)

  mu_fixed = 0*mu + 1



  
out = mvn_emf(y, mu, R0, mu_fixed
              , tol = 1e-6, maxiter = as.integer(1000))
out$maxiter





# compare with avail case -------------------------------------------------

SigAC = cov(t(xobs), use = 'pairwise.complete.obs')




err_em = (out$estR - Sig0)
err_ac = (SigAC - Sig0)

as.numeric(mu0)
round(out$estE,2)

# manual
# tempmu = rowMeans(x)
tempmu = mu
meanS = x %*% t(x)/N
Sigman = meanS  - tempmu %*% t(tempmu)


Sig0
out$estR
SigAC
Sigman

qlist = c(0.25, 0.5, 0.75)

data.frame(
  err_em = quantile(abs(err_em),qlist)
 ,  err_ac = quantile(abs(err_ac), qlist)
) %>% 
  mutate(
    em_improvement = err_ac - err_em 
  )


# x[, 1]
# x[,1] %*% t(x[,1])
# 
# x[, 2]
# x[,2] %*% t(x[,2])






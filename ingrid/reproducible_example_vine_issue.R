# Reproducible example for Ingrid ---------------------------------------------
library(VineCopula)
library(texmex)
library(tidyverse)

# function to estimate alpha and beta dependence parameters

# function to calculate observed residuals


# 1. reproduce a particular dataset simulated from a vine copula --------------
#i <- 806
i <- 317
set.seed(i*12)
m1w <- vcw3
m1s <- vcs1
m2w <- vcw2
m2s <- vcs1
# 1a. simulate from both summer and winter vine copulas
xw <- RVineSim(N=nrow(obsresw),RVM=m1w)
xs <- RVineSim(N=nrow(obsress),RVM=m1s)

# 2. refit model for summer and winter separately -----------------------------
refitm1w <- RVineSeqEst(data = xw,RVM=m1w)
refitm1s <- RVineSeqEst(data = xs,RVM=m1s)
# 2a. record the likelihood under null hypothesis
ll_m1w <- refitm1w$pair.logLik
ll_m1s <- refitm1s$pair.logLik
ll_m1 <- ll_m1w+ll_m1s
refitw <- RVineCopSelect(data=xw,Matrix = m2w$Matrix, trunclevel = tree,selectioncrit = "logLik")
refits <- RVineCopSelect(data=xs,Matrix = m2s$Matrix, trunclevel = tree,selectioncrit = "logLik")
# 2b. record likelihood of the more complex model
ll_m2w <- refitw$pair.logLik
ll_m2s <- refits$pair.logLik
f_m2w <- refitw$family
f_m2s <- refits$family
ll_m2 <- ll_m2w+ll_m2s

# 3. verify bivariate link estimation for the three links in the first tree ---
b2s <- BiCopSelect(u1=xs[,1],u2=xs[,2],selectioncrit = "logLik")
b1s <- BiCopSelect(u1=xs[,1],u2=xs[,2],familyset=c(214))
BiCopCompare(u1=xs[,1],u2=xs[,4])
b2s
b1s
b2s$logLik
b1s$logLik
# for 3,4 link, model 1 (summer) gives higher log-likelihood

# using vine copulas should give us the same answer
x <- matrix(c(2,1,0,1),ncol=2)
m2 <- VineCopula::RVineCopSelect(data=xs[,3:4],selectioncrit = "logLik",Matrix=x)
m2a <- VineCopula::RVineCopSelect(data=xs[,3:4],selectioncrit = "logLik",Matrix=x,familyset = c(20))
m1 <- VineCopula::RVineSeqEst(data=xs[,3:4],RVM = m2a)
m2$logLik
m2a$logLik
m1$logLik
# in this case, m2 has a better loglikelihood than m1 as expected

m2 <- VineCopula::RVineCopSelect(data=xs,selectioncrit = "logLik",Matrix=m1s$Matrix)
m1 <- VineCopula::RVineSeqEst(data=xs,RVM = m1s)
m2$logLik
m1$logLik
m2$pair.logLik
m1$pair.logLik
m2
m1


x <- matrix(c(2,1,0,1),ncol=2)
# get a copula for the tree
m <- VineCopula::RVineStructureSelect(data=xs[,c(1,2)],selectioncrit = "logLik")
m2 <- VineCopula::RVineCopSelect(data=xs[,c(1,2)],selectioncrit = "logLik",Matrix=x)
m2a <- VineCopula::RVineCopSelect(data=xs[,c(1,2)],selectioncrit = "logLik",Matrix=x,familyset = c(214))
m1 <- VineCopula::RVineSeqEst(data=xs[,c(1,2)],RVM = m2a)
m2$logLik
m2a$logLik
m1$logLik
m$logLik

# in this case, m2 has a better loglikelihood than m1 as expected
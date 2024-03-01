library(VineCopula) ## Vine copula package (there is another package called rvinecopulib, that is similar)
library(network) ## needed for plotting


d <- 5 #number of variables
n <- 10000 #number of data points


## Specify copula families for all d*(d-1)/2 pairs with Gumbel copulas
## in the first tree, Clayton in the second, Gaussian in the third
## and Student's t in the fourth
family <- c(rep(4,4),rep(3,3),rep(1,2),2)

## Specify first (par) and potentially second parameter (par2) of all
## pair-copulas
par <- c(rep(3,4),rep(2,3),rep(0.5,2),-0.9)
par2 <- c(rep(0,4),rep(0,3),rep(0,2),3)

## Specify the R-vine object corresponding
## to a D-vine with the order (in the first tree): 1-2-3-4-5
RVineMatrix <- D2RVine(1:5,family=family,par=par,par2=par2)


## Simulate n observations from the vine copula
u.sim <- RVineSim(n,RVineMatrix)


## Fit parameters of an R-vine with known structure an d copula families
fit.1 <- RVineSeqEst(u.sim,RVM=RVineMatrix)

## Find structure, pair copula families and fit parameters of an
## R-vine with unknown structure (using the method by Dissmann et al., 2013)
fit.2 <- RVineStructureSelect(u.sim)

## Plot the trees of the R-vine with corresponding copula families:
plot(fit.1,edge.labels = "family")

plot(fit.2,edge.labels = "family")



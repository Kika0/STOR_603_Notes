#library(copula)
library(tidyverse)
library(latex2exp)
library(gridExtra)
library(VineCopula)
library(texmex) # for pollution data

# set theme defaults to be black rectangle
theme_set(theme_bw())
theme_replace(
  panel.spacing = unit(2, "lines"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA) )

# simulate from the copula
cop_refit <- function(i) {
  set.seed(i*12)
logistic <- BiCop(family="104",par=1/2,par2=0)
N <- 1000
u <- copula::rCopula(copula = logistic,n = N)
# record the likelihoods
refit <- fitCopula(copula=logistic,data=u)
ll <- refit@loglik
lltrue <- loglikCopula(param=1/2,u = u,copula = logistic)
return(data.frame(ll=ll,lltrue=lltrue))
}
Nrep <- 100
tmp <- do.call(rbind,lapply(1:Nrep,FUN=cop_refit))

# visualise the results
tmp <- tmp %>% mutate(W=2*(ll-lltrue),Wstar=-2*ll+2+2*lltrue)

pll <- ggplot(tmp) + geom_histogram(aes(x=ll)) + xlab("Log likelihood of fitted parameter")
plltrue <- ggplot(tmp) + geom_histogram(aes(x=lltrue)) + xlab("Log likelihood of true parameter")
pW <- ggplot(tmp) + geom_histogram(aes(x=W)) + xlab("Likelihood ratio test")
pchi <- ggplot(data.frame(rchi=rchisq(1000,df=8))) + geom_histogram(aes(x=rchi)) + xlab(TeX("Random sample from $\\chi^2_1$ distribution"))
pAIC <- ggplot(tmp) + geom_histogram(aes(x=Wstar)) + xlab("AIC difference")

ggsave(pll,filename = "../Documents/p1.png") 
ggsave(plltrue,filename = "../Documents/p2.png") 
ggsave(grid.arrange(pW,pchi,nrow=1),filename = "../Documents/p3.png") 
ggsave(pAIC,filename = "../Documents/p4.png") 

sum(tmp$Wstar>0)

# similar analysis for the vine copula models ----
v <- 0.7
j <- 1
# transform to Laplace margins
winter_lap <- as.data.frame((winter %>% apply(c(2),FUN=row_number))/(nrow(winter)+1)) %>% apply(c(1,2),FUN=unif_laplace_pit) %>% as.data.frame()
summer_lap <- as.data.frame((summer %>% apply(c(2),FUN=row_number))/(nrow(summer)+1)) %>% apply(c(1,2),FUN=unif_laplace_pit) %>% as.data.frame()
pew <-  par_est(df=winter_lap,v=v,given=j,margin = "Normal", method = "sequential2")
obsresw <- (observed_residuals(df = winter_lap,given = j,v = v,a = pew$a,b=pew$b) %>% apply(c(2),FUN=row_number))/(nrow(winter_lap)*(1-v)+1)

pes <-  par_est(df=summer_lap,v=v,given=j,margin = "Normal", method = "sequential2")
obsress <- (observed_residuals(df = summer_lap,given = j,v = v,a = pes$a,b=pes$b) %>% apply(c(2),FUN=row_number))/(nrow(summer_lap)*(1-v)+1)
# run the models from least to most restricted
# 1. separate estimates for winter and summer residuals
vcw1 <- VineCopula::RVineStructureSelect(data=obsresw,selectioncrit = "AIC")
vcs1 <- VineCopula::RVineStructureSelect(data=obsress,selectioncrit = "AIC")

# 2. winter with imposed summer tree structure
vcw2 <- VineCopula::RVineCopSelect(data=obsresw,Matrix=vcs1$Matrix)

# 3. winter with imposed summer tree and family
vcw3 <- VineCopula::RVineSeqEst(data=obsresw,RVM = vcs1)

# 4. joined data
#rvine_join <- VineCopula::RVineStructureSelect(data=rbind(obsresw,obsress)) 
# 4. joined data with imposed summer structure
rvine_join <- VineCopula::RVineSeqEst(data=rbind(obsresw,obsress),RVM=vcw3) 

cop_refit <- function(i,m1w,m1s,m2w,m2s,level="l1") {
  set.seed(i*12)
  # 1. simulate from both vine copulas
  xw <- RVineSim(N=nrow(obsresw),RVM=m1w)
  xs <- RVineSim(N=nrow(obsress),RVM=m1s)
  # 2. refit model for summer and winter separately
  if (level=="l1") {
    refitm1 <- RVineSeqEst(data=rbind(xw,xs),RVM = m1w)
    ll_m1w <- NA
    ll_m1s <- NA
    ll_m1 <- refitm1$logLik
    refitw <- RVineSeqEst(data=xw,RVM = m2w)
    refits <- RVineSeqEst(data=xs,RVM = m2s)
  } else if (level=="l2") {
  refitm1w <- RVineSeqEst(data = xw,RVM=m1w)
  refitm1s <- RVineSeqEst(data = xs,RVM=m1s)
  # 3a. record the likelihood under null hypothesis
  ll_m1w <- refitm1w$logLik
  ll_m1s <- refitm1s$logLik
  ll_m1 <- ll_m1w+ll_m1s
  refitw <- RVineCopSelect(data=xw,Matrix = m2w$Matrix)
  refits <- RVineCopSelect(data=xs,Matrix = m2s$Matrix)
  } else if (level=="l3") {
    # 3a. fit the simpler model M1
    refitm1w <- RVineCopSelect(data = xw,Matrix=m1w$Matrix)
    refitm1s <- RVineCopSelect(data = xs,Matrix=m1s$Matrix)
    # 3a. record the likelihood of the simpler model
    ll_m1w <- refitm1w$logLik
    ll_m1s <- refitm1s$logLik
    ll_m1 <- ll_m1w+ll_m1s
    # 3b. fit the more complex model M2
    refitw <- RVineStructureSelect(data=xw)
    refits <- RVineStructureSelect(data=xs)
  }
  #ll_m1 <- RVineLogLik(data = rbind(xw,xs),RVM=rvine_join)$loglik
  # 3b. record likelihood of the more complex model
   ll_m2w <- refitw$logLik
   ll_m2s <- refits$logLik
   ll_m2 <- ll_m2w+ll_m2s
  return(data.frame(ll_m1w,ll_m1s,ll_m1,ll_m2w,ll_m2s,ll_m2))
}

Nrep <- 50
compare_levels <- function(Nrep=50,l="l3") {
start_time <- Sys.time() 
#model comparison for each level
if (l=="l1") {
tmp <- do.call(rbind,lapply(1:Nrep,FUN=cop_refit,m1w=rvine_join,m1s=rvine_join,m2w=vcw3,m2s=vcs1,level=l))
}
if (l=="l2") {
tmp <- do.call(rbind,lapply(1:Nrep,FUN=cop_refit,m1w=vcw3,m1s=vcs1,m2w=vcw2,m2s=vcs1,level=l))
}
if (l=="l3") {
tmp <- do.call(rbind,lapply(1:Nrep,FUN=cop_refit,m1w=rvine_join,m1s=rvine_join,m2w=vcw1,m2s=vcs1,level=l))
}
end_time <- Sys.time() 
end_time - start_time

# calculate the AIC
npar1 <- 8
npar2 <- 8
tmp <- tmp 
tmp <- tmp %>% mutate(W=2*(ll_m2-ll_m1),AIC1=-2*ll_m1+2*npar1,AIC2=-2*ll_m2+2*npar2)
tmp <- tmp %>% mutate(AICdiff=AIC2-AIC1)

# improve the plots
p1 <- ggplot(tmp %>% dplyr::select(c(ll_m2,ll_m1)) %>% pivot_longer(everything(),names_to = "Model")) + geom_density(aes(x=value,fill=Model),alpha=0.5)  +  scale_fill_manual(name="Model",
                      labels=c(TeX("$M_1$"),
                               TeX("$M_2$")),
                      values=c("black","#C11432")) +
  xlab("Log-likelihood")

pW <- ggplot(tmp) + geom_density(aes(x=W),fill="black",alpha=0.5) + xlab("Likelihood ratio test")

pAIC <- ggplot(tmp) + geom_density(aes(x=AICdiff),fill="black",alpha=0.5) + xlab("AIC difference")

#p1 <- ggplot(tmp %>% dplyr::select(c(ll,lltrue)) %>% pivot_longer(everything(),names_to = "log_likelihood")) + geom_density(aes(x=value,col=log_likelihood))
#pW <- ggplot(tmp) + geom_density(aes(x=W)) + xlab("Likelihood ratio test")
#pAIC <- ggplot(tmp) + geom_histogram(aes(x=AICdiff)) + xlab("AIC difference")

ggsave(p1,filename = paste0("../Documents/p1m",l,".png")) 
ggsave(pW,filename = paste0("../Documents/pW",l,".png")) 
ggsave(pAIC,filename = paste0("../Documents/pAICm",l,".png")) 
}

compare_levels(l="l3",Nrep = 100)

# Analyse for one link: bivariate study ------------------------------------
# simulate a sample from bivariate logistic

BiCopS
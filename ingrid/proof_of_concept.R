library(copula)
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
logistic <- evCopula(family = "tawn",param=c(1/2))
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
pchi <- ggplot(data.frame(rchi=rchisq(1000,df=1))) + geom_histogram(aes(x=rchi)) + xlab(TeX("Random sample from $\\chi^2_1$ distribution"))
pAIC <- ggplot(tmp) + geom_histogram(aes(x=Wstar)) + xlab("AIC difference")

ggsave(pll,filename = "../Documents/p1.png") 
ggsave(plltrue,filename = "../Documents/p2.png") 
ggsave(grid.arrange(pW,pchi,nrow=1),filename = "../Documents/p3.png") 
ggsave(pAIC,filename = "../Documents/p4.png") 

sum(tmp$Wstar>0)

# similar analysis for the vine copula models ----
v <- 0.7
# transform to Laplace margins
winter_lap <- as.data.frame((winter %>% apply(c(2),FUN=row_number))/(nrow(winter)+1)) %>% apply(c(1,2),FUN=unif_laplace_pit) %>% as.data.frame()
summer_lap <- as.data.frame((summer %>% apply(c(2),FUN=row_number))/(nrow(summer)+1)) %>% apply(c(1,2),FUN=unif_laplace_pit) %>% as.data.frame()

# run the models from least to most restricted
# 1. separate estimates for winter and summer residuals
vcw1 <- VineCopula::RVineStructureSelect(data=obsresw,selectioncrit = "AIC")
vcs1 <- VineCopula::RVineStructureSelect(data=obsress,selectioncrit = "AIC")

# 2. winter with imposed summer structure
vcw2 <- VineCopula::RVineCopSelect(data=obsresw,Matrix=vcs1$Matrix)

# 3. winter with imposed summer
vcw3 <- VineCopula::RVineSeqEst(data=obsresw,RVM = vcs1)

# 4. joined data
rvine_join <- VineCopula::RVineStructureSelect(data=rbind(obsresw,obsress)) 

cop_refit <- function(i,m1,m2) {
  set.seed(i*12)
  # simulate from both vine copulas
  x1 <- RVineSim(N=nrow(obsresw),RVM=m1)
  # record the likelihoods
   ll <- RVineLogLik(data=x1,RVM=m1)$loglik
   lltrue <- RVineLogLik(data=x1,RVM=m2)$loglik
  return(data.frame(ll=ll,lltrue=lltrue))
}

Nrep <- 500
# select level of imposed structure for saving the plots
l <- "l2"
start_time <- Sys.time() 
tmp <- do.call(rbind,lapply(1:Nrep,FUN=cop_refit,m1=vcw3,m2=vcw2))
end_time <- Sys.time() 
end_time - start_time

# calculate the AIC
npar1 <- 9
npar2 <- 9
tmp <- tmp %>% mutate(W=2*(ll-lltrue),AIC1=-2*ll+2*npar1,AIC2=-2*lltrue+2*npar2)
tmp <- tmp %>% mutate(AICdiff=AIC1-AIC2)

# improve the plots
p1 <- ggplot(tmp %>% dplyr::select(c(lltrue,ll)) %>% pivot_longer(everything(),names_to = "Model")) + geom_density(aes(x=value,fill=Model),alpha=0.5)  +  scale_fill_manual(name="Model",
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


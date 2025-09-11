# repeat analysis with rvinecopulib package
library(tidyverse)
library(texmex)
library(rvinecopulib)
library(latex2exp)

file.sources = list.files(pattern="*helpers.R")
sapply(file.sources,source,.GlobalEnv)

# set theme defaults to be black rectangle
theme_set(theme_bw())
theme_replace(
  panel.spacing = unit(2, "lines"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA) )

# calculate observed residuals --------------------------------------------
v <- 0.7
j <- 1
# transform to Laplace margins
winter_lap <- as.data.frame((winter %>% apply(c(2),FUN=row_number))/(nrow(winter)+1)) %>% apply(c(1,2),FUN=unif_laplace_pit) %>% as.data.frame()
summer_lap <- as.data.frame((summer %>% apply(c(2),FUN=row_number))/(nrow(summer)+1)) %>% apply(c(1,2),FUN=unif_laplace_pit) %>% as.data.frame()
pew <-  par_est(df=winter_lap,v=v,given=j,margin = "Normal", method = "sequential2")
obsresw <- (observed_residuals(df = winter_lap,given = j,v = v,a = pew$a,b=pew$b) %>% apply(c(2),FUN=row_number))/(nrow(winter_lap)*(1-v)+1)

pes <-  par_est(df=summer_lap,v=v,given=j,margin = "Normal", method = "sequential2")
obsress <- (observed_residuals(df = summer_lap,given = j,v = v,a = pes$a,b=pes$b) %>% apply(c(2),FUN=row_number))/(nrow(summer_lap)*(1-v)+1)

# fit vine copula models for each level ---------------------------------
# 1. separate estimates for winter and summer residuals
vcw1 <- rvinecopulib::vinecop(data=obsresw,selcrit = "aic")
vcs1 <- rvinecopulib::vinecop(data=obsress,selcrit = "aic")

# 2. winter with imposed summer tree structure
vcw2 <- rvinecopulib::vinecop(data=obsresw, structure = vcs1$structure,selcrit = "aic")

# 3. winter with imposed summer tree and family
vcw3 <- rvinecopulib::vinecop(data=obsresw,vinecop_object = vcs1)

# 4. joined data with imposed summer structure and families
rvine_join <- rvinecopulib::vinecop(data=rbind(obsresw,obsress),vinecop_object = vcs1) 


# rewrite function with rvinecopulib functions
cop_refit <- function(i,m1w,m1s,m2w,m2s,level="l1") {
  set.seed(i*13)
  # 1. simulate from both vine copulas
  xw <- rvinecopulib::rvinecop(n = nrow(obsresw),vine = m1w)
  xs <- rvinecopulib::rvinecop(n = nrow(obsress),vine=m1s)
  # 2. refit model for summer and winter separately
  if (level=="l1") {
    refitm1 <- rvinecopulib::vinecop(data=rbind(xw,xs),vinecop_object = m1w,presel = FALSE)
    ll_m1w <- NA
    ll_m1s <- NA
    ll_m1 <- refitm1$loglik
    npar1 <- refitm1$npars
    refitw <- rvinecopulib::vinecop(data=xw,vinecop_object = m2w, presel = FALSE)
    refits <- rvinecopulib::vinecop(data=xs,vinecop_object = m2s, presel = FALSE)
  } else if (level=="l2") {
    refitm1w <- rvinecopulib::vinecop(data = xw,vinecop_object=m1w, presel = FALSE)
    refitm1s <- rvinecopulib::vinecop(data = xs,vinecop_object=m1s, presel = FALSE)
    # 3a. record the likelihood under null hypothesis
    ll_m1w <- refitm1w$loglik
    ll_m1s <- refitm1s$loglik
    ll_m1 <- ll_m1w+ll_m1s
    npar1 <- refitm1w$npars + refitm1s$npars
    refitw <- rvinecopulib::vinecop(data=xw,structure = m2w$structure, selcrit = "loglik", presel = FALSE,family_set = c("onepar","twopar","threepar"))
    refits <- rvinecopulib::vinecop(data=xs,structure = m2s$structure, selcrit = "loglik", presel = FALSE,family_set = c("onepar","twopar","threepar"))
  } else if (level=="l3") {
    # 3a. fit the simpler model M1
    refitm1w <- rvinecopulib::vinecop(data = xw, structure = m1w$structure, selcrit = "loglik", presel = FALSE,family_set = c("onepar","twopar","threepar"))
    refitm1s <- rvinecopulib::vinecop(data = xs, structure = m1s$structure, selcrit = "loglik", presel = FALSE,family_set = c("onepar","twopar","threepar"))
    # 3a. record the likelihood of the simpler model
    ll_m1w <- refitm1w$loglik
    ll_m1s <- refitm1s$loglik
    ll_m1 <- ll_m1w+ll_m1s
    npar1 <- refitm1w$npars + refitm1s$npars
    # 3b. fit the more complex model M2
    refitw <- rvinecopulib::vinecop(data=xw, selcrit = "loglik", presel = FALSE,family_set = c("onepar","twopar","threepar"))
    refits <- rvinecopulib::vinecop(data=xs, selcrit = "loglik", presel = FALSE,family_set = c("onepar","twopar","threepar"))
  }
  # 3b. record likelihood of the more complex model
  ll_m2w <- refitw$loglik
  ll_m2s <- refits$loglik
  ll_m2 <- ll_m2w+ll_m2s
  # check positive log-likleihood difference in level 1
  npar2 <- refitw$npars + refits$npars
  if (level == "l1") {
    ll_m1pair <- summary(refitm1)$loglik
    ll_m2pair <- summary(refits)$loglik + summary(refitw)$loglik
    if (sum(ll_m2pair[1:3]<ll_m1pair[1:3]-10^(-6))>=1) {return(i)}
  }
  return(data.frame(ll_m1w,ll_m1s,ll_m1,ll_m2w,ll_m2s,ll_m2,npar1,npar2))
}

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
  print(end_time - start_time)
  
  # calculate the AIC
  tmp <- tmp 
  tmp <- tmp %>% mutate(W=2*(ll_m2-ll_m1),AIC1=-2*ll_m1+2*npar1,AIC2=-2*ll_m2+2*npar2)
  tmp <- tmp %>% mutate(AICdiff=AIC2-AIC1)
  
  # plot log-likelihood comparison of M1 and M2, likelihood ratio test and AIC difference
  p1 <- ggplot(tmp %>% dplyr::select(c(ll_m2,ll_m1)) %>% pivot_longer(everything(),names_to = "Model")) + geom_density(aes(x=value,fill=Model),alpha=0.5)  +  scale_fill_manual(name="Model",
                                                                                                                                                                                labels=c(TeX("$M_1$"),                                                                                                                                                                                     TeX("$M_2$")),
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
  return(tmp)
}

tmp <- compare_levels(l="l1",Nrep = 100) # setting Nrep = 1000 (time to compute cca 2mins) leads to two such replications when this occurs
# both have negative difference in joe copula tree 2 edge 2


# see what is going on:l1
m1w <- rvine_join # imposed tree, family and parameters
m1s <- rvine_join
m2w <- vcw3 # imposed tree and family
m2s <- vcs1
# explore an example of negative log-likelihood difference
tmp[tmp$W<0,]
i <- 27 # pick the index to reproduce the example
# i <- 288 # shows an iteration with negative log-likelihood difference in tree 1 edge 3
set.seed(i*13)
# 1. simulate from both vine copulas
xw <- rvinecopulib::rvinecop(n = nrow(obsresw),vine = m1w)
xs <- rvinecopulib::rvinecop(n = nrow(obsress),vine=m1s)
# 2. refit model for summer and winter separately

  refitm1 <- rvinecopulib::vinecop(data=rbind(xw,xs),vinecop_object = m1w,presel = FALSE)
  ll_m1w <- NA
  ll_m1s <- NA
  ll_m1 <- refitm1$loglik
  npar1 <- refitm1$npars
  refitw <- rvinecopulib::vinecop(data=xw,vinecop_object = m2w, presel = FALSE)
  refits <- rvinecopulib::vinecop(data=xs,vinecop_object = m2s, presel = FALSE)

# 3b. record likelihood of the more complex model
ll_m2w <- refitw$loglik
ll_m2s <- refits$loglik
ll_m2 <- ll_m2w+ll_m2s

summary(refitm1)
summary(refitw)
summary(refits)
refitm1$loglik
refitw$loglik
refits$loglik
refitm1
refitw
refits
data.frame(llm1=summary(refitm1)$loglik,llm2=summary(refits)$loglik + summary(refitw)$loglik,llm2s= summary(refits)$loglik,llm2w = summary(refitw)$loglik)
# line 5 corresponds to tree 2 edge 2 joe copula


#library(copula)
library(tidyverse)
library(latex2exp)
library(gridExtra)
library(VineCopula)
library(texmex) # for pollution data
library(copula)
library(xtable)

# set theme defaults to be black rectangle
theme_set(theme_bw())
theme_replace(
  panel.spacing = unit(2, "lines"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA) )

# 0. proof of concept on the bivariate copula --------------------------
# simulate from the copula
cop_refit <- function(i) {
  set.seed(i*12)
logistic <- copula::evCopula(family = "huslerReiss",param=1/2)
reprod <- BiCop(family= c(214), par=2.47,par2=0.73)
N <- round(nrow(winter)*0.3)
u <- BiCopSim(N=N,family=c(214),par=2.47,par2=0.73)
# record the likelihoods
m1 <- BiCopSelect(u1=u[,1],u2=u[,2],selectioncrit = "logLik",familyset = c(214))
m1r <- BiCopSelect(u1=u[,1],u2=u[,2],selectioncrit = "logLik",familyset = c(214),rotations=FALSE)
m2 <- BiCopSelect(u1=u[,1],u2=u[,2],selectioncrit = "logLik")
llm1 <- m1$logLik
f1 <- m1$familyname
llm1r <- m1r$logLik
f1r <- m1r$familyname
llm2 <- m2$logLik
f2 <- m2$familyname
return(data.frame(llm1=llm1,family1=f1,llm1r=llm1r,family1rotfalse=f1r,llm2=llm2,family2=f2))
}
Nrep <- 1000
tmp <- do.call(rbind,lapply(1:Nrep,FUN=cop_refit))

# visualise the results
tmp <- tmp %>% mutate(m2m1=(llm2-llm1),m2m1r=llm2-llm1r)
summary(tmp)
tmp %>% filter(m2m1r<0) %>% view()

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

# 1. similar analysis for the vine copula models --------------------------------
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

# plot log-likelihood comparison of M1 and M2, likelihood ratio test and AIC difference
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
# simulate a sample from bivariate Clayton
Clay <- BiCop(family=3,par=1/2)
N <- nrow(obsresw)
x <- BiCopSim(N=N,family=3,par=1/2)
# fit a copula with a fixed family
m2 <- BiCopSelect(u1=x[,1],u2=x[,2])
m1 <- BiCopSelect(u1=x[,1],u2=x[,2],familyset=c(3))
# record the likelihood
m1$logLik
m2$logLik
# explore all possible copulas
BiCopCompare(u1=x[,1],u2=x[,2])

# repeat in a loop
clay_comp <- function(i) {
x <- BiCopSim(N=N,family=4,par=2)
# fit a copula with a fixed family
m2 <- BiCopSelect(u1=x[,1],u2=x[,2])
m1 <- BiCopSelect(u1=x[,1],u2=x[,2],familyset=c(4))
# record the likelihood
return(c(m1$logLik,m2$logLik))
}
Nsim <- 1000
X <- t(sapply(1:Nsim,clay_comp))
X <- as.data.frame(X)
names(X) <- c("ll_m1","ll_m2")
tmp <- X %>% mutate(W=2*(ll_m1-ll_m2))
pW <- ggplot(tmp[tmp$W!=0,]) + geom_density(aes(x=W),fill="black",alpha=0.5) + xlab("Likelihood ratio test")

# as expected
# perform considering different families


cop_refit_l2study <- function(i,m1w,m1s,m2w,m2s,level="l1",tree=1) {
  set.seed(i*12)
  # 1. simulate from both vine copulas
  xw <- RVineSim(N=nrow(obsresw),RVM=m1w)
  xs <- RVineSim(N=nrow(obsress),RVM=m1s)
  # 2. refit model for summer and winter separately
    refitm1w <- RVineSeqEst(data = xw,RVM=m1w)
    refitm1s <- RVineSeqEst(data = xs,RVM=m1s)
    # 3a. record the likelihood under null hypothesis
    ll_m1w <- refitm1w$pair.logLik
    ll_m1s <- refitm1s$pair.logLik
    ll_m1 <- ll_m1w+ll_m1s
    refitw <- RVineCopSelect(data=xw,Matrix = m2w$Matrix, trunclevel = tree,selectioncrit = "logLik")
    refits <- RVineCopSelect(data=xs,Matrix = m2s$Matrix, trunclevel = tree,selectioncrit = "logLik")
  # 3b. record likelihood of the more complex model
  ll_m2w <- refitw$pair.logLik
  ll_m2s <- refits$pair.logLik
  f_m2w <- refitw$family
  f_m2s <- refits$family
  ll_m2 <- ll_m2w+ll_m2s
  return(list(ll_m1w,ll_m1s,ll_m1,ll_m2w,ll_m2s,ll_m2,f_m2w,f_m2s))
}

# testing inside the function
# ll_m2-ll_m1
# 
# m1w <- vcw3
# m1s <- vcs1
# m2w <- vcw2
# m2s <- vcs1

# explore bivariate link
b2s <- BiCopSelect(u1=xs[,3],u2=xs[,4],selectioncrit = "logLik")
b1s <- BiCopSelect(u1=xs[,3],u2=xs[,4],familyset=c(20))
b2w <- BiCopSelect(u1=xw[,3],u2=xw[,4],selectioncrit="logLik")
b1w <- BiCopSelect(u1=xw[,3],u2=xw[,4],familyset=c(20))
b1s$logLik
b1w$logLik
b2s$logLik
b2w$logLik
b1s$logLik+b1w$logLik

# try only with the first tree
# rerun the models from least to most restricted with truncation in T1
# 1. separate estimates for winter and summer residuals
vcw1 <- VineCopula::RVineStructureSelect(data=obsresw,selectioncrit = "AIC",trunclevel = 1)
vcs1 <- VineCopula::RVineStructureSelect(data=obsress,selectioncrit = "AIC",trunclevel = 1)

# 1a. separate estimates with iterative algorithm
vcw1i <- VineCopula::RVineMLE(data=obsresw,RVM = vcw1)
vcs1i <- VineCopula::RVineStructureSelect(data=obsress,selectioncrit = "AIC",trunclevel = 1)


# 2. winter with imposed summer tree structure
vcw2 <- VineCopula::RVineCopSelect(data=obsresw,Matrix=vcs1$Matrix,trunclevel=1)

# 3. winter with imposed summer tree and family
vcw3 <- VineCopula::RVineSeqEst(data=obsresw,RVM = vcs1)

# 4. joined data
#rvine_join <- VineCopula::RVineStructureSelect(data=rbind(obsresw,obsress)) 
# 4. joined data with imposed summer structure
rvine_join <- VineCopula::RVineSeqEst(data=rbind(obsresw,obsress),RVM=vcw3) 

# explore truncated model simulation
ll_m1w
ll_m1s
ll_m2w
ll_m2s
ll_m2-ll_m1

# run Nrep times
Nrep <- 100
X <- lapply(1:Nrep,FUN=cop_refit_l2study,m1w=vcw3,m1s=vcs1,m2w=vcw2,m2s=vcs1,tree=1)
# unlist to get the likelihoods of each link of both models
m1w <- as.data.frame(t(sapply(X=X,FUN=function(x)x[[1]][4,1:3],simplify=TRUE)))
m1s <- as.data.frame(t(sapply(X=X,FUN=function(x)x[[2]][4,1:3],simplify=TRUE)))
m1sw <- as.data.frame(t(sapply(X=X,FUN=function(x)x[[3]][4,1:3],simplify=TRUE)))
m2w <- as.data.frame(t(sapply(X=X,FUN=function(x)x[[4]][4,1:3],simplify=TRUE)))
m2s <- as.data.frame(t(sapply(X=X,FUN=function(x)x[[5]][4,1:3],simplify=TRUE)))
m2sw <- as.data.frame(t(sapply(X=X,FUN=function(x)x[[6]][4,1:3],simplify=TRUE)))
fam2w <- as.data.frame(t(sapply(X=X,FUN=function(x)x[[7]][4,1:3],simplify=TRUE)))
fam2s <- as.data.frame(t(sapply(X=X,FUN=function(x)x[[8]][4,1:3],simplify=TRUE)))
names(m1w) <- names(m1s) <- names(m1sw) <- c("NO2_NO","NO2_PM10","SO2_PM10")
names(m2w) <- names(m2s) <- names(m2sw) <- c("NO2_NO","NO2_PM10","SO2_PM10")
names(fam2w) <- names(fam2s) <- c("NO2_NO","NO2_PM10","SO2_PM10")

# plot log-likelihood difference
tmpf <- rbind(fam2w,fam2s) %>% pivot_longer(cols = c("NO2_NO","NO2_PM10","SO2_PM10"),names_to = "link",values_to = "family_M2") %>% add_row(link=rep(NA,3*Nrep),family_M2=rep(NA,3*Nrep))
tmp1 <- rbind(m1w %>% mutate(season="winter"),m1s %>% mutate(season="summer"),m1sw %>% mutate(season="both")) %>% pivot_longer(cols = c("NO2_NO","NO2_PM10","SO2_PM10"),names_to = "link",values_to = "llM1")
tmp2 <- rbind(m2w,m2s,m2sw) %>% pivot_longer(cols = c("NO2_NO","NO2_PM10","SO2_PM10"),names_to = "link",values_to = "llM2")
tmp <- cbind(tmp1,tmp2 %>% dplyr::select(-link),tmpf %>% dplyr::select(-link)) %>% mutate(W=2*(llM2-llM1))

p1 <- ggplot(tmp) + geom_histogram(aes(x=W),fill="black",alpha=0.5) + xlab("Likelihood ratio test") + facet_wrap(c("season","link"))
# plot only negative values
p2 <- ggplot(tmp[tmp$W<(-10^(-6)),]) + geom_histogram(aes(x=W),fill="black",alpha=0.5) + xlab("Likelihood ratio test") + facet_wrap(c("season","link"))
ggsave(p1,filename="../Documents/ratiotest_linkstudy.png")
ggsave(p2,filename="../Documents/ratiotest_linkstudy_negative.png")

# count zero likelihood difference
xtable(tmp %>% filter(abs(W)<10^(-6)) %>% count(season,link))
# count negative values
xtable(tmp %>% filter((W)<(-10^(-6))) %>% count(season,link))

# print the log-likelihoods of negative values
tmp %>% filter(W<(-10^(-6))) %>% view()

# plot a histogram of families
p1 <- ggplot(tmp) + geom_bar(aes(x=factor(family_M2)),fill="black",alpha=0.5) + xlab("Likelihood ratio test") + facet_wrap(c("season","link"))

# what families are picked for negative values?
p1 <- ggplot(tmp %>% filter(W<(-10^(-6)),season !=c("both"))) + geom_bar(aes(x=factor(family_M2)),fill="black",alpha=0.5) + xlab("Counts of families for negative values of log-likelihood difference") + facet_wrap(c("link","season"),nrow=1)
p1
ggsave(p1,filename="../Documents/negativeis.png",height=3,width=10)
which(tmp$W<(-10^(-6))) 


# explore one index
tmp <- tmp %>% mutate(iteration=rep(rep(1:Nrep,each=3),3))
indeces <- tmp$iteration[(tmp$W<(-10^(-6)))]  
indeces
#i <- indeces[1]
i <- 806
set.seed(i*12)
m1w <- vcw3
m1s <- vcs1
m2w <- vcw2
m2s <- vcs1
# 1. simulate from both vine copulas
xw <- RVineSim(N=nrow(obsresw),RVM=m1w)
xs <- RVineSim(N=nrow(obsress),RVM=m1s)
# 2. refit model for summer and winter separately
refitm1w <- RVineSeqEst(data = xw,RVM=m1w)
refitm1s <- RVineSeqEst(data = xs,RVM=m1s)
# 3a. record the likelihood under null hypothesis
ll_m1w <- refitm1w$pair.logLik
ll_m1s <- refitm1s$pair.logLik
ll_m1 <- ll_m1w+ll_m1s
refitw <- RVineCopSelect(data=xw,Matrix = m2w$Matrix, trunclevel = tree,selectioncrit = "logLik")
refits <- RVineCopSelect(data=xs,Matrix = m2s$Matrix, trunclevel = tree,selectioncrit = "logLik")
# 3b. record likelihood of the more complex model
ll_m2w <- refitw$pair.logLik
ll_m2s <- refits$pair.logLik
f_m2w <- refitw$family
f_m2s <- refits$family
ll_m2 <- ll_m2w+ll_m2s

# explore the output
refitm1w
refitw
ll_m2w
ll_m1w
refitm1s
refits
ll_m2s
ll_m1s
# explore bivariate link of mismatch
b2s <- BiCopSelect(u1=xs[,1],u2=xs[,2],selectioncrit = "logLik")
b1s <- BiCopSelect(u1=xs[,1],u2=xs[,2],familyset=c(214))
b2s$logLik
b1s$logLik
# look at one more pair to verify
b2w <- BiCopSelect(u1=xw[,1],u2=xw[,4],selectioncrit = "logLik")
b1w <- BiCopSelect(u1=xw[,1],u2=xw[,4],familyset=c(214))
b2w$logLik
b1w$logLik

# check where correct
b2s <- BiCopSelect(u1=xs[,3],u2=xs[,4],selectioncrit = "logLik")
b1s <- BiCopSelect(u1=xs[,3],u2=xs[,4],familyset=c(214))
b2s$logLik
b1s$logLik
ll_m2s
ll_m1s

# consider all possible matrices
# start with the winter matrix
vcs1$Matrix

# extract all possible permutations
permutations <- function(n){
  if(n==1){
    return(matrix(1))
  } else {
    sp <- permutations(n-1)
    p <- nrow(sp)
    A <- matrix(nrow=n*p,ncol=n)
    for(i in 1:n){
      A[(i-1)*p+1:p,] <- cbind(i,sp+(sp>=i))
    }
    return(A)
  }
}
all_tree1s <- permutations(4)
# remove symmetric ones
all_tree1s <- all_tree1s[c(1:9,11,13,15),]
# keep reverse for now as a check
giventree_fitfamilies <- function(data=obsresw,tree1=c(1,2,3,4)) {
  # define matrix with default order
  dvine_default <- matrix(data=c(1,0,0,0,4,2,0,0,3,4,3,0,2,3,4,4),nrow=4,byrow = TRUE)
  # order the elements
  dvine_mat <- matrix(case_match(as.numeric(dvine_default),1~tree1[1],2~tree1[2],3~tree1[3],4~tree1[4],0~0),nrow = 4)
  # estimate the family and parameters for dvine_mat tree structure
 y <-  VineCopula::RVineCopSelect(data=data,Matrix=dvine_mat,selectioncrit = "logLik")
 return(list(vc=y,ll=y$logLik))
}

dvine_winter <- apply(all_tree1s,MARGIN = c(1),FUN = giventree_fitfamilies,data=obsresw)
# print log-likelihood of each permutation
loglikw <- sapply(X=dvine_winter,FUN=function(x)x[[2]],simplify=TRUE)

# compare with winter fitted vine copula as a check
vcw1 <- VineCopula::RVineStructureSelect(data=obsresw,selectioncrit = "logLik")
dvine_winter[[2]]

# repeat for summer
dvine_summer <- apply(all_tree1s,MARGIN = c(1),FUN = giventree_fitfamilies,data=obsress)
# print log-likelihood of each permutation
logliks <- sapply(X=dvine_summer,FUN=function(x)x[[2]],simplify=TRUE)

which.max(loglikw+logliks)

dvine_winter[[8]]
dvine_summer[[8]]

# check structure for each conditioning pollutant variable -------------------
# print best fitting tree structure for summer and winter

# compare with joint fit and separate fits
tree_select <- function(j,winterlap=winter_lap,summerlap=summer_lap) {
  # calculate the observed residuals
  v <- 0.7
  pew <-  par_est(df=winterlap,v=v,given=j,margin = "Normal", method = "sequential2")
  obsresw <- (observed_residuals(df = winterlap,given = j,v = v,a = pew$a,b=pew$b) %>% apply(c(2),FUN=row_number))/(nrow(winterlap)*(1-v)+1)
  
  pes <-  par_est(df=summer_lap,v=v,given=j,margin = "Normal", method = "sequential2")
  obsress <- (observed_residuals(df = summerlap,given = j,v = v,a = pes$a,b=pes$b) %>% apply(c(2),FUN=row_number))/(nrow(summerlap)*(1-v)+1)
  # fit overall vine copula
  dvine_winter <- apply(all_tree1s,MARGIN = c(1),FUN = giventree_fitfamilies,data=obsresw)
  # print log-likelihood of each permutation
  loglikw <- sapply(X=dvine_winter,FUN=function(x)x[[2]],simplify=TRUE)
  
  # repeat for summer
  dvine_summer <- apply(all_tree1s,MARGIN = c(1),FUN = giventree_fitfamilies,data=obsress)
  # print log-likelihood of each permutation
  logliks <- sapply(X=dvine_summer,FUN=function(x)x[[2]],simplify=TRUE)
  
 i <-  which.max(loglikw+logliks)
  
 vc_ws <-  dvine_winter[[i]]
 
 # compare with winter fitted vine copulas
vcw <- VineCopula::RVineStructureSelect(data=obsresw,selectioncrit = "logLik")
 vcwd <- dvine_winter[[which.max(loglikw)]]
 vcs <- VineCopula::RVineStructureSelect(data=obsress,selectioncrit = "logLik")
 vc_join <- VineCopula::RVineStructureSelect(data=rbind(obsresw,obsress)) 
  return(list(vcw,vcwd,vcs,vc_ws,vc_join))
}

x <- tree_select(j=2)
x[[1]]
x[[1]]$logLik
x[[2]]

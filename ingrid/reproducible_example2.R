#library(copula)
library(tidyverse)
library(VineCopula)
library(texmex) # for pollution data
library(xtable)
library(rvinecopulib)

# set theme defaults to be black rectangle
theme_set(theme_bw())
theme_replace(
  panel.spacing = unit(2, "lines"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA) )

# proof of concept on the bivariate copula --------------------------
# 1. try with one link
cop_refit <- function(i) {
  set.seed(i*12)
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
Nrep <- 10
tmp <- do.call(rbind,lapply(1:Nrep,FUN=cop_refit))

# repeat with another link
cop_refit <- function(i) {
  set.seed(i*12)
  N <- round(nrow(winter)*0.3)
  u <- BiCopSim(N=N,family=c(5),par=3.08)
  # record the likelihoods
  m1 <- BiCopSelect(u1=u[,1],u2=u[,2],selectioncrit = "logLik",familyset = c(5))
  m1r <- BiCopSelect(u1=u[,1],u2=u[,2],selectioncrit = "logLik",familyset = c(5),rotations=FALSE)
  m2 <- BiCopSelect(u1=u[,1],u2=u[,2],selectioncrit = "logLik")
  llm1 <- m1$logLik
  f1 <- m1$familyname
  llm1r <- m1r$logLik
  f1r <- m1r$familyname
  llm2 <- m2$logLik
  f2 <- m2$familyname
  t <- cor(u,method = "kendall")[1,2]
  return(data.frame("tau"=t,llm1=llm1,family1=f1,llm1r=llm1r,family1rotfalse=f1r,llm2=llm2,family2=f2))
}

tmp1 <- do.call(rbind,lapply(1:Nrep,FUN=cop_refit))

# visualise the results
tmp1 <- tmp1 %>% mutate(m2m1=(llm2-llm1),m2m1r=llm2-llm1r)
summary(tmp1)
library(xtable)
xtable(tmp1 %>% filter(m2m1r<0))
xtable(tmp1 %>% filter(m2m1r<0))
tmp1[tmp$family1==c("indep"),]

# repeat with rvinecopulib package
cop_refit <- function(i) {
  set.seed(i*12)
  N <- round(nrow(winter)*0.3)
  u <- BiCopSim(N=N,family=c(5),par=3.08)
  # record the likelihoods
  m1 <- bicop(data=u,selcrit = "loglik",family_set = "frank")
  m1r <- bicop(data=u,selcrit = "loglik",family_set = "frank",allow_rotations=FALSE)
  m2 <- bicop(data=u,selcrit = "loglik",family_set = c("onepar","twopar"))
  llm1 <- m1$loglik
  f1 <- m1$family
  llm1r <- m1r$loglik
  f1r <- m1r$family
  llm2 <- m2$loglik
  f2 <- m2$family
  t <- cor(u,method = "kendall")[1,2]
  return(data.frame("tau"=t,llm1=llm1,family1=f1,llm1r=llm1r,family1rotfalse=f1r,llm2=llm2,family2=f2))
}


tmp <- do.call(rbind,lapply(1:Nrep,FUN=cop_refit))

# visualise the results
tmp <- tmp %>% mutate(m2m1=(llm2-llm1),m2m1r=llm2-llm1r)
summary(tmp)
library(xtable)
xtable(tmp[tmp1$m2m1r<0,])


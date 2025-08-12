# Reproducible example for Ingrid ---------------------------------------------
library(VineCopula)
library(texmex)
library(tidyverse)
library(LaplacesDemon)

# functions to estimate alpha and beta dependence parameters ------------------
#' Calculate negative log-likelihood of Normal regression
#' 
#' Conditioning on Y_1 being extreme to model Y_2 (conditional model in bivariate case)
#'
#' @param theta A set of 4 parameters: a,b,mu,sig.
#' @param df A dataset with column names of paste0("Y",number).
#' @param given A numeric specifying column name of cond. variable Y1.
#' @param sim A numeric specifying column name of other variable Y2.
#' @return A numeric negative log-likelihood.
#' @export
#'
#' @examples
Y_NLL <- function(theta,df=Y_given1extreme,given=1,sim=2,a_hat=NULL,b_hat=NULL,b_max=1) {
  if (is.null(a_hat)==FALSE) {
    a <- a_hat
  } else {a <- theta[1]}
  if (is.null(b_hat)==FALSE) {
    b <- b_hat
  } else {b <- theta[length(theta)-2]}
  mu <- theta[length(theta)-1]
  sig <- theta[length(theta)]
  Y1 <- df %>% dplyr::select(paste0("Y",given)) %>% pull()
  Y2 <- df %>% dplyr::select(paste0("Y",sim)) %>% pull()
  if (a<(-1) | a>1 | b<0 | b>=b_max | sig<0) {
    log_lik <- (10^6) # low log-likelihood outside Keef bounds
  }
  else {
    log_lik <- sum(log(Y1^b *sig*sqrt(2*pi)) + ((Y2-a*Y1-mu*Y1^b)^2/(2*(Y1^b*sig)^2))  )
  }
  return(log_lik)
}

#' NLL for AGG with different both scale and shape for lower and upper tail
#'
#' @param x A numerical vector of data.
#' @param theta A vector of parameters c(mu,sigl,sigu,deltal,deltau).
#'
#' @return A number of negative log-likelihood.
#' @export
#'
#' @examples NLL_AGG(x=rnorm(50),theta=c(0,1,1,2,2))
NLL_AGG <- function(x,theta) {
  mu <- theta[1]
  sigl <- theta[2]
  sigu <- theta[3]
  deltal <- theta[4]
  deltau <- theta[5]
  z <- c()
  if(sigl<=0 | sigu<=0 | deltal<=0 |deltau<=0 ){return(10e10)}
  C_AGG <-  (sigl/deltal*gamma(1/deltal) + sigu/deltau*gamma(1/deltau)  )^(-1)
  for (i in 1:length(x)) {
    if (x[i]<mu) {
      #   z[i] <- C_AGG*exp(-abs((x[i]-mu)/sigl)^deltal)
      z[i] <- log(C_AGG)-((mu-x[i])/sigl)^deltal 
    }
    else 
      #   z[i] <- C_AGG*exp(-abs((x[i]-mu)/sigu)^deltau)
    {z[i] <- log(C_AGG)-((x[i]-mu)/sigu)^deltau }
  }
  return(-sum(z))
}

# keef constraints for beta
keef_constraint1 <- function(b,a,Y1,Y2) {
  if (b<0 | b>1) {return(10^6)}
  v <- max(Y1)
  ZmAI <- max((Y2-a*Y1)/(Y1^b))
  ZmAD <- max(Y2-Y1)
  return((1-b*ZmAI*v^(b-1)  -a)^2)
}

keef_constraint2 <- function(b,a,Y1,Y2) {
  if (b<0 | b>1) {return(10^6)}
  v <- max(Y1)
  ZmAI <- max((Y2-a*Y1)/(Y1^b))
  ZmAD <- max(Y2-Y1)
  return((1-v^(b-1)*ZmAI+v^(-1)*ZmAD  -a)^2)
}


# generate a table of parameter estimates conditional on (given) each of the specified vector of variables
par_est <- function(df=sims,v=0.99,given=c(1),margin="Normal",method="sequential2", a=NULL, keef_constraints=0) {
  lika <- likb <- likmusig <- a_hat <- b_hat <- bmax <- mu_hat <- sig_hat <- res_var <- c()
  names(df) <- paste0("Y",1:ncol(df))
  d <- ncol(df)
  for (j in given) {
    Y_given1extreme <- df %>% filter(df[,j]>quantile(df[,j],v))
    res <- c(1:d)[-j]
    init_par <- c()
    init_lik <- c()
    for (i in 2:d) {
      # optimise using the initial parameters
      Y1 <- Y_given1extreme[,j]
      Y2 <- Y_given1extreme[,res[i-1]]
      if (method=="sequential2") {
        init_para <- c(0.8,0,1)
        opta <- optim(par=init_para,fn = Y_NLL,df=Y_given1extreme,given=j,sim=res[i-1],b_hat=0,control = list(maxit=2000))
        a_hat <- append(a_hat,opta$par[1])
        lika <- append(lika,opta$value)
        if (1 %in% keef_constraints) {
          b_max1 <- optim(par=0.8,fn = keef_constraint1,a=a_hat[length(a_hat)],Y1=Y1,Y2=Y2,control = list(maxit=2000),lower=0,upper=1, method = "Brent")$par
        } else {b_max1 <- 1}
        if (2 %in% keef_constraints) {
          b_max2 <- optim(par=0.8,fn = keef_constraint2,a=a_hat[length(a_hat)],Y1=Y1,Y2=Y2,control = list(maxit=2000),lower=0, upper=1, method = "Brent")$par
        } else {b_max2 <- 1} 
        b_max <- min(b_max1,b_max2)
        bmax <- append(bmax,b_max)
        init_parb <- c(b_max/2,0,1)
        optb <- optim(par=init_parb,fn = Y_NLL,df=Y_given1extreme,given=j,sim=res[i-1],a_hat=opta$par[1],b_max=b_max,control = list(maxit=2000), method = "Nelder-Mead")
        b_hat <- append(b_hat,optb$par[length(optb$par)-2])
        #     optmusig <- optim(par=init_parb,fn = Y_NLL,df=Y_given1extreme,given=j,sim=res[i-1],a_hat=opta$par[1],b_hat=optb$par[1],control = list(maxit=2000))
        mu_hat <- append(mu_hat,optb$par[length(optb$par)-1])
        sig_hat <- append(sig_hat,optb$par[length(optb$par)])  
        likmusig <- append(likmusig,optb$value)
        likb <- append(likb,optb$value)
      }
      res_var <- append(res_var,res[i-1])
    }
  }
  if (margin=="Normal") {
    if (method %in% c("sequential2")) {
      par_sum <- data.frame("lika" = lika,"likb"=likb,"likmusig"=likmusig,
                            "a" = a_hat, "b" = b_hat,"b_max"=bmax,
                            "mu" = mu_hat,
                            "sig" = sig_hat,
                            "given" = rep(given,each=(d-1)), "res" = res_var)  }
  }
  
  
  return(par_sum)
}

# function to calculate observed residuals
observed_residuals <- function(df=sims,given=1,v=0.99,a=NULL,b=NULL) {
  j <- given
  a_hat <- b_hat <- res_var <- c()
  tmp_z <- tmp_z1 <- c()
  df_orig <- df
  names(df) <- paste0("Y",1:ncol(df))
  d <- ncol(df)
  Y_given1extreme <- df %>% filter(df[,j]>quantile(df[,j],v))
  nv <- nrow(Y_given1extreme)
  res <- c(1:d)[-j]
  init_par <- c()
  for (i in 2:d) {
    if (is.numeric(a) & is.numeric(b)) {
      a_hat <- a[i-1]
      b_hat <- b[i-1]
      res_var <- append(res_var,rep(paste0("Z",res[i-1]),nv))
      Y1 <- Y_given1extreme[,j]
      Y2 <- Y_given1extreme[,res[i-1]]
      tmp_z <- append(tmp_z,(Y2-a_hat*Y1/(Y1^b_hat)))
    }
  }
  Z <- data.frame(res_var,tmp_z) %>% mutate(res_var=factor(res_var,levels=paste0("Z",res))) %>% group_by(res_var) %>% 
    mutate(row = row_number()) %>%
    tidyr::pivot_wider(names_from = res_var, values_from = tmp_z) %>% 
    dplyr::select(-row)
  if (sum(names(df_orig)!=paste0("Y",1:ncol(df_orig)))==ncol(df_orig)) {
    names(Z) <- names(df_orig)[-j]
  }
  return(Z)
}

# 0. fit a vine copula to the observed residuals ------------------------------
v <- 0.7
j <- 1 # change j to repeat analysis cond. on other variables

# transform to Laplace margins
winter_lap <- as.data.frame((winter %>% apply(c(2),FUN=row_number))/(nrow(winter)+1)) %>% apply(c(1,2),FUN=qlaplace) %>% as.data.frame()
summer_lap <- as.data.frame((summer %>% apply(c(2),FUN=row_number))/(nrow(summer)+1)) %>% apply(c(1,2),FUN=qlaplace) %>% as.data.frame()
pew <-  par_est(df=winter_lap,v=v,given=j,margin = "Normal", method = "sequential2")
obsresw <- (observed_residuals(df = winter_lap,given = j,v = v,a = pew$a,b=pew$b) %>% apply(c(2),FUN=row_number))/(nrow(winter_lap)*(1-v)+1)
pes <-  par_est(df=summer_lap,v=v,given=j,margin = "Normal", method = "sequential2")
obsress <- (observed_residuals(df = summer_lap,given = j,v = v,a = pes$a,b=pes$b) %>% apply(c(2),FUN=row_number))/(nrow(summer_lap)*(1-v)+1)
# 01. separate estimates for winter and summer residuals
vcw1 <- VineCopula::RVineStructureSelect(data=obsresw,selectioncrit = "AIC")
vcs1 <- VineCopula::RVineStructureSelect(data=obsress,selectioncrit = "AIC")
# 02. winter with imposed summer tree structure
vcw2 <- VineCopula::RVineCopSelect(data=obsresw,Matrix=vcs1$Matrix)
# 03. winter with imposed summer tree and family
vcw3 <- VineCopula::RVineSeqEst(data=obsresw,RVM = vcs1)
# 04. joined data with imposed summer structure
rvine_join <- VineCopula::RVineSeqEst(data=rbind(obsresw,obsress),RVM=vcw3)

# 1. reproduce a particular dataset simulated from a vine copula --------------
i <- 94
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
tree <- 1 # for simplicity truncate at the first level
refitw <- RVineCopSelect(data=xw,Matrix = m2w$Matrix, trunclevel = tree,selectioncrit = "logLik")
refits <- RVineCopSelect(data=xs,Matrix = m2s$Matrix, trunclevel = tree,selectioncrit = "logLik")
# 2b. record likelihood of the more complex model
ll_m2w <- refitw$pair.logLik
ll_m2s <- refits$pair.logLik
f_m2w <- refitw$family
f_m2s <- refits$family
ll_m2 <- ll_m2w+ll_m2s

# print pair log-likelihoods -------------------------------------------------
ll_m1w
ll_m2w
# link 1,2 for winter shows higher likelihood for Model 1

# 3. verify bivariate link estimation for the link in the first tree ---
b2w <- BiCopSelect(u1=xw[,1],u2=xw[,2],selectioncrit = "logLik")
b1w <- BiCopSelect(u1=xw[,1],u2=xw[,2],familyset=c(214))
BiCopCompare(u1=xw[,1],u2=xw[,2])
b2w
b1w
b2w$logLik
b1w$logLik
# for 1,2 link, model 1 (winter) gives higher log-likelihood (see above vine copula approach)
# however, bivariate links do not identify 214 (rotated Tawn type 2 180) as in BiCopCompare
# instead, a different rotation is selected

# using vine copulas should give us the same answer
x <- matrix(c(2,1,0,1),ncol=2)
m2 <- VineCopula::RVineCopSelect(data=xw[,1:2],selectioncrit = "logLik",Matrix=x)
m2a <- VineCopula::RVineCopSelect(data=xw[,1:2],selectioncrit = "logLik",Matrix=x,familyset = c(214))
m1 <- VineCopula::RVineSeqEst(data=xw[,1:2],RVM = m2a)
m2$logLik
m2a$logLik
m1$logLik
# in this case, m2 has a better log-likelihood than m1 as expected
# by the algorithm design, the answer should be the same as in refitw and refitm1w

m2 <- VineCopula::RVineCopSelect(data=xw,selectioncrit = "logLik",Matrix=m1s$Matrix)
m1 <- VineCopula::RVineSeqEst(data=xw,RVM = m1s)
m2$logLik
m1$logLik
m2$pair.logLik
m1$pair.logLik
m2
m1

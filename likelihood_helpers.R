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
Y_likelihood <- function(theta,df=Y_given_1_extreme,given=1,sim=2,a_hat=NULL,b_hat=NULL,b_max=1) {
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
    log_lik <- (-10^6) # low log-likelihood outside Keef bounds
  }
  else {
    log_lik <- sum(-log(Y1^b *sig*sqrt(2*pi)) + (-(Y2-a*Y1-mu*Y1^b)^2/(2*(Y1^b*sig)^2))  )
  }
  return(log_lik)
}

#' Calculate negative log-likelihood of Normal regression
#' 
#' Conditioning on Y_1 being extreme to model Y_2 (conditional model in bivariate case)
#'
#' @param theta A set of 4 parameters: a,b,mu,sig,delta.
#' @param df A dataset with column names of paste0("Y",number).
#' @param given A numeric specifying column name of cond. variable Y1.
#' @param sim A numeric specifying column name of other variable Y2.
#' @return A numeric negative log-likelihood.
#' @export
#'
#' @examples Y_likelihoodGG(theta=c(0.8,0.2,0.1,2),df=data.frame(Y1=rnorm(5000),Y2=rnorm(5000)) %>% filter(Y1>1))
Y_likelihoodGG <- function(theta,df=Y_given1extreme,given=1,sim=2,a_hat=NULL,b_hat=NULL,b_max=1) {
  if (is.null(a_hat)==FALSE) {
    a <- a_hat
  } else {a <- theta[1]}
  if (is.null(b_hat)==FALSE) {
    b <- b_hat
  } else {b <- theta[length(theta)-3]}
  mu <- theta[length(theta)-2]
  sig <- theta[length(theta)-1]
  delta <- theta[length(theta)]
  Y1 <- df %>% dplyr::select(paste0("Y",given)) %>% pull()
  Y2 <- df %>% dplyr::select(paste0("Y",sim)) %>% pull()
  if (a<(-1) | a>1 | b<0 | b>=b_max | sig<0 | delta<0) {
    log_lik <- (-10^6) # low log-likelihood outside Keef bounds
  }
  else {
    #log_lik <- sum(log(delta)-log(Y1^b *sig*gamma(1/delta)) + (-abs(Y2-a*Y1-mu*Y1^b)^delta/(Y1^b*sig)^delta))  )
    log_lik <- sum(dgnorm(x=(Y2-a*Y1)/(Y1^b),mu=mu,alpha=sig,beta=delta,log = TRUE))
  }
  return(log_lik)
}

#' NLL for Gaussian
#'
#' @param x A numerical vector of data.
#' @param theta A vector of parameters c(mu,sig,delta).
#'
#' @return A number of negative log-likelihood.
#' @export
#'
#' @examples NLL_GenGaus(x=rnorm(50),theta=c(0,1,2))
NLL_norm <- function(x,theta) {
  mu <- theta[1]
  sig <- theta[2]
  if(sig<=0){return(10e10)}
  return(-sum(dnorm(x,mean = mu,sd=sig,log=TRUE)))
}

#' NLL for generalized Gaussian
#'
#' @param x A numerical vector of data.
#' @param theta A vector of parameters c(mu,sig,delta).
#'
#' @return A number of negative log-likelihood.
#' @export
#'
#' @examples NLL_GenGaus(x=rnorm(50),theta=c(0,1,2))
NLL_GenGaus <- function(x,theta) {
  mu <- theta[1]
  sig <- theta[2]
  delta <- theta[3]
  if(sig<=0 | delta<=0){return(10e10)}
  return(-sum(gnorm::dgnorm(x,mu=mu,alpha=sig,beta=delta,log=TRUE)))
}

#' NLL for AGG with different shape for lower and upper tail
#'
#' @param x A numerical vector of data.
#' @param theta A vector of parameters c(mu,sig,deltal,deltau).
#'
#' @return A number of negative log-likelihood.
#' @export
#'
#' @examples NLL_AGGdelta(x=rnorm(50),theta=c(0,1,2,2))
NLL_AGGdelta <- function(x,theta) {
  mu <- theta[1]
  sig <- theta[2]
  deltal <- theta[3]
  deltau <- theta[4]
  z <- c()
  if(sig<=0  | deltal<=0 |deltau<=0 ){return(10e10)}
  C_AGG <-  (sig/deltal*gamma(1/deltal) + sig/deltau*gamma(1/deltau)  )^(-1)
  for (i in 1:length(x)) {
    if (x[i]<mu) {
      #   z[i] <- C_AGG*exp(-abs((x[i]-mu)/sig)^deltal)
      z[i] <- log(C_AGG)-((mu-x[i])/sig)^deltal 
    }
    else 
      #   z[i] <- C_AGG*exp(-abs((x[i]-mu)/sig)^deltau)
      z[i] <- log(C_AGG)-((x[i]-mu)/sig)^deltau 
  }
  return(-sum(z))
}

#' NLL for AGG with different scale for lower and upper tail
#'
#' @param x A numerical vector of data.
#' @param theta A vector of parameters c(mu,sigl,sigu,delta).
#'
#' @return A number of negative log-likelihood.
#' @export
#'
#' @examples NLL_AGGsig(x=rnorm(50),theta=c(0,1,1,2))
NLL_AGGsig <- function(x,theta) {
  mu <- theta[1]
  sigl <- theta[2]
  sigu <- theta[3]
  delta <- theta[4]
  z <- c()
  if(sigl<=0 | sigu<=0 | delta<=0 ){return(10e10)}
  C_AGG <-  ((sigl + sigu)/delta*gamma(1/delta)  )^(-1)
  for (i in 1:length(x)) {
    if (x[i]<mu) {
      #   z[i] <- C_AGG*exp(-abs((x[i]-mu)/sigl)^deltal)
      z[i] <- log(C_AGG)-((mu-x[i])/sigl)^delta 
    }
    else 
      #   z[i] <- C_AGG*exp(-abs((x[i]-mu)/sigu)^deltau)
      z[i] <- log(C_AGG)-((x[i]-mu)/sigu)^delta
  }
  return(-sum(z))
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
NLL_AGG <- function(x,theta,mu_hat=NULL,sigl_hat=NULL,sigu_hat=NULL,deltal_hat=NULL,deltau_hat=NULL) {  
  if (is.null(mu_hat)==FALSE) {
    theta <- append(theta,mu_hat,after=0)
}
  if (is.null(sigl_hat)==FALSE) {
    theta <- append(theta,sigl_hat,after=1)
  } 
  if (is.null(sigu_hat)==FALSE) {
    theta <- append(theta,sigu_hat,after=2)
  }
  if (is.null(deltal_hat)==FALSE) {
    theta <- append(theta,deltal_hat,after=3)
  } 
  if (is.null(deltau_hat)==FALSE) {
    theta <- append(theta,deltau_hat,after=4)
  } 
  mu <- theta[1]
  sigl <- theta[2]
  sigu <- theta[3]
  deltal <- theta[4]
  deltau <- theta[5]
  z <- c()
  if(sigl<=0 | sigu<=0 | deltal<=1 |deltau<=1 ){return(10e10)}
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


# More likelihood functions (mostly not relevant)--------------------------------------------
NLL_exp_norm_noise <- function(d,x,theta) {
  phi <- theta[1]
  sd <- theta[2]
 return(-sum(dnorm(x,mean = exp(-phi*d),sd=sd,log = TRUE)))
}

NLL_expalpha_HT <- function(phi,df = Y_given1extreme, d1j. = d1j,mu1=as.numeric(unlist(mu[,1])),sig1=as.numeric(unlist(sig[,1])),d.=d,given.=given,res.=res,beta1=as.numeric(unlist(b[,1]))) {
  nv <- nrow(df)
  mu1 <- rep(mu1,each=nv)
  sig1 <- rep(sig1,each=nv)
  beta1 <- rep(beta1,each=nv)
  Y1 <- rep(as.numeric(unlist(df[,given.])),d.-1)
  Yj <- as.numeric(unlist(df[,res.]))
  dij. <- rep(d1j.,each=nv)
  # if (a<(-1) | a>1 ) {
  #   log_lik <- (-10^6) # low log-likelihood outside bounds
  # }
  log_lik <- sum(log(sig1*Y1^beta1) + (Yj-exp(-phi*dij.)*Y1-mu1*Y1^beta1)^2/(2*(sig1*Y1^beta1)^2))
  
    # log_lik <- sum(log(sig1) + (Yj-exp(-phi*dij.)*Y1-mu1)^2/(2*sig1^2))
  return(log_lik)
}

NLL_exp_sigmau <- function(phi,x = Z, d1j, mu1=as.numeric(unlist(mu_agg[,1])),sigl1=as.numeric(unlist(sigl[,1])),deltal=2,deltau=2,d.=d) {
  N <- nrow(x)
  mu1 <- rep(mu1, each = N)
  sigl1 <- rep(sigl1, each = N)
  z <- as.numeric(unlist(x))
  dij. <- rep(d1j, each = N)
  log_lik <- rep(NA,length(z))
  sigu <- phi[1]*(1-exp(-(phi[2]*dij.)))
  C_AGG <-  (sigl1/deltal*gamma(1/deltal) + sigu/deltau*gamma(1/deltau)  )^(-1)
  log_lik[x<mu1] <- log(C_AGG[x<mu1])-((mu1[x<mu1]-x[x<mu1])/sigl1[x<mu1])^deltal 
  log_lik[x>=mu1] <- log(C_AGG[x>=mu1])-((x[x>=mu1]-mu1[x>=mu1])/sigu[x>=mu1])^deltau 
  return(-sum(log_lik))
}

NLL_AGG_deltas <- function(theta,x = Z,mu1=as.numeric(unlist(mu_agg[,1])),sigl1=as.numeric(unlist(sigl[,1])),sigu1=as.numeric(unlist(sigu[,1]))) {
  N <- nrow(x)
  deltal <- theta[1]
  deltau <- theta[2]
  if(deltal<=1 |deltau<=1 ){return(10e10)}
  mu1 <- rep(mu1,each=N)
  sigl1 <- rep(sigl1,each=N)
  sigu1 <- rep(sigu1,each=N)
  z <- as.numeric(unlist(x))
  log_lik <- rep(NA,length(z))
  # calculate the constant
  C_AGG <-  (sigl1/deltal*gamma(1/deltal) + sigu1/deltau*gamma(1/deltau)  )^(-1)
 log_lik[z<mu1] <-  log(C_AGG[z<mu1])-((mu1[z<mu1]-z[z<mu1])/sigl1[z<mu1])^deltal 
 log_lik[z>=mu1] <- log(C_AGG[z>=mu1])-((z[z>=mu1]-mu1[z>=mu1])/sigu1[z>=mu1])^deltau
  return(-sum(log_lik))
}

NLL_expalpha_twophi <- function(theta,df = Y_given1extreme, d1j. = d1j,SN.=SN,mu1=as.numeric(unlist(mu[,1])),sig1=as.numeric(unlist(sig[,1])),d.=d,given.=given,res.=res,beta1=as.numeric(unlist(b[,1]))) {
  phi1 <- theta[1] # in the same region as the conditioning site
  phi0 <- theta[2] # in a different region as the conditioning site
  # if(phi1<=0 | phi0<=0 ){return(10e10)}
  nv <- nrow(df)
  mu1 <- rep(mu1,each=nv)
  sig1 <- rep(sig1,each=nv)
  beta1 <- rep(beta1,each=nv)
  Y1 <- rep(as.numeric(unlist(df[,given.])),d.-1)
  Yj <- as.numeric(unlist(df[,res.]))
  dij <- rep(d1j.,each=nv)
  SN. <- rep(SN.,each=nv)
  log_lik <- sum(log(sig1*Y1^beta1) + (Yj-exp(-(phi1*SN.+phi0*!SN.)*dij)*Y1-mu1*Y1^beta1)^2/(2*(Y1^beta1*sig1)^2))
  return(log_lik)
}

NLL_AGGdelta_onestep <- function(x,theta,a_hat=NULL,b_hat=NULL) {
  if (is.null(a_hat)==FALSE) {
    a <- a_hat
  } else {theta <- append(theta,a,after=0)}
  if (is.null(b_hat)==FALSE) {
    b <- b_hat
  } else {theta <- append(theta,b,after=1)}
  mu <- theta[1]
  sig <- theta[2]
  deltal <- theta[3]
  deltau <- theta[4]
  Y1 <- x[,1]
  Y2 <- x[,2]
  obs_res <- (Y2-a*Y1)/(Y1^b)
  if(sig<=0 | deltal<=0 |deltau<=0 | a<(-1) | a>1 | b<0 | b>=1){return(10e10)}
  z <- c()
  C <- (1/deltal*gamma(1/deltal) +1/deltau*gamma(1/deltau) )
  for (i in 1:length(obs_res)) {
    if (obs_res[i]<mu) {
      z[i] <- -log(C)-log(sig)-log(Y1[i]^b)-(abs(obs_res[i]-mu)/sig)^deltal
    }
    else {
      z[i] <- -log(C)-log(sig)-log(Y1[i]^b)-(abs(obs_res[i]-mu)/sig)^deltau
    }
  }
  return(-sum(z))
}

# density function for AGG with different shape for lower and upper tail
dgnormsk <- function(x,mu,sig,deltal,deltau) {
  z <- c()
  C <- (1/deltal*gamma(1/deltal) +1/deltau*gamma(1/deltau) )^(-1)
  for (i in 1:length(x)) {
    if (x[i]<0) {
      z[i] <- C/sig*exp(-abs((x[i]-mu)/sig)^deltal)
    }
    z[i] <- C/sig*exp(-abs((x[i]-mu)/sig)^deltau)
  }
  return(z)
}

# density function for AGG with different scale for lower and upper tail
dgnormsksig <- function(x,theta) {
  mu <- theta[1]
  sigl <- theta[2]
  sigu <- theta[3]
  delta <- theta[4]
  z <- c()
  C_AGG <-  (sigl/delta*gamma(1/delta) + sigu/delta*gamma(1/delta)  )^(-1)
  for (i in 1:length(x)) {
    if (x[i]<mu) {
      z[i] <- C/sigl*exp(-abs((x[i]-mu)/sigl)^delta)
    }
    z[i] <- C/sigu*exp(-abs((x[i]-mu)/sigu)^delta)
  }
  return(z)
}

#' Calculate initial negative log-likelihood of Normal regression
#' 
#' Fixing b=0, calculate the initial values for other 3 parameters
#'
#' @param theta A set of 3 parameters: a,mu,sig
#' @param df A dataset with column names of paste0("Y",number).
#' @param given A numeric specifying column name of cond. variable Y1.
#' @param sim A numeric specifying column name of other variable Y2.#'
#' @return A numeric of negative log-likelihood.
#' @export
#'
#' @examples
Y_likelihood_initial <- function(theta,df=Y_given_1_extreme,given=1,sim=2) {
  a <- theta[1]
  b <- 0
  mu <- theta[2]
  sig <- theta[3]
  Y1 <- df %>% dplyr::select(paste0("Y",given)) %>% pull()
  Y2 <- df %>% dplyr::select(paste0("Y",sim)) %>% pull()
  if (a<(-1) | a>1 | b<0 | b>=1) {
    log_lik <- (-10^6) # low log-likelihood outside bounds
  }
  else {
    log_lik <- sum(-log(Y1^b *sig) + (-(Y2-a*Y1-mu*Y1^b)^2/(2*(Y1^b*sig)^2))  )
  }
  return(log_lik)
}

Y_likelihood_fix_ab <- function(theta,a=1,b=0,df=Y_given_1_extreme,given=1,sim=2) {
  mu <- theta[1]
  sig <- theta[2]
  Y1 <- df %>% dplyr::select(paste0("Y",given)) %>% pull()
  Y2 <- df %>% dplyr::select(paste0("Y",sim)) %>% pull()
  log_lik <- sum(-log(Y1^b *sig) + (-(Y2-a*Y1-mu*Y1^b)^2/(2*(Y1^b*sig)^2))  )
  return(log_lik)
}

#' Calculate negative log-likelihood for Normal regression with Keef constraints
#' 
#' Currently under development of including constraints of Keef et al. (2013) on a,b.
#'
#' @param theta A set of 4 parameters: a,b,mu,sig.
#' @param df A dataset with column names of paste0("Y",1:d).
#' @param given A numeric specifying column name of cond. variable Y1.
#' @param sim A numeric specifying column name of other variable Y2.
#' @param v A numeric quantile threshold.
#'
#' @return A numeric negative log-likelihood.
#' @export
#'
#' @examples
Y_likelihood_constrained <- function(theta,df=Y_given_1_extreme,given=1,sim=2,v=0.99) {
  a <- theta[1]
  b <- theta[2]
  mu <- theta[3]
  sig <- theta[4]
  Y1 <- df %>% dplyr::select(paste0("Y",given)) %>% pull()
  Y2 <- df %>% dplyr::select(paste0("Y",sim)) %>% pull()
  # positive residual quantile
  q <- 0.9
  # zp <- quantile(Y2-Y1,q)
  zp <- max(Y2-Y1)
  # residual quantile
  # z <- quantile( (Y1-a*Y1)/(Y1^b),q)
  z <- max((Y1-a*Y1)/(Y1^b))
  # negative residual quantile
  # zn <- quantile(Y2+Y1,q)
  zn <- max(Y2+Y1)
  
  if (  ( (a<=min(1,1-b*z*v^(b-1),1-v^(b-1)*z+v^(-1)*zp) )|(  ((1-b*z*v^(b-1))<a & a<=1) & (( (1-b^(-1))*(b*z)^(1/(1-b)) *(1-a)^(-b/(1-b)) +zp)>0  ) ) )&
        ( (a<=min(1,1+b*z*v^(b-1),1+v^(b-1)*z-v^(-1)*zn) )|(  ((1+b*z*v^(b-1))<(-a) & (-a)<=1) & (( (1-b^(-1))*(-b*z)^(1/(1-b)) *(1+a)^(-b/(1-b)) -zn)>0  ) ) ) 
  ) {
    
    log_lik <- sum(-log(Y1^b *sig) + (-(Y2-a*Y1-mu*Y1^b)^2/(2*(Y1^b*sig)^2))  )
  }
  else {
    log_lik <- (-10^6)
  }
  return(log_lik)
}

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
Y_likelihood <- function(theta,df=Y_given_1_extreme,given=1,sim=2,a_hat=NULL,b_hat=NULL) {
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
  if (a<(-1) | a>1 | b<0 | b>=1) {
    log_lik <- (-10^6) # low log-likelihood outside bounds
  }
  else {
    log_lik <- sum(-log(Y1^b *sig*sqrt(2*pi)) + (-(Y2-a*Y1-mu*Y1^b)^2/(2*(Y1^b*sig)^2))  )
  }
  return(log_lik)
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
#' @param df A dataset with column names of paste0("Y",number).
#' @param given A numeric specifying column name of cond. variable Y1.
#' @param sim A numeric specifying column name of other variable Y2.
#' @param v A numeric quantile threshold.
#'
#' @return Anumeric negative log-likelihood.
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

DLLL <- function(x,theta) {
  mu <- theta[1]
  sig <- theta[2]
  delta <- theta[3]
  a <- theta[4]
  b <- theta[5]
  Y1 <- x[,1]
  Y2 <- x[,2]
  obs_res <- (Y2-a*Y1)/(Y1^b)
  if(sig<=0 | delta<=0 | a<(-1) | a>1 | b<0 | b>=1){return(10e10)}
  -sum((dgnorm(obs_res,mu=mu,alpha=sig,beta=delta,log=TRUE)-log(Y1^b)))
}

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
    if (x[i]<0) {
      z[i] <- C/sigl*exp(-abs((x[i]-mu)/sigl)^delta)
    }
    z[i] <- C/sigu*exp(-abs((x[i]-mu)/sigu)^delta)
  }
  return(z)
}

# density with different both scale and shape for lower and upper tail
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
    if (x[i]<0) {
   #   z[i] <- C_AGG/sigl*exp(-abs((x[i]-mu)/sigl)^deltal)
   z[i] <- log(C_AGG)-log(sigl)-abs((x[i]-mu)/sigl)^deltal 
    }
 #   z[i] <- C_AGG/sigu*exp(-abs((x[i]-mu)/sigu)^deltau)
    z[i] <- log(C_AGG)-log(sigu)-abs((x[i]-mu)/sigu)^deltau 
  }
  return(-sum(z))
}

DLLLsk <- function(x,theta,a_hat=NULL,b_hat=NULL) {
  if (is.null(a_hat)==FALSE) {
    a <- a_hat
  } else {a <- theta[5]}
  if (is.null(b_hat)==FALSE) {
    b <- b_hat
  } else {b <- theta[6]}
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
    if (obs_res[i]<0) {
      z[i] <- -log(C)-log(sig)-log(Y1[i]^b)-(abs(obs_res[i]-mu)/sig)^deltal
    }
    else {
      z[i] <- -log(C)-log(sig)-log(Y1[i]^b)-(abs(obs_res[i]-mu)/sig)^deltau
    }
  }
  return(-sum(z))
}

#' Calculate negative log-likelihood for AGG distribution
#'
#' @param x A vector of observed values.
#' @param theta A vector of 3 parameters: mu,sig,delta.
#'
#' @return A negative log-likelihood.
#' @export
#'
#' @examples 
#' DLLL2step(x=rgnorm(50),theta=c(0,1,1))
DLLL2step <- function(x,theta) {
  mu <- theta[1]
  sig <- theta[2]
  delta <- theta[3]
  if(sig<=0 | delta<=0 ){return(10e10)}
  return(-sum(dgnorm(x,mu=mu,alpha=sig,beta=delta,log=T)))
}

NLL_exp_norm_noise <- function(d,x,theta) {
  phi <- theta[1]
  sd <- theta[2]
 return(-sum(dnorm(x,mean = exp(-phi*d),sd=sd,log = TRUE)))
}


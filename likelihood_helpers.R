# conditioning on Y_1 being extreme to model Y_2 (conditional model in bivariate case)
# calculates MLE a_hat, b_hat, mu_hat and sig_hat
#' Calculate negative log-likelihood of Normal regression
#'
#' @param theta A set of 4 parameters: a,b,mu,sig
#' @param df A dataset with column names of paste0("Y",number)
#' @param given A numeric specifying column name of cond. variable Y1
#' @param sim A numeric specifying column name of other variable Y2
#'
#' @return
#' @export
#'
#' @examples
Y_likelihood <- function(theta,df=Y_given_1_extreme,given=1,sim=2) {
  a <- theta[1]
  b <- theta[2]
  mu <- theta[3]
  sig <- theta[4]
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

Y_likelihood_initial <- function(theta,df=Y_given_1_extreme,given=1,sim=2) {
  a <- theta[1]
  b <- 0
  mu <- theta[2]
  sig <- theta[3]
  Y1 <- df %>% dplyr::select(paste0("Y",given)) %>% pull()
  Y2 <- df %>% dplyr::select(paste0("Y",sim)) %>% pull()
  if (a<(-1) | a>1 | b<0 | b>=1) {
    log_lik <- (-10^6)
  }
  else {
    log_lik <- sum(-log(Y1^b *sig) + (-(Y2-a*Y1-mu*Y1^b)^2/(2*(Y1^b*sig)^2))  )
  }
  return(log_lik)
}

Y_likelihood_fix_ab <- function(theta,a=1,b=0,df=Y_given_1_extreme,given=1,sim=2) {
  mu <- theta[1]
  sig <- theta[2]
  Y1 <- df %>% dplyr::select(starts_with("Y") & contains(as.character(given))) %>% pull()
  Y2 <- df %>% dplyr::select(starts_with("Y") & contains(as.character(sim))) %>% pull()
  log_lik <- sum(-log(Y1^b *sig) + (-(Y2-a*Y1-mu*Y1^b)^2/(2*(Y1^b*sig)^2))  )
  return(log_lik)
}

Y_likelihood_constrained <- function(theta,df=Y_given_1_extreme,given=1,sim=2) {
  v <- 0.99
  a <- theta[1]
  b <- theta[2]
  mu <- theta[3]
  sig <- theta[4]
  Y1 <- df %>% dplyr::select(starts_with("Y") & contains(as.character(given))) %>% pull()
  Y2 <- df %>% dplyr::select(starts_with("Y") & contains(as.character(sim))) %>% pull()
  #lik <-  prod(1/(Y_1^b *sig)*exp(-(Y_2-a*Y_1-mu*Y_1^b)^2/(2*(Y_1^b*sig)^2)) )
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


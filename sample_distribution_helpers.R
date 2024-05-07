#' Sample from trivariate logistic EVD
#' 
#' Generate random samples from trivariate logistic EVD distribution with X-Y-Z dependence
#' 
#' @param N A number of generated points.
#' @param dep A vector of dependence parameters for each link.
#' @return A dataframe with N rows and 3 columns.
#' @export
#' @examples generate_dep_X_Y_Y_Z(N=10,dep=c(1/2,1/2))
generate_dep_X_Y_Y_Z <- function(N,dep=c(1/2,1/2)) {
  
  dat <- tibble(x1=numeric(),x2=numeric(),x3=numeric(),R=numeric(),u=numeric(),q=numeric(),w1=numeric(),w2=numeric(),w3=numeric())
  # set.seed(12)
  x_y <- evd::rbvevd(N,dep=dep[1],model="log")
  x <- exp(x_y[,1])
  Y <- exp(x_y[,2])
  a <- dep[2]
  # generate z
  to_opt <- function(z) {
    (  (  y^(-(1/a)+1)*(y^(-1/a)+z^(-1/a))^(a-1)*exp(-(y^(-1/a)+z^(-1/a))^a)*exp(1/y)  )-Unif)^2
  }
  z <- c()
  for (i in 1:nrow(x_y)){
    Unif <- runif(1) # generate U
    y <- Y[i]
    z[i] <- optim(par=1,fn=to_opt,lower=0,upper=10^6,method="Brent")$par
  }
  sims <- data.frame(X_1=x,X_2=Y,X_3=z)
  return(sims)
}

#' Sample from trivariate X-logistic-Y-Normal-Z distribution
#'
#' @param N A number of generated points.
#' @param dep A vector of dependence parameters for each link.
#'
#' @return A dataframe with N rows and 3 columns.
#' @export
#'
#' @examples generate_dep_log_norm(N=10,dep=c(1/2,1/2))
generate_dep_log_norm <- function(N,dep=c(1/2,1/2)) {
  dat <- tibble(x1=numeric(),x2=numeric(),x3=numeric(),R=numeric(),u=numeric(),q=numeric(),w1=numeric(),w2=numeric(),w3=numeric())
  # set.seed(12)
  x_y <- evd::rbvevd(N,dep=dep[1],model="log")
  x <- exp(x_y[,1])
  Y <- exp(x_y[,2])
  r <- dep[2]
  # generate z
  U <- runif(N)
  y_N <- qnorm(exp(-1/Y))
  z_N <- qnorm(1-U)*(1-r^2)^(1/2) + r*y_N
  z <- (-log(pnorm(z_N)))^(-1)
  sims <- data.frame(X_1=x,X_2=Y,X_3=z)
  return(sims)
}

#' Sample from trivariate X-Normal-Y-Normal-Z distribution
#'
#' @param N A number of generated points.
#' @param dep A vector of dependence parameters for each link. 
#'
#' @return A dataframe with N rows and 3 columns.
#' @export
#'
#' @examples generate_dep_norm_norm(N,dep=c(1/2,1/2))
generate_dep_norm_norm <- function(N,dep=c(1/2,1/2)) {
  dat <- tibble(x1=numeric(),x2=numeric(),x3=numeric(),R=numeric(),u=numeric(),q=numeric(),w1=numeric(),w2=numeric(),w3=numeric())
  # set.seed(12)
  x_y <- mvrnorm(n=N,mu=c(0,0),Sigma = matrix(c(1,dep[1],dep[1],1),ncol=2))
  x <- (-log(pnorm(x_y[,1])))^(-1)
  Y <- (-log(pnorm(x_y[,2])))^(-1)
  r <- dep[2]
  # generate z
  U <- runif(N)
  y_N <- qnorm(exp(-1/Y))
  z_N <- qnorm(1-U)*(1-r^2)^(1/2) + r*y_N
  z <- (-log(pnorm(z_N)))^(-1)
  sims <- data.frame(X_1=x,X_2=Y,X_3=z)
  return(sims)
}

#' Sample from standard Fréchet
#'
#' @param N A number of samples.
#'
#' @return A vector of length N.
#' @export
#'
#' @examples generate_Y(N=10)
generate_Y <- function(N) {
  U <- runif(N)
  Y <- (-log(U))^(-1)
  sims <- data.frame(Y1=Y)
  return(sims)
}

#' Sample from X|Y logistic distribution
#'
#' Take last column of a dataframe (Y) and add new column X sample from distribution
#' conditional on Y (both have standard Fréchet margins).
#'
#' @param sims A dataframe of variables as columns.
#' @param dep A dependence parameter a of bivariate logistic distribution.
#'
#' @return A dataframe with ncol(sims)+1 columns.
#' @export
#'
#' @examples link_log(sims=generate_Y(N=10),dep=1/2)
link_log <- function(sims,dep=1/2) {
  Y <- sims %>% pull(-1)
  N <- length(Y)
  a <- dep
  # generate x
  to_opt <- function(x) {
    (  (  y^(-(1/a)+1)*(y^(-1/a)+x^(-1/a))^(a-1)*exp(-(y^(-1/a)+x^(-1/a))^a)*exp(1/y)  )-Unif)^2
  }
  x <- c()
  for (i in 1:N){
    Unif <- runif(1) # generate U
    y <- Y[i]
    x[i] <- optim(par=1,fn=to_opt,lower=0,upper=10^6,method="Brent")$par
  }
  sims <- sims %>% mutate(X=x)
  names(sims) <- paste0("Y",1:ncol(sims)) # rename columns
  return(sims)
}

#' Sample from X|Y normal distribution
#'
#' Take last column of a dataframe (Y) and add new column X sample from distribution
#' conditional on Y (both have standard Fréchet margins).
#'
#' @param sims A dataframe of variables as columns.
#' @param dep A dependence parameter a of bivariate normal distribution.
#'
#' @return A dataframe with ncol(sims)+1 columns.
#' @export
#'
#' @examples link_norm(sims=generate_Y(N=10),dep=1/2)
link_norm <- function(sims,dep=1/2) {
  Y <- sims %>% pull(-1) 
  N <- length(Y)
  r <- dep
  # generate x
  U <- runif(N)
  y_N <- qnorm(exp(-1/Y)) # PIT from standard Fréchet to uniform
  x_N <- qnorm(1-U)*(1-r^2)^(1/2) + r*y_N
  x <- (-log(pnorm(x_N)))^(-1)
  sims <- sims %>% mutate(X=x)
  names(sims) <- paste0("Y",1:ncol(sims)) # rename columns
  return(sims)
}

#' Transform from Fréchet to Laplace margins
#'
#' @param x A number sampled from Fréchet distribution.
#'
#' @return x transformed to Laplace margin.
#' @export
#'
#' @examples frechet_laplace_pit(x=generate_Y(N=1))
frechet_laplace_pit <- function(x) {
  if (exp(-1/x)<0.5) { y <- log(2*exp(-1/x)) }
  else { y <- -log(2*(1-exp(-1/x))) }
  return(y)
}

#' Transform back from Laplace to Fréchet margins
#'
#' @param x A number sampled from Laplace distribution.
#'
#' @return x transformed to Fréchet margin.
#' @export
#'
#' @examples laplace_frechet_pit(y=frechet_laplace_pit(x=generate_Y(N=1)))
laplace_frechet_pit <- function(y) {
  if (y<0) { x <- (log(2)-y)^(-1) }
  else { x <- -(log(1-exp(-(y+log(2)))))^(-1) }
  return(x)
}

#' Transform from normal to Laplace margins
#'
#' @param x A number sampled from normal distribution.
#'
#' @return x transformed to Laplace margin.
#' @export
#'
#' @examples norm_laplace_pit(x=rnorm(1))
norm_laplace_pit <-  function(x) {
  if (pnorm(x)<0.5) { y <-log(2*pnorm(x)) }
  else { y <- -log(2*(1-pnorm(x))) }
  return(y)
}

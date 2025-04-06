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
  names(sims)[ncol(sims)] <- paste0("Y",ncol(sims)) # rename columns
  return(sims)
}

#' Sample from X|Y logistic distribution (neater formula)
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
link_log1 <- function(sims,dep=1/2) {
  Y <- sims %>% pull(-1)
  N <- length(Y)
  a <- dep
  # generate x
  to_opt <- function(x) {
    z <- x/y
    (  (  (1+z^(-1/a))^(a-1)*  exp(y^(-1)*( 1-(1+z^(-1/a))^a ) ) )-Unif)^2
  }
  x <- c()
  for (i in 1:N){
    Unif <- runif(1) # generate U
    y <- Y[i]
    x[i] <- optim(par=1,fn=to_opt,lower=0,upper=10^6,method="Brent")$par
  }
  sims <- sims %>% mutate(X=x)
  names(sims)[ncol(sims)] <- paste0("Y",ncol(sims)) # rename columns
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
  names(sims)[ncol(sims)] <- paste0("Y",ncol(sims)) # rename columns
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

#' Transform from Uniform to Laplace margins
#'
#' @param x A number sampled from U(0,1) distribution.
#'
#' @return x transformed to Laplace margin
#' @export
#'
#' @examples unif_laplace_pit(x=runif(1))
unif_laplace_pit <-  function(x) {
  if (x<0.5) { y <-log(2*x) }
  else { y <- -log(2*(1-x)) }
  return(y)
}

#' CDF of AGG margins for residuals
#'
#' @param x A number sampled from U(0,1) distribution.
#'
#' @return x transformed to AGG margins
#' @export
#'
#' @examples F_AGG(x=runif(1))
F_AGG <- function(x,theta) {
  mu <- theta[1]
  sigl <- theta[2]
  sigu <- theta[3]
  deltal <- theta[4]
  deltau <- theta[5]
  C_AGG <-  (sigl/deltal*gamma(1/deltal) + sigu/deltau*gamma(1/deltau)  )^(-1)
  if (x<mu) { y <- C_AGG*sigl/deltal*as.numeric(pracma::gammainc(x=((mu-x)/sigl)^deltal, a=1/deltal))[2] }
  else { y <- C_AGG*sigl/deltal *gamma(1/deltal) + C_AGG*sigu/deltau*as.numeric(pracma::gammainc(x=((x-mu)/sigu)^deltau, a=1/deltau)[1]) }
  return(y)
}

AGG_density <- function(x,theta) {
  mu <- theta[1]
  sigl <- theta[2]
  sigu <- theta[3]
  deltal <- theta[4]
  deltau <- theta[5]
  C_AGG <-  (sigl/deltal*gamma(1/deltal) + sigu/deltau*gamma(1/deltau)  )^(-1)
  y <- c()
  for (i in seq(x)) {
    if (x[i]<mu) {
      y[i] <- C_AGG*exp(-abs((mu-x[i])/sigl)^deltal)
    }
    else {y[i] <- C_AGG*exp(-((x[i]-mu)/sigu)^deltau)}
  }
  return(y)
}

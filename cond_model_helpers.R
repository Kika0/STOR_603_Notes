# generate from trivariate logistic EVD distribution with X-Y-Z dependence
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

# transform from Fréchet to Laplace
frechet_laplace_pit <- function(x) {
  if (exp(-1/x)<0.5) {
    y <-log(2*exp(-1/x))
  }
  else {
    y <--log(2*(1-exp(-1/x)))
  }
  return(y)
}

# transform back from Laplace to Fréchet margins
laplace_frechet_pit <- function(y) {
  if (y<0) {
    x <- (log(2)-y)^(-1)
  }
  else {
    x <- -(log(1-exp(-(y+log(2)))))^(-1)
  }
  return(x)
}

# conditioning on Y_1 being extreme to model Y_2 (conditional model in bivariate case)
# calculates MLE a_hat, b_hat, mu_hat and sig_hat
Y_2_likelihood <- function(theta,df=Y_given_1_extreme,given=1,sim=2) {
  a <- theta[1]
  b <- theta[2]
  mu <- theta[3]
  sig <- theta[4]
  Y1 <- df %>% select(starts_with("Y") & contains(as.character(given))) %>% pull()
  Y2 <- df %>% select(starts_with("Y") & contains(as.character(sim))) %>% pull()
  #lik <-  prod(1/(Y_1^b *sig)*exp(-(Y_2-a*Y_1-mu*Y_1^b)^2/(2*(Y_1^b*sig)^2)) )
  log_lik <- sum(-log(Y1^b *sig) + (-(Y2-a*Y_1-mu*Y1^b)^2/(2*(Y1^b*sig)^2))  )
  return(log_lik)
}
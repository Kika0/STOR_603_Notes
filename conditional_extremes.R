library(evd)
library(tidyverse)
library(latex2exp)
library(gridExtra)

generate_dep_X_Y_Y_Z <- function(N,dep=c(1/2,1/2)) {
  dat <- tibble(x1=numeric(),x2=numeric(),x3=numeric(),R=numeric(),u=numeric(),q=numeric(),w1=numeric(),w2=numeric(),w3=numeric())
  set.seed(12)
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
    z[i] <- optim(par=1,fn=to_opt)$par
  }
  sims <- data.frame(X_1=x,X_2=Y,X_3=z)
  return(sims)
}

# generate trivariate sample
sims <- generate_dep_X_Y_Y_Z(N=50000)

ggplot(sims %>% pivot_longer(everything())) + geom_density(aes(x=value),,stat="density") + facet_wrap(~name)

frechet_laplace_pit <- function(x) {
  if (exp(-1/x)<0.5) {
    y <-log(2*exp(-1/x))
  }
  else {
    y <--log(2*(1-exp(-1/x)))
  }
  return(y)
}

# PIT to Laplace
sims <- sims %>% mutate(Y_1=as.numeric(map_if(X_1,.p=function(x){exp(-1/x)<0.5},
                           .f=function(x) {log(2*exp(-1/x))},
                           .else=function(x){-log(2*(1-exp(-1/x)))}))) %>% 
  mutate(Y_2=as.numeric(map_if(X_2,.p=function(x){exp(-1/x)<0.5},
                    .f=function(x) {log(2*exp(-1/x))},
                    .else=function(x){-log(2*(1-exp(-1/x)))}))) %>%
  mutate(Y_3=as.numeric(map_if(X_3,.p=function(x){exp(-1/x)<0.5},
                    .f=function(x) {log(2*exp(-1/x))},
                    .else=function(x){-log(2*(1-exp(-1/x)))})))

ggplot(sims %>% select(Y_1,Y_2,Y_3) %>% pivot_longer(everything())) + geom_density(aes(x=value),,stat="density") + facet_wrap(~name)

grid.arrange(ggplot(sims) + geom_point(aes(x=Y_1,y=Y_2),alpha=0.5),
             ggplot(sims) + geom_point(aes(x=Y_2,y=Y_3),alpha=0.5),
             ggplot(sims) + geom_point(aes(x=Y_1,y=Y_3),alpha=0.5),ncol=3)

# filter for Y_1 being extreme
v <- 0.9
Y_given_1_extreme <- sims %>% filter(Y_1>quantile(Y_1,v))

Y_2_likelihood <- function(theta,df=Y_given_1_extreme) {
  a <- theta[1]
  b <- theta[2]
  mu <- theta[3]
  sig <- theta[4]
  Y_1 <- df$Y_1
  Y_2 <- df$Y_2
 #lik <-  prod(1/(Y_1^b *sig)*exp(-(Y_2-a*Y_1-mu*Y_1^b)^2/(2*(Y_1^b*sig)^2)) )
  log_lik <- sum(-log(Y_1^b *sig) + (-(Y_2-a*Y_1-mu*Y_1^b)^2/(2*(Y_1^b*sig)^2))  )
 return(log_lik)
}

optim(par=c(1,0,0,1),fn = Y_2_likelihood,df=Y_given_1_extreme,control = list(fnscale=-1))

# plot the values inferenced on

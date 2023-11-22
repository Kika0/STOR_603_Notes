library(evd)
library(tidyverse)
library(latex2exp)
library(gridExtra)

# set theme defaults to be black rectangle
theme_set(theme_bw())
theme_replace(
  panel.spacing = unit(2, "lines"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA) )

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
sims <- sims %>% mutate(Y_1=as.numeric(map(.x=X_1,.f=frechet_laplace_pit))) %>% 
  mutate(Y_2=as.numeric(map(.x=X_2,.f=frechet_laplace_pit))) %>%
  mutate(Y_3=as.numeric(map(.x=X_3,.f=frechet_laplace_pit)))

ggplot(sims %>% select(Y_1,Y_2,Y_3) %>% pivot_longer(everything())) + geom_density(aes(x=value),,stat="density") + facet_wrap(~name)

grid.arrange(ggplot(sims) + geom_point(aes(x=Y_1,y=Y_2),alpha=0.5),
             ggplot(sims) + geom_point(aes(x=Y_2,y=Y_3),alpha=0.5),
             ggplot(sims) + geom_point(aes(x=Y_1,y=Y_3),alpha=0.5),ncol=3)

# filter for Y_1 being extreme -----
v <- 0.9
Y_given_1_extreme <- sims %>% filter(Y_1>quantile(Y_1,v))
Y_not_1_extreme <- sims %>% filter(Y_1<quantile(Y_1,v))

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

opt <- optim(par=c(1,0,0,1),fn = Y_2_likelihood,df=Y_given_1_extreme,control = list(fnscale=-1))
a_hat <- opt$par[1]
b_hat <- opt$par[2]
mu_hat <- opt$par[3]
sig_hat <- opt$par[4]
Y_1 <- Y_given_1_extreme[,4]
# plot the values inferenced on
# generate from Normal distribution
N <- 50000
Y_2_sim <- rnorm(n=length(Y_1),mean = a_hat*Y_1 + Y_1^b_hat*mu_hat,sd=sig_hat*Y_1^b_hat )

Y_given_1_extreme <- Y_given_1_extreme %>% mutate(Y_2_sim=Y_2_sim)
ggplot(Y_given_1_extreme %>% select(Y_1,Y_2_sim,Y_2,Y_3) %>% pivot_longer(everything())) + geom_density(aes(x=value),,stat="density") + facet_wrap(~name)

# transform back to Fréchet margins
laplace_frechet_pit <- function(y) {
  if (y<0) {
    x <- (log(2)-y)^(-1)
  }
  else {
    x <- -(log(1-exp(-(y+log(2)))))^(-1)
  }
  return(x)
}

Y_2_sim <- bind_rows(Y_not_1_extreme %>% mutate(Y_2_sim=Y_2),Y_given_1_extreme) %>% 
  mutate(X_2_sim=rep(0,N))

X_2_simdf <- Y_2_sim %>% 
  mutate(X_2=as.numeric(map(.x=Y_2,.f=laplace_frechet_pit))) %>%
  mutate(X_2_sim=as.numeric(map(.x=Y_2_sim,.f=laplace_frechet_pit)))

X_2_simdf <- X_2_simdf %>% mutate(v=c(rep("below_threshold",45000),rep("above_threshold",5000)))

grid.arrange(ggplot(X_2_simdf) + geom_point(aes(x=X_1,y=X_2,col=v),alpha=0.5),
             ggplot(X_2_simdf) + geom_point(aes(x=X_2,y=X_3,col=v),alpha=0.5),
             ggplot(X_2_simdf) + geom_point(aes(x=X_1,y=X_2_sim,col=v),alpha=0.5),
             ggplot(X_2_simdf) + geom_point(aes(x=X_2_sim,y=X_3,col=v),alpha=0.5),ncol=2)

grid.arrange(ggplot(X_2_simdf) + geom_point(aes(x=Y_1,y=Y_2,col=v),alpha=0.5),
             ggplot(X_2_simdf) + geom_point(aes(x=Y_2,y=Y_3,col=v),alpha=0.5),
             ggplot(X_2_simdf) + geom_point(aes(x=Y_1,y=Y_2_sim,col=v),alpha=0.5),
             ggplot(X_2_simdf) + geom_point(aes(x=Y_2_sim,y=Y_3,col=v),alpha=0.5),ncol=2)
# plot only extremes
grid.arrange(ggplot(X_2_simdf %>% filter(v=="above_threshold")) + geom_point(aes(x=X_1,y=X_2),alpha=0.5),
             ggplot(X_2_simdf%>% filter(v=="above_threshold")) + geom_point(aes(x=X_2,y=X_3),alpha=0.5),
             ggplot(X_2_simdf%>% filter(v=="above_threshold")) + geom_point(aes(x=X_1,y=X_2_sim),alpha=0.5),
             ggplot(X_2_simdf%>% filter(v=="above_threshold")) + geom_point(aes(x=X_2_sim,y=X_3),alpha=0.5),ncol=2)

grid.arrange(ggplot(X_2_simdf %>% filter(v=="above_threshold")) + geom_point(aes(x=Y_1,y=Y_2),alpha=0.5),
             ggplot(X_2_simdf%>% filter(v=="above_threshold")) + geom_point(aes(x=Y_2,y=Y_3),alpha=0.5),
             ggplot(X_2_simdf%>% filter(v=="above_threshold")) + geom_point(aes(x=Y_1,y=Y_2_sim),alpha=0.5),
             ggplot(X_2_simdf%>% filter(v=="above_threshold")) + geom_point(aes(x=Y_2_sim,y=Y_3),alpha=0.5),ncol=2)

# repeat the simulation procedure for X_3 given X_1 is extreme
v <- 0.9
Y_given_1_extreme <- sims %>% filter(Y_1>quantile(Y_1,v))
Y_not_1_extreme <- sims %>% filter(Y_1<quantile(Y_1,v))

Y_3_likelihood <- function(theta,df=Y_given_1_extreme) {
  a <- theta[1]
  b <- theta[2]
  mu <- theta[3]
  sig <- theta[4]
  Y_1 <- df$Y_1
  Y_2 <- df$Y_3
  #lik <-  prod(1/(Y_1^b *sig)*exp(-(Y_2-a*Y_1-mu*Y_1^b)^2/(2*(Y_1^b*sig)^2)) )
  log_lik <- sum(-log(Y_1^b *sig) + (-(Y_2-a*Y_1-mu*Y_1^b)^2/(2*(Y_1^b*sig)^2))  )
  return(log_lik)
}

opt <- optim(par=c(1,0,0,1),fn = Y_3_likelihood,df=Y_given_1_extreme,control = list(fnscale=-1))
a_hat <- opt$par[1]
b_hat <- opt$par[2]
mu_hat <- opt$par[3]
sig_hat <- opt$par[4]
Y_1 <- Y_given_1_extreme[,4]
# plot the values inferenced on
# generate from Normal distribution
N <- 50000
Y_3_sim <- rnorm(n=length(Y_1),mean = a_hat*Y_1 + Y_1^b_hat*mu_hat,sd=sig_hat*Y_1^b_hat )

Y_given_1_extreme <- Y_given_1_extreme %>% mutate(Y_3_sim=Y_3_sim)
ggplot(Y_given_1_extreme %>% select(Y_1,Y_3_sim,Y_2,Y_3) %>% pivot_longer(everything())) + geom_density(aes(x=value),,stat="density") + facet_wrap(~name)

# transform back to Fréchet margins
laplace_frechet_pit <- function(y) {
  if (y<0) {
    x <- (log(2)-y)^(-1)
  }
  else {
    x <- -(log(1-exp(-(y+log(2)))))^(-1)
  }
  return(x)
}

Y_3_sim <- bind_rows(Y_not_1_extreme %>% mutate(Y_3_sim=Y_3),Y_given_1_extreme) %>% 
  mutate(X_3_sim=rep(0,N))

X_3_simdf <- Y_3_sim %>% 
  mutate(X_3=as.numeric(map(.x=Y_3,.f=laplace_frechet_pit))) %>%
  mutate(X_3_sim=as.numeric(map(.x=Y_3_sim,.f=laplace_frechet_pit)))

X_3_simdf <- X_3_simdf %>% mutate(v=c(rep("below_threshold",45000),rep("above_threshold",5000)))

grid.arrange(ggplot(X_3_simdf) + geom_point(aes(x=X_1,y=X_2,col=v),alpha=0.5),
             ggplot(X_3_simdf) + geom_point(aes(x=X_2,y=X_3,col=v),alpha=0.5),
             ggplot(X_3_simdf) + geom_point(aes(x=X_1,y=X_3,col=v),alpha=0.5),
             ggplot(X_3_simdf) + geom_point(aes(x=X_1,y=X_3_sim,col=v),alpha=0.5),ncol=2)

grid.arrange(ggplot(X_3_simdf) + geom_point(aes(x=Y_1,y=Y_2,col=v),alpha=0.5),
             ggplot(X_3_simdf) + geom_point(aes(x=Y_2,y=Y_3,col=v),alpha=0.5),
             ggplot(X_3_simdf) + geom_point(aes(x=Y_1,y=Y_3_sim,col=v),alpha=0.5),
             ggplot(X_3_simdf) + geom_point(aes(x=Y_2,y=Y_3_sim,col=v),alpha=0.5),ncol=2)
# plot only extremes
grid.arrange(ggplot(X_3_simdf %>% filter(v=="above_threshold")) + geom_point(aes(x=X_1,y=X_2),alpha=0.5),
             ggplot(X_3_simdf%>% filter(v=="above_threshold")) + geom_point(aes(x=X_2,y=X_3),alpha=0.5),
             ggplot(X_3_simdf%>% filter(v=="above_threshold")) + geom_point(aes(x=X_1,y=X_3_sim),alpha=0.5),
             ggplot(X_3_simdf%>% filter(v=="above_threshold")) + geom_point(aes(x=X_2,y=X_3_sim),alpha=0.5),ncol=2)

grid.arrange(ggplot(X_3_simdf %>% filter(v=="above_threshold")) + geom_point(aes(x=Y_1,y=Y_2),alpha=0.5),
             ggplot(X_3_simdf%>% filter(v=="above_threshold")) + geom_point(aes(x=Y_2,y=Y_3),alpha=0.5),
             ggplot(X_3_simdf%>% filter(v=="above_threshold")) + geom_point(aes(x=Y_1,y=Y_3_sim),alpha=0.5),
             ggplot(X_3_simdf%>% filter(v=="above_threshold")) + geom_point(aes(x=Y_2,y=Y_3_sim),alpha=0.5),ncol=2)

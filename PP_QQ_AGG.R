# libraries and functions ----
#library(evd)
library(tidyverse)
#library(latex2exp)
library(gridExtra)
library(pracma)
#library(viridis)
#library(MASS) #use dplyr::select to avoid function conflict
#library(xtable)
#library(gnorm)
file.sources = list.files(pattern="*helpers.R")
sapply(file.sources,source,.GlobalEnv)
# set theme defaults to be black rectangle
theme_set(theme_bw())
theme_replace(
  panel.spacing = unit(2, "lines"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA) )

# try Nrep macroreplications ----
N <- 50000
v <- 0.99
Nv <- N*(1-v)
u1 <- l1 <-  1:Nv/(Nv+1) # x-axis of PP plot
bf1 <- data.frame(x=1:Nv)
bf2 <- data.frame(x=1:Nv)
bf1num <- bf2num <- numeric()
Nrep <- 100
# store a,b,mu,sig
a <- b <- mu <- sig <- numeric()
for (i in 1:Nrep) {
  p1 <- p2 <- c()
  # sample data
  set.seed(i*12)
  sim2 <- generate_Y(N=N) %>% link_log(dep=1/2) %>%
    apply(c(1,2),FUN=frechet_laplace_pit) %>% as.data.frame()
  pe <- par_est(df=sim2,v=v,given=1,margin = "AGGsigdelta", method = "sequential2")
  # calculate p values for the observed residuals
  Z2p <- 1:Nv/(Nv+1)
  Y2 <- Z2p # observed residuals vector
  Z2 <- sort(as.numeric(as.data.frame(obs_res)[,1]))
  # append to a,b,mu,sig
  a <- append(a,pe$a[1])
  b <- append(b,pe$b[1])
  mu <- append(mu,pe$mu[1])
  sig <- append(sig,pe$sig[1])
  obs_res <- as.data.frame(observed_residuals(df = sim2,given = 1,v = v,a=pe$a[1],b=pe$b[1]))
  Z2p <- 1:Nv/(Nv+1)
  Y2 <- Z2p # observed residuals vector
  Y1 <- c()
  Z2sort <- sort(as.numeric(obs_res[,1])) # sorted observed residuals
  opt <- optim(fn=NLL_AGGsigdelta,x=Z2sort,par=c(0,1,1,1.2,1.8),control=list(maxit=2000),method = "BFGS")
  Y1 <- sapply(1:Nv,function(i){ F_AGG(x=Z2sort[i],theta = opt$par)})
  bf1 <- cbind(bf1,Y1)
  bf2 <- cbind(bf2,Y2)
}
bf1num <- as.numeric(unlist(bf1[,2:(Nrep+1)]))
bf2num <- as.numeric(unlist(bf2[,2:(Nrep+1)]))
Uup <- Ulow <- c()
for (i in 1:Nv) {
  Uup[i] <- quantile(bf1num[round(bf2num,5)==round(u1[i],5)],p=0.975)
  Ulow[i] <- quantile(bf1num[round(bf2num,5)==round(u1[i],5)],p=0.025)
}

# compare bootstrap and beta distribution for obtaining tolerance bounds
set.seed(11)
v <- 0.99
sim2 <- generate_Y(N=N) %>% link_log(dep=1/2) %>%
  apply(c(1,2),FUN=frechet_laplace_pit) %>% as.data.frame()
# fit PP plot to the observed residuals
pe <- par_est(df=sim2,v=v,given=1,margin = "AGGsigdelta", method = "two_step")
# calculate observed residuals
obs_res <- as.data.frame(observed_residuals(df = sim2,given = 1,v = v,a=pe$a[1],b=pe$b[1])) 
# calculate empirical p values for the observed residuals
Z2p <- (1:Nv)/(Nv+1)
# PP plot
Z2sort <- sort(as.numeric(obs_res[,1]))
opt <- optim(fn=NLL_AGGsigdelta,x=Z2sort,par=c(0,1,1,1.2,1.8),control=list(maxit=2000),method = "BFGS")
Z2fit <- sapply(1:Nv,function(i){ F_AGG(x=Z2sort[i],theta = opt$par)})

# compare also with more samples
names(bf1) <- names(bf2) <- c("remove",paste0("rep",1:Nrep))
tmp <- cbind(bf1[,2:(Nrep+1)] %>% pivot_longer(everything(),names_to = "samp",values_to="y"),
                 (bf2[,2:(Nrep+1)] %>% pivot_longer(everything(),names_to = "samp",values_to="x"))[,2]) %>% 
  mutate(samp=factor(samp))

# QQ plots
to_opt <- function(x,i) {
  return( (F_AGG(x, theta = c(mu,sigl,sigu,deltal,deltau))-(i/(Nv+1)))^2  )  
}

Z2Q <- sapply(1:Nv, function(i){optim(par=1,fn=to_opt,i=i)$par})

Qmin <- min(Z2sort,Z2Q)-0.2
Qmax <- max(Z2sort,Z2Q)+0.2

ggplot(data.frame(obs_res=Z2sort,Z=Z2Q)) + geom_point(aes(y=obs_res,x=Z)) + coord_fixed() + 
  geom_segment(data=data.frame(x1=Qmin,x2=Qmax,y1=Qmin,y2=Qmax),mapping=aes(x=x1,y=y1,xend=x2,yend=y2),linetype="dashed") 



### plotting PP and QQ ----
# comparison of bootstrap and beta tolerance bounds for PP plots
p1 <- PP_plot(observed = Z2p, simulated = Z2fit, tol_bounds = "bootstrap", title = "Bootstrap")
p2 <- PP_plot(observed = Z2p, simulated = Z2fit, tol_bounds ="beta_dist", title = "Beta distribution")
grid.arrange(p1,p2,ncol=2)

# macroreplications PP plots
Uup <- Ulow <- c()
for (i in 1:Nv) {
  Uup[i] <- quantile(bf1num[round(bf2num,5)==round(u1[i],5)],p=0.975)
  Ulow[i] <- quantile(bf1num[round(bf2num,5)==round(u1[i],5)],p=0.025)
}
p3 <- ggplot(tmp) + geom_point(aes(x=x,y=y,col=samp))+ theme(legend.position = "none")  + coord_fixed() + ggtitle("100 simulations") + xlab("Model") + ylab("Empirical")
p4 <- PP_plot(observed = Z2p, simulated = Z2fit, Uup = Uup, Ulow = Ulow, tol_bounds ="custom", title= "100 simulations tolerance bounds")
grid.arrange(p3,p4,ncol=2)

# comparison of bootstrap and beta tolerance bounds for QQ plots  
# calculate tolerance bounds for the beta distribution
Uup <- Ulow <- c()
Uup <- sapply(1:Nv, function(i){optim(fn=function(x,i) {
  return( (F_AGG(x, theta = c(mu,sigl,sigu,deltal,deltau))-qbeta(0.975, i, Nv+1-i))^2  )  
},par=1,i=i)$par})
Ulow <- sapply(1:Nv, function(i){optim(fn=function(x,i) {
  return( (F_AGG(x, theta = c(mu,sigl,sigu,deltal,deltau))-qbeta(0.025, i, Nv+1-i))^2  )  
},par=1,i=i)$par})
p5 <- QQ_plot(observed = Z2sort, simulated = Z2Q, tol_bounds = "bootstrap", title = "Bootstrap")
p6 <- QQ_plot(observed = Z2sort, simulated = Z2Q, Uup = Uup, Ulow = Ulow, tol_bounds ="custom", title = "Beta distribution")
grid.arrange(p5,p6,ncol=2) 

# plot also histograms of parameter estimates

# plot these for different threshold



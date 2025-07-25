# libraries and functions ----
#library(evd)
library(tidyverse)
library(gridExtra)
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

PP_QQ_save <- function(v=0.99,Nv=50) {
# try Nrep macroreplications ----
N <- round(Nv/(1-v))
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
  pe <- par_est(df=sim2,v=v,given=1,margin = "Normal", method = "sequential2")
  # calculate p values for the observed residuals
  Z2p <- 1:Nv/(Nv+1)
  Y2 <- Z2p # observed residuals vector
  # append to a,b,mu,sig
  a <- append(a,pe$a[1])
  b <- append(b,pe$b[1])
  mu <- append(mu,pe$mu[1])
  sig <- append(sig,pe$sig[1])
  obs_res <- as.data.frame(observed_residuals(df = sim2,given = 1,v = v,a=pe$a[1],b=pe$b[1]))
  Z2p <- 1:Nv/(Nv+1)
  Z2 <- sort(as.numeric(as.data.frame(obs_res)[,1]))
  Y2 <- Z2p # observed residuals vector
  Y1 <- c()
  Z2sort <- sort(as.numeric(obs_res[,1])) # sorted observed residuals
  opt <- optim(fn=NLL_AGG,x=Z2sort,par=c(mean(Z2sort),sd(Z2sort),sd(Z2sort),1.2,1.8),control=list(maxit=2000),method = "Nelder-Mead")
  Y1 <- pAGG(x=Z2sort,theta = opt$par)
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
set.seed(1234)
sim2 <- generate_Y(N=N) %>% link_log(dep=1/2) %>%
  apply(c(1,2),FUN=frechet_laplace_pit) %>% as.data.frame()
# fit PP plot to the obser <- ed residuals
pe <- par_est(df=sim2,v=v,given=1,margin = "Normal", method = "sequential2")
# calculate observed residuals
obs_res <- as.data.frame(observed_residuals(df = sim2,given = 1,v = v,a=pe$a[1],b=pe$b[1])) 
# calculate empirical p values for the observed residuals
Z2p <- (1:Nv)/(Nv+1)
# PP plot
Z2sort <- sort(as.numeric(obs_res[,1]))
opt <- optim(fn=NLL_AGG,x=Z2sort,par=c(mean(Z2sort),sd(Z2sort),sd(Z2sort),1.2,1.8),control=list(maxit=2000),method = "SANN")
Z2fit <- pAGG(x=Z2sort,theta = opt$par)

# compare different optim methods
methods_optim <- c("Nelder-Mead","BFGS","CG","SANN")
tr1 <- data.frame(y=as.numeric(),x=as.numeric(),Method_optim=as.character())
nllv <- c()
for (i in 1:4) {
  opt <- optim(fn=NLL_AGG,x=Z2sort,par=c(mean(Z2sort),1,1,1.2,1.8),control=list(maxit=2000),method = methods_optim[i])
  nllv[i] <- opt$value
  tr1 <- rbind(tr1,data.frame(y=AGG_density(x=seq(min(Z2sort),max(Z2sort),length.out=Nv*10),theta = opt$par),x=seq(min(Z2sort),max(Z2sort),length.out=Nv*10),Method_optim=methods_optim[i]))
}
pl <- ggplot(data.frame(x=Z2sort,y=Z2p)) + geom_density(aes(x=x))
pl1 <- pl + geom_line(data=tr1 %>% mutate(Method_optim=factor(Method_optim)),aes(x=x,y=y,col=Method_optim)) +
  xlab("Observed residuals fitted density") + ylab("Density") + ggtitle("Kernel smoothed and fitted density")
# ggsave(pl1,filename = paste0("method_optim_v",vi,".png"))
# compare also with more samples
names(bf1) <- names(bf2) <- c("remove",paste0("rep",1:Nrep))
tmp <- cbind(bf1[,2:(Nrep+1)] %>% pivot_longer(everything(),names_to = "samp",values_to="y"),
                 (bf2[,2:(Nrep+1)] %>% pivot_longer(everything(),names_to = "samp",values_to="x"))[,2]) %>% 
  mutate(samp=factor(samp))

# QQ plots
to_opt <- function(x,i,theta) {
  return( (pAGG(x, theta = theta)-Z2p[i])^2  )  
}

Z2Q <- qAGG(p=Z2p,theta = opt$par)
Qmin <- min(Z2sort,Z2Q)-0.2
Qmax <- max(Z2sort,Z2Q)+0.2

### plotting PP and QQ ----
# comparison of bootstrap and beta tolerance bounds for PP plots
p1 <- PP_plot(observed = Z2p, simulated = Z2fit, tol_bounds = "bootstrap", title = "Bootstrap")
p2 <- PP_plot(observed = Z2p, simulated = Z2fit, tol_bounds ="beta_dist", title = "Beta distribution")
ggsave(grid.arrange(p1,p2,ncol=2),filename =paste0("../Documents/PP_QQ_AGG/bootstrap_beta_PP",Nv,"_v",vi,".png"))

# macroreplications PP plots
Uup <- Ulow <- c()
for (i in 1:Nv) {
  Uup[i] <- quantile(bf1num[round(bf2num,5)==round(u1[i],5)],p=0.975)
  Ulow[i] <- quantile(bf1num[round(bf2num,5)==round(u1[i],5)],p=0.025)
}
p3 <- ggplot(tmp) + geom_point(aes(x=x,y=y,col=samp),size=0.1)+ theme(legend.position = "none")  + coord_fixed() + ggtitle("100 simulations") + xlab("Model") + ylab("Empirical")
p4 <- PP_plot(observed = Z2p, simulated = Z2fit, Uup = Uup, Ulow = Ulow, tol_bounds ="custom", title= paste0("100 simulations tolerance bounds"))
ggsave(grid.arrange(p3,p4,ncol=2), filename = paste0("../Documents/PP_QQ_AGG/100simul_PP",Nv,"_v",vi,"Nelder-Mead",".png"))

# comparison of bootstrap and beta tolerance bounds for QQ plots  
# calculate tolerance bounds for the beta distribution
Uup <- Ulow <- c()
Uup <- sapply(1:Nv, function(i){optim(fn=function(x,i) {
  return( (pAGG(x, theta = opt$par)-qbeta(0.975, i, Nv+1-i))^2  )  
},par=1,i=i)$par})
Ulow <- sapply(1:Nv, function(i){optim(fn=function(x,i) {
  return( (pAGG(x, theta = opt$par)-qbeta(0.025, i, Nv+1-i))^2  )  
},par=1,i=i)$par})
p5 <- QQ_plot(observed = Z2sort, simulated = Z2Q, tol_bounds = "bootstrap", title = "Bootstrap")
p6 <- QQ_plot(observed = Z2sort, simulated = Z2Q, Uup = Uup, Ulow = Ulow, tol_bounds ="custom", title = "Beta distribution")
ggsave(grid.arrange(p5,p6,ncol=2),filename = paste0("../Documents/PP_QQ_AGG/bootstrap_beta_QQ",Nv,"_v",vi,".png") )
print(nllv)
}

# plot these for different threshold and fixed number of exceedances -----------------------------
v <- c(0.8,0.9,0.99,0.999)
for (vi in v) {
  PP_QQ_save(v=vi,Nv=50)
}


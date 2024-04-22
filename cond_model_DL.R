# libraries and functions ----
library(evd)
library(tidyverse)
library(latex2exp)
library(gridExtra)
library(viridis)
library(MASS) #use dplyr::select to avoid function conflict
library(gnorm)
source("cond_model_helpers.R")

# set theme defaults to be black rectangle
theme_set(theme_bw())
theme_replace(
  panel.spacing = unit(2, "lines"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA) )

ggplot() + xlim(-10,20) + geom_function(fun=G_laplace,args=list(a=0.25),col="#C11432")+ geom_function(fun=G_laplace,args=list(a=1/2),col="#009ada")+ geom_function(fun=G_laplace,args=list(a=0.75),col="#66A64F")+ xlab(TeX("$z$")) + ylab(TeX("$G(z)$"))
ggplot() + xlim(-10,20) + geom_function(fun=g_laplace,args=list(a=0.25),col="#C11432")+ geom_function(fun=g_laplace,args=list(a=1/2),col="#009ada")+ geom_function(fun=g_laplace,args=list(a=0.75),col="#66A64F") + ylim(0,1) + xlab(TeX("$z$")) + ylab(TeX("$g(z)$"))

# generate trivariate sample ----
a_hat <- c()
b_hat <- c()
mu_hat <- c()
sig_hat <- c()
likl1 <- c()
likl2 <- c()
mu2 <- c()
sig2 <- c()
delta2 <- c()

for (i in 1:100) {
set.seed(123*i)
N <- 50000
sims <- generate_dep_X_Y_Y_Z(N=N,dep = c(0.9,1/2))
# PIT to Laplace
sims <- sims %>% mutate(Y1=as.numeric(map(.x=X_1,.f=frechet_laplace_pit))) %>% 
  mutate(Y2=as.numeric(map(.x=X_2,.f=frechet_laplace_pit))) %>%
  mutate(Y3=as.numeric(map(.x=X_3,.f=frechet_laplace_pit)))
df <- (sims %>% dplyr::select(starts_with("Y")))[,1:2]
Z <- c()
j <- given
cond_colours <- c("#C11432","#66A64F","#009ADA")
d <- ncol(df)
Y_given_1_extreme <- df %>% filter(df[,j]>quantile(df[,j],v))
res <- c(1:d)[-j]
# get initial parameters
init_opt <- optim(par=c(1,0,1), fn=Y_likelihood_initial,df=Y_given_1_extreme,given=j,sim=res[1],control = list(fnscale=-1))$par
init_par <- c(init_opt[1],0.2,init_opt[2],init_opt[3])
# optimise using the initial parameters
opt <- optim(par=init_par,fn = Y_likelihood,df=Y_given_1_extreme,given=j,sim=res[1],control = list(fnscale=-1))
a_hat <- append(a_hat,opt$par[1])
b_hat <- append(b_hat,opt$par[2])
mu_hat <- append(mu_hat,opt$par[3])
sig_hat <- c(sig_hat,opt$par[4])
likl1 <- append(likl1,opt$value)

# calculate observed residuals Z ----
Y1 <- Y_given_1_extreme[,j]
Y2 <- Y_given_1_extreme[,res[1]]
tmp_z2 <- (Y2-a_hat[i]*Y1)/(Y1^b_hat[i])
opt <- optim(DLLL2step,x=tmp_z2,par=c(0,1,1))
# optimise over the parameters
opt <- optim(DLLL,x=data.frame(Y1,Y2),par=c(0,1,1.5,0.8,0.8))
mu2 <- append(mu2,opt$par[1])
sig2 <- append(sig2,opt$par[2])
delta2 <-  append(delta2,opt$par[3])
likl2 <- append(likl2,opt$value)
}
 
# plot
plot(density(a_hat))
plot(density(b_hat))
plot(density(mu_hat))
plot(density(sig_hat))
plot(density(likl1))
plot(density(likl2))
plot(density(mu2))
plot(density(sig2))
plot(density(delta2))
# sample residual Z_star ----
Z_star <- rgnorm(n=1000,mu=mu,alpha=sig,beta=delta)

U <- runif(1000)
Y_1_gen <- -log(2*(1-0.999)) + rexp(1000)
Gen_Y_1 <- data.frame(Y_1=Y_1_gen)

# for each Y, generate a residual and calculate Y_2
Y_1 <- Gen_Y_1$Y_1
Y_2 <- a_hat*Y_1 + Y_1^b_hat *Z_star
Gen_Y_1 <- Gen_Y_1 %>% mutate(Y_2=Y_2) %>% mutate(sim=rep("model",1000))
names(Gen_Y_1) <- c(paste0("Y",j),paste0("Y",res[1]),"sim")
# generate Y_1 (extrapolate so above largest observed value)

#plot
Y1 <- Y_given_1_extreme[,j]
Y2 <- Y_given_1_extreme[,res[1]]
tmp <- data.frame(Y1,Y2) %>% mutate(sim=rep("data",50))
names(tmp) <- c(paste0("Y",j),paste0("Y",res[1]),"sim")
thres <- frechet_laplace_pit( qfrechet(0.999))
v <- 0.99
l <- min(Gen_Y_1 %>% dplyr::select(-sim),Y1,Y2)
u <- max(Gen_Y_1 %>% dplyr::select(-sim),Y1,Y2)
Gen_orig <- rbind(Gen_Y_1,tmp)
p1 <- ggplot(Gen_orig) +  annotate("rect",xmin=thres, ymin=thres, xmax=Inf,ymax=Inf, alpha=0.25,fill="#C11432") + geom_point(aes(x=Y1,y=Y2,col=sim),alpha=0.5) + geom_vline(xintercept=thres,col=cond_colours[j],linetype="dashed") +geom_abline(slope=0,intercept=thres,col=cond_colours[j],linetype="dashed") +
  scale_color_manual(values = c("data"="black","model" = cond_colours[j])) + xlab(TeX("$Y_1$")) +ylab(TeX("$Y_2$")) +
  xlim(c(l,u)) + ylim(c(l,u)) 

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

ggplot() + xlim(-10,20) + geom_function(fun=G_laplace,args=list(a=0.25),col="#C11432")+ geom_function(fun=G_laplace,args=list(a=1/2),col="#009ADA")+ geom_function(fun=G_laplace,args=list(a=0.75),col="#66A64F")+ xlab(TeX("$z$")) + ylab(TeX("$G(z)$"))
ggplot() + xlim(-10,20) +
  geom_function(fun=g_laplace,args=list(a=0.25),col="#C11432")+
  geom_function(fun=g_laplace,args=list(a=1/2),col="#009ADA")+
  geom_function(fun=g_laplace,args=list(a=0.75),col="#66A64F") + ylim(0,1) + xlab(TeX("$z$")) + ylab(TeX("$g(z)$"))

# compare different methods to model univariate residual ----
# two-step DL (Gaussian regression in the first step)
a_hat <- b_hat <- mu_hat <- sig_hat <- likl1a <- conv1a <- c()
mu2 <- sig2 <- delta2 <- likl1b <- conv1b <- c()
# one-step DL
a1 <- b1 <- mu1 <- sig1 <- delta1 <- likl2 <- conv2 <- c()
# asymmetric DL
a <- b <- mu <- sig <- deltal <- deltau <- likl3 <- conv3 <- c()

for (i in 1:200) {
set.seed(12*i)
N <- 50000
sims <- generate_dep_X_Y_Y_Z(N=N,dep = c(1/2,1/2))
# PIT to Laplace
sims <- sims %>% mutate(Y1=as.numeric(map(.x=X_1,.f=frechet_laplace_pit))) %>% 
  mutate(Y2=as.numeric(map(.x=X_2,.f=frechet_laplace_pit))) %>%
  mutate(Y3=as.numeric(map(.x=X_3,.f=frechet_laplace_pit)))
df <- (sims %>% dplyr::select(starts_with("Y")))[,1:2]
j <- 1
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
likl1a <- append(likl1a,opt$value)
conv1a <- append(conv1a,opt$convergence)

# calculate observed residuals Z ----
Y1 <- Y_given_1_extreme[,j]
Y2 <- Y_given_1_extreme[,res[1]]
tmp_z2 <- (Y2-opt$par[1]*Y1)/(Y1^opt$par[2])
opt <- optim(DLLL2step,x=tmp_z2,par=c(0,1,1),control = list(maxit=1000))
mu2 <- append(mu2,opt$par[1])
sig2 <- append(sig2,opt$par[2])
delta2 <-  append(delta2,opt$par[3])
likl1b <- append(likl1b,opt$value)
conv1b <- append(conv2,opt$convergence)
# optimise over all of the parameters
opt <- optim(fn=DLLL,x=data.frame(Y1,Y2),par=c(0,1,1.5,0.8,0.6),control=list(maxit=1000))
mu1 <- append(mu1,opt$par[1])
sig1 <- append(sig1,opt$par[2])
delta1 <-  append(delta1,opt$par[3])
a1 <- append(a1,opt$par[4])
b1 <-  append(b1,opt$par[5])
likl2 <- append(likl2,opt$value)
conv2 <- append(conv2,opt$convergence)
opt <- optim(fn=DLLLsk,x=data.frame(Y1,Y2),par=c(0,1,1.5,1.5,0.8,0.3),control=list(maxit=2000))
mu <- append(mu,opt$par[1])
sig <- append(sig,opt$par[2])
deltal <-  append(deltal,opt$par[3])
deltau <-  append(deltau,opt$par[4])
a <- append(a,opt$par[5])
b <-  append(b,opt$par[6])
likl3 <- append(likl3,opt$value)
conv3 <- append(conv3,opt$convergence)
}
c(sum(conv1a),sum(conv1b),sum(conv2),sum(conv3))

teststat1 <- 2 * (as.numeric(-likl1a)-as.numeric(likl2))
p.val1 <- pchisq(teststat1, df = 1, lower.tail = FALSE)
sum(p.val1<0.05)
teststat <- 2 * (as.numeric(likl2)-as.numeric(likl3))
p.val <- pchisq(teststat, df = 1, lower.tail = FALSE)
sum(p.val<0.05) 
# plot the estimated parameters
tmp_df <- data.frame(a_hat,b_hat,mu_hat,sig_hat,likl1a,likl1b,mu2=mu2,sig2=sig2,delta2,
                     a1,b1,mu1,sig1,delta1,likl2,a,b,mu,sig,deltal,deltau,likl3)
grid.arrange(ggplot(tmp_df %>% filter(conv1a==0,conv2==0,conv3==0)) + geom_density(aes(x=a_hat),col="#C11432",linetype="dashed")+ geom_density(aes(x=a1),col="#66A64F")+ geom_density(aes(x=a),col="#009ADA")+ xlab(TeX("$\\hat{\\alpha}$")),
ggplot(tmp_df) + geom_density(aes(x=b_hat),col="#C11432",linetype="dashed")+ geom_density(aes(x=b1),col="#66A64F")+ geom_density(aes(x=b),col="#009ADA")+ xlab(TeX("$\\hat{\\beta}$")),
ggplot(tmp_df) + geom_density(aes(x=mu_hat),col="#C11432",linetype="dashed") + geom_density(aes(x=mu2),col="#C11432")+ geom_density(aes(x=mu1),col="#66A64F")+ geom_density(aes(x=mu),col="#009ADA")+ xlab(TeX("$\\hat{\\mu}$")),
ggplot(tmp_df) + geom_density(aes(x=sig_hat),col="#C11432",linetype="dashed") + geom_density(aes(x=sig2),col="#C11432")+ geom_density(aes(x=sig1),col="#66A64F")+ geom_density(aes(x=sig),col="#009ADA")+ xlab(TeX("$\\hat{\\sigma}$")),
ggplot(tmp_df) + geom_density(aes(x=-likl1a),col="#C11432",linetype="dashed") + geom_density(aes(x=likl1b),col="#C11432",alpha=0.5)+ geom_density(aes(x=likl2),col="#66A64F",alpha=0.5)+ geom_density(aes(x=likl3),col="#009ADA")+xlab("negative log-likelihood"),
ggplot(tmp_df) + geom_density(aes(x=delta2),col="#C11432")+ geom_density(aes(x=delta1),col="#66A64F")+ geom_density(aes(x=deltal),col="#009ADA",size=1.5)+ geom_density(aes(x=deltau),col="#009ADA")+xlim(c(0,4))+ xlab(TeX("$\\hat{\\delta}$")),ncol=2)


ggplot(tmp_df%>% filter(conv1a==0,conv2==0,conv3==0) %>% mutate(diff12=-(likl2+likl1a),diff23=(likl2-likl3)))+
  geom_density(aes(x=diff12),linetype="dashed")+
  geom_density(aes(x=diff23)) + xlab("Difference in log-likelihood")

# plot the densities G(z) for one of the simulations
ggplot() + xlim(-10,10) + geom_density(mapping=aes(tmp_z2),linetype="dashed")+
  geom_function(fun=dgnorm,args=list(mu=mu1,alpha=sig1,beta=delta1),col="#66A64F")+
  geom_function(fun=dgnormsk,args=list(mu=mu,sig=sig,deltal=deltal,deltau=deltau),col="#009ADA") +
  geom_function(fun=g_laplace,args=list(a=1/2),col="black")+
  geom_function(fun=dgnorm,args=list(mu=mu2,alpha=sig2,beta=delta2),col="#C11432") + ylim(0,0.75) + xlab(TeX("$z$")) + ylab(TeX("$g(z)$"))

# sample residual Z_star ----
Z_star <- rgnorm(n=1000,mu=mu,alpha=sig,beta=delta)

# plot residuals
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

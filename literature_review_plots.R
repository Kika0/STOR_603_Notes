library(gridExtra)
library(tidyverse)
library(latex2exp) # latex expressions for plot labels
library(LaplacesDemon)
library(rvinecopulib)

theme_set(theme_bw())
theme_replace(
  panel.spacing = unit(2, "lines"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA) )

# set folder name to save plots in
folder_name <- "../Documents/literature_review_plots/"

# 1. chi measure for gaussian and logistic copula
# for loop for diff values of rho
U <- Chi <- Rho <- c()
edge <- 10^(-8)
u <- seq(0+edge,1-edge,length.out=10000)
rho_vec <- seq(-0.9,0.9,0.1)
for (i in 1:length(rho_vec)) {
  rho <- rho_vec[i]
  # calculate the copula using VGAM package
  c_uu <- rvinecopulib::pbicop(u=matrix(c(u,u),ncol=2),family="gaussian",parameters=rho)
  # calculate chi(u)
  chi <- 2- (log(c_uu))/(log(u))
  u <- u[chi>-10^6]
  chi <- chi[chi>-10^6]
  # subset by Frechet bounds (lower bound is sufficient)
  lower_bound1 <- chi>(2-(log(2*u-1))/(log(u)))
  lower_bound <- lower_bound1 %>% replace(is.na(.), TRUE)
  # add items to a list
  U <- append(U,u[lower_bound])
  Chi <- append(Chi,chi[lower_bound])
  Rho <- append(Rho,rep(rho_vec[i],sum(lower_bound)))
}
df <- data.frame(U,Chi,Rho=as.factor(Rho))
# subset to remove \chi(u)<-1 to match Coles' paper
df <- df[Chi>=-1,]
# ggplot requires prepared data before plotting
# ,col= guide_legend(title= TeX("$\rho$"))
p1 <- ggplot(df) + geom_line(aes(x=U,y=Chi,fill=Rho))+ xlab(TeX("$u$")) + ylab(TeX("$\\chi(u)$")) +  guides(colour = guide_legend(reverse=T)) + labs(fill=TeX("$\rho$")) + geom_segment(aes(x = 0, y = 1, xend = 1, yend = 1), linetype="dashed")
p1

# repeat with a logistic copula
# check which is the copula
library(evd)
library(LaplacesDemon)
tmp <- generate_Y(N=10000) %>% link_log(dep=1/2) 
tmp <- apply(tmp,MARGIN=c(2),FUN= function(i) { evd::pfrechet(i)})
rvinecopulib::bicop(data=tmp,family_set="parametric") # gumbel copula with 1/alpha parameter

U <- Chi <- Rho <- c()
edge <- 10^(-8)
u <- seq(0+edge,1-edge,length.out=10000)
rho_vec <- seq(0.1,0.9,0.1)
for (i in 1:length(rho_vec)) {
  rho <- rho_vec[i]
  # calculate the copula using VGAM package
  c_uu <- rvinecopulib::pbicop(u=matrix(c(u,u),ncol=2),family="gumbel",parameters=c(1/rho))
  plot
  # calculate chi(u)
  chi <- 2- (log(c_uu))/(log(u))
  u <- u[chi>-10^6]
  chi <- chi[chi>-10^6]
  # subset by Frechet bounds (lower bound is sufficient)
  lower_bound1 <- chi>(2-(log(2*u-1))/(log(u)))
  lower_bound <- lower_bound1 %>% replace(is.na(.), TRUE)
  # add items to a list
  U <- append(U,u[lower_bound])
  Chi <- append(Chi,chi[lower_bound])
  Rho <- append(Rho,rep(rho_vec[i],sum(lower_bound)))
}
df <- data.frame(U,Chi,Rho=as.factor(Rho))
# subset to remove \chi(u)<-1 to match Coles' paper
df <- df[Chi>=-1,]
# ggplot requires prepared data before plotting
# ,col= guide_legend(title= TeX("$\rho$"))
p2 <- ggplot(df) + geom_line(aes(x=U,y=Chi,fill=Rho))+ xlab(TeX("$u$")) + ylab(TeX("$\\chi(u)$")) +  guides(colour = guide_legend(reverse=T)) + labs(fill=TeX("$\rho$")) + geom_segment(aes(x = 0, y = 1, xend = 1, yend = 1), linetype="dashed")
p2

p <- grid.arrange(p1,p2,ncol=2)
ggsave(p,filename=paste0(folder_name,"chi_u_illust.png"),width=10,height=5)

# 2. illustration of different bivariate copulas
N <- 10000
size_point <- 0.1
# simulate random independent
tmpu <- data.frame("x"=runif(N),"y"=runif(N))
tmpl <- apply(tmpu,MARGIN=c(2),FUN=LaplacesDemon::qlaplace)
p1 <- ggplot(tmpu) + geom_point(aes(x=x,y=y),size=size_point)
p2 <- ggplot(tmpl) + geom_point(aes(x=x,y=y),size=size_point)
# simulate random Gaussian
tmpu <- rvinecopulib::rbicop(n=N,family=c("gaussian"),parameters=c(1/2)) %>% as.data.frame()
names(tmpu) <- c("x","y")
tmpl <- apply(tmpu,MARGIN=c(2),FUN=LaplacesDemon::qlaplace)
p3 <- ggplot(tmpu) + geom_point(aes(x=x,y=y),size=size_point)
p4 <- ggplot(tmpl) + geom_point(aes(x=x,y=y),size=size_point)
# simulate random logistic
tmpu <- rvinecopulib::rbicop(n=N,family=c("gumbel"),parameters=c(2)) %>% as.data.frame()
names(tmpu) <- c("x","y")
tmpl <- apply(tmpu,MARGIN=c(2),FUN=LaplacesDemon::qlaplace)
p5 <- ggplot(tmpu) + geom_point(aes(x=x,y=y),size=size_point)
p6 <- ggplot(tmpl) + geom_point(aes(x=x,y=y),size=size_point)
p <- grid.arrange(p1,p3,p5,p2,p4,p6,ncol=3)
ggsave(p,filename=paste0(folder_name,"copulas_illust.png"),width=10,height=7)

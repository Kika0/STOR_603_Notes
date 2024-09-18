# libraries and functions ----
library(VineCopula)
library(evd)
library(ismev)
library(tidyverse)
library(latex2exp)
library(gridExtra)
library(xtable) # for latex tables
library(GGally) # for ggpairs function
library(MASS) # use dplyr::select to avoid function conflict
library(texmex)
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

# simulate from the model ----
set.seed(11)
N <- 500
v <- 0.99
sims <- generate_Y(N=N) %>% link_log(dep=1/2) %>%
  link_log(dep=1/2) %>% link_log(dep=1/2) %>% link_log(dep=1/2) %>%
  apply(c(1,2),FUN=frechet_laplace_pit) %>% as.data.frame()
# repeat PP plots for summer and winter dataset

# fit a vine to all the data
fit <- RVineStructureSelect(sims %>% apply(c(2),FUN=row_number)/(nrow(sims)+1),
                             trunclevel = 3, indeptest = TRUE)
sim <- RVineSim(N=N,RVM=fit)
# transform back to laplace margins
# sl <- as.data.frame(sim %>% apply(c(1,2),FUN=unif_laplace_pit)) %>% as.data.frame()
sl <- generate_Y(N=N) %>% link_log(dep=1/2) %>%
  link_log(dep=1/2) %>% link_log(dep=1/2) %>% link_log(dep=1/2) %>%
  apply(c(1,2),FUN=frechet_laplace_pit) %>% as.data.frame()
grid.arrange(PP_plot(observed=as.numeric(sims[,1]),simulated = as.numeric(sl[,1]),title = TeX("$Y_1$")),
             PP_plot(observed=as.numeric(sims[,2]),simulated = as.numeric(sl[,2]),title = TeX("$Y_2$")),
             PP_plot(observed=as.numeric(sims[,3]),simulated = as.numeric(sl[,3]),title = TeX("$Y_3$")),
             PP_plot(observed=as.numeric(sims[,4]),simulated = as.numeric(sl[,4]),title = TeX("$Y_4$")),
             PP_plot(observed=as.numeric(sims[,5]),simulated = as.numeric(sl[,5]),title = TeX("$Y_5$")),ncol=2)

# compare all observed and simulated data for different subsets K
# K={1},{1,2},{1,2,3},{1,2,3,4},{1,2,3,4,5}
observed2 <- as.data.frame(sims %>% dplyr::select(1:2) %>% apply(c(2),FUN=row_number)/(nrow(sims)+1)) %>% apply(MARGIN=1,FUN = max)
# simulated2 <- as.data.frame(sim) %>% dplyr::select(1:2)  %>% apply(MARGIN=1,FUN = max)
simulated2 <- as.data.frame(sl %>% dplyr::select(1:2) %>% apply(c(2),FUN=row_number)/(nrow(sims)+1)) %>% apply(MARGIN=1,FUN = max)
observed3 <- as.data.frame(sims %>% dplyr::select(1:3) %>% apply(c(2),FUN=row_number)/(nrow(sims)+1)) %>% apply(MARGIN=1,FUN = max)
# simulated3 <-  as.data.frame(sim) %>% dplyr::select(1:3)  %>% apply(MARGIN=1,FUN = max)
simulated3 <- as.data.frame(sl %>% dplyr::select(1:3) %>% apply(c(2),FUN=row_number)/(nrow(sims)+1)) %>% apply(MARGIN=1,FUN = max)

observed4 <- as.data.frame(sims %>% dplyr::select(1:4) %>% apply(c(2),FUN=row_number)/(nrow(sims)+1)) %>% apply(MARGIN=1,FUN = max)
# simulated4 <-  as.data.frame(sim) %>% dplyr::select(1:4)  %>% apply(MARGIN=1,FUN = max)
simulated4 <-  as.data.frame(sl %>% dplyr::select(1:4) %>% apply(c(2),FUN=row_number)/(nrow(sims)+1)) %>% apply(MARGIN=1,FUN = max)
observed5 <- as.data.frame(sims %>% dplyr::select(1:5) %>% apply(c(2),FUN=row_number)/(nrow(sims)+1)) %>% apply(MARGIN=1,FUN = max)
# simulated5 <-  as.data.frame(sim) %>% dplyr::select(1:5)  %>% apply(MARGIN=1,FUN = max)
simulated5 <-  as.data.frame(sl %>% dplyr::select(1:5) %>% apply(c(2),FUN=row_number)/(nrow(sims)+1)) %>% apply(MARGIN=1,FUN = max)

datmax <- grid.arrange(PP_plot(observed = observed2,simulated = simulated2, title = TeX("$K=\\{1,2\\}$")),
                       PP_plot(observed = observed3,simulated = simulated3, title = TeX("$K=\\{1,2,3\\}$")),
                       PP_plot(observed = observed4,simulated = simulated4, title = TeX("$K=\\{1,2,3,4\\}$")),
                       PP_plot(observed = observed5,simulated = simulated5, title = TeX("$K=\\{1,2,3,4,5\\}$")),ncol=2)

# explore residuals transformed to uniform margins
# ggpairs((observed_residuals(df = sims,given = 1,v = 0.99) %>% apply(c(2),FUN=row_number))/(n_v+1))
# transform to Uniform margins and fit a vine
j <- 1
fit3 <- RVineStructureSelect((observed_residuals(df = sims,given = j,v = v) %>% apply(c(2),FUN=row_number))/(nrow(sims)*(1-v)+1),
                        trunclevel = 3, indeptest = TRUE)
fit3
# fit2 <- RVineStructureSelect((observed_residuals(df = sims,given = 1,v = v) %>% apply(c(2),FUN=row_number))/(nrow(sims)*(1-v)+1),
#                              trunclevel = 2, indeptest = FALSE)
# fit2
# fit1 <- RVineStructureSelect((observed_residuals(df = sims,given = 1,v = v) %>% apply(c(2),FUN=row_number))/(nrow(sims)*(1-v)+1),
#                              trunclevel = 1, indeptest = FALSE)
# fit1
# RVineClarkeTest(data=(observed_residuals(df = sims,given = 2,v = 0.99) %>% apply(c(2),FUN=row_number))/(nrow(sims)*(1-v)+1),
#                 RVM1 = fit1,RVM2 = fit2)
# 
# plot(fit3,edge.labels = "family")
# simulate from the copula
N_sim <- N*(1-v)
Zsim <- RVineSim(N=N_sim,RVM=fit3)
# transform back residuals to original margins
# can use kernel smoothed distribution as initial step
obs_res <- observed_residuals(df = sims,given = 1,v = v) 
to_opt <- function(z) {
  return( (mean(pnorm((z-obs_res[,k] %>% pull())/density(obs_res[,k] %>% pull())$bw)) - Zsim[i,k])^2)
}
Z <- Zsim
for (i in 1:nrow(Zsim)) {
  for (k in 1:ncol(Zsim)) {
    Z[i,k] <- optim(fn=to_opt,par=1)$par
  }
}
pairs(Z)
ggpairs(Z)
pairs(obs_res)
obsr <- (obs_res %>% 
    apply(c(2),FUN=row_number)) 
for (i in 1:nrow(obsr)) {
  for (k in 1:ncol(obsr)) {
    obsr[i,k] <- obsr[i,k]/(nrow(sims)*(1-v)+1)
  }
}

obsr %>% 
  ggpairs()

rbind(obs_res %>% as.data.frame() %>% mutate(res=rep("data",N*(1-v))),
      Z %>% as.data.frame() %>%
        mutate(res=rep("model",N*(1-v)))) %>% 
  ggpairs(columns = 1:4,ggplot2::aes(color=res,alpha=0.5), upper = list(continuous = wrap("cor", size = 2.5))) +
  scale_color_manual(values = c("data"="black","model" = "#C11432")) + scale_fill_manual(values = c("data"="black","model" = "#C11432"))

# simulate
v_sim <- v
Z_star <- as.data.frame(Z)
U <- runif(N_sim)
Y1_gen <- -log(2*(1-v_sim)) + rexp(N_sim)
Gen_Y1 <- data.frame(Y1=Y1_gen)

# for each Y, generate a residual and calculate Y_2
Y1 <- Gen_Y1$Y1
Y2 <- a_hat*Y1 + Y1^b_hat *Z_star[,1]
Y3 <-  a_hat*Y1 + Y1^b_hat *Z_star[,2]
Y4 <- a_hat*Y1 + Y1^b_hat *Z_star[,3]
Y5 <-  a_hat*Y1 + Y1^b_hat *Z_star[,4]
res <- c(1:d)[-j]
Gen_Y1 <- Gen_Y1 %>% mutate(Y2=Y2,Y3=Y3,Y4=Y4,Y5=Y5) %>% mutate(sim="model")
names(Gen_Y1) <- c(paste0("Y",j),paste0("Y",res[1]),paste0("Y",res[2]),paste0("Y",res[3]),paste0("Y",res[4]),"sim")
# generate Y_1 (extrapolate so above largest observed value)

#plot
d <- ncol(sims)
Y_given_1_extreme <- sims %>% filter(sims[,j]>quantile(sims[,j],v))
tmp <- data.frame(Y_given_1_extreme[,j],
                Y_given_1_extreme[,res[1]],
                Y_given_1_extreme[,res[2]],
                Y_given_1_extreme[,res[3]],
                Y_given_1_extreme[,res[4]]) %>%
  mutate(sim=rep("data",nrow(Y_given_1_extreme)))
names(tmp) <- c(paste0("Y",j),paste0("Y",res[1]),paste0("Y",res[2]),paste0("Y",res[3]),paste0("Y",res[4]),"sim")
Gen_orig <- rbind(Gen_Y1,tmp)
ggpairs(Gen_orig,columns = 1:5,ggplot2::aes(color=sim,alpha=0.5), upper = list(continuous = wrap("cor", size = 2.5))) +
  scale_color_manual(values = c("data"="black","model" = "#C11432")) + scale_fill_manual(values = c("data"="black","model" = "#C11432"))

# model diagnostics:PP plot ----
# start with observed and simulated residuals
tmpz <- rbind(obs_res %>% as.data.frame() %>% mutate(res=rep("data",N*(1-v))),
      Z %>% as.data.frame() %>%
        mutate(res=rep("model",N*(1-v))))
observed1 <- tmpz %>% filter(res=="data") %>% dplyr::select(1) %>% pull()
simulated1 <- tmpz %>% filter(res=="model") %>% dplyr::select(1) %>% pull()
observed2 <- tmpz %>% filter(res=="data") %>% dplyr::select(2) %>% pull()
simulated2 <- tmpz %>% filter(res=="model") %>% dplyr::select(2) %>% pull()
observed3 <- tmpz %>% filter(res=="data") %>% dplyr::select(3) %>% pull()
simulated3 <- tmpz %>% filter(res=="model") %>% dplyr::select(3) %>% pull()
observed4 <- tmpz %>% filter(res=="data") %>% dplyr::select(4) %>% pull()
simulated4 <- tmpz %>% filter(res=="model") %>% dplyr::select(4) %>% pull()

reseach <- grid.arrange(PP_plot(observed = observed1,simulated = simulated1, title = TeX("$Z_2$")),
             PP_plot(observed = observed2,simulated = simulated2, title = TeX("$Z_3$")),
             PP_plot(observed = observed3,simulated = simulated3, title = TeX("$Z_4$")),
             PP_plot(observed = observed4,simulated = simulated4, title = TeX("$Z_5$")),ncol=2)
ggsave(filename="plots/reseach.png",reseach,width=10,height=10)
# PP_plot(observed=rnorm(500),simulated=rnorm(500))
observed <- Gen_orig %>% filter(sim=="data") %>% dplyr::select(1) %>% pull()
simulated <- Gen_orig %>% filter(sim=="model") %>% dplyr::select(1) %>% pull()
PP_plot(observed = observed,simulated = simulated)
# compare all observed and simulated residuals for different subsets K
# K={1},{1,2},{1,2,3},{1,2,3,4}
observed2 <- as.data.frame(obs_res %>% dplyr::select(1:2) %>% apply(c(2),FUN=row_number)/(nrow(sims)*(1-v)+1)) %>% apply(MARGIN=1,FUN = max)
simulated2 <- as.data.frame(Zsim) %>% dplyr::select(1:2) %>% apply(MARGIN=1,FUN=max) 
observed3 <- as.data.frame(obs_res %>% dplyr::select(1:3) %>% apply(c(2),FUN=row_number)/(nrow(sims)*(1-v)+1)) %>% apply(MARGIN=1,FUN = max)
simulated3 <- as.data.frame(Zsim) %>% dplyr::select(1:3) %>% apply(MARGIN=1,FUN=max) 
observed4 <- as.data.frame(obs_res %>% dplyr::select(1:4) %>% apply(c(2),FUN=row_number)/(nrow(sims)*(1-v)+1)) %>% apply(MARGIN=1,FUN = max)
simulated4 <- as.data.frame(Zsim) %>% dplyr::select(1:4) %>% apply(MARGIN=1,FUN=max) 

resmax <- grid.arrange(PP_plot(observed = observed2,simulated = simulated2, title = TeX("$K=\\{1,2\\}$")),
             PP_plot(observed = observed3,simulated = simulated3, title = TeX("$K=\\{1,2,3\\}$")),
             PP_plot(observed = observed4,simulated = simulated4, title = TeX("$K=\\{1,2,3,4\\}$")),ncol=3)

ggsave(filename = "plots/resmax.png",resmax,width=10,height=3)

# compare observed and simulated data for each variable
observed1 <- Gen_orig %>% filter(sim=="data") %>% dplyr::select(1) %>% pull()
simulated1 <- Gen_orig %>% filter(sim=="model") %>% dplyr::select(1) %>% pull()
observed2 <- Gen_orig %>% filter(sim=="data") %>% dplyr::select(2) %>% pull()
simulated2 <- Gen_orig %>% filter(sim=="model") %>% dplyr::select(2) %>% pull()
observed3 <- Gen_orig %>% filter(sim=="data") %>% dplyr::select(3) %>% pull()
simulated3 <- Gen_orig %>% filter(sim=="model") %>% dplyr::select(3) %>% pull()
observed4 <- Gen_orig %>% filter(sim=="data") %>% dplyr::select(4) %>% pull()
simulated4 <- Gen_orig %>% filter(sim=="model") %>% dplyr::select(4) %>% pull()
observed5 <- Gen_orig %>% filter(sim=="data") %>% dplyr::select(5) %>% pull()
simulated5 <- Gen_orig %>% filter(sim=="model") %>% dplyr::select(5) %>% pull()
dateach <- grid.arrange(PP_plot(observed = observed1,simulated = simulated1, title = TeX("$Y_1$")),
             PP_plot(observed = observed2,simulated = simulated2, title = TeX("$Y_2$")),
             PP_plot(observed = observed3,simulated = simulated3, title = TeX("$Y_3$")),
             PP_plot(observed = observed4,simulated = simulated4, title = TeX("$Y_4$")),
             PP_plot(observed = observed5,simulated = simulated5, title = TeX("$Y_5$")),ncol=3)
ggsave(filename = "plots/dateach.png",dateach,width=10,height=7)
# compare all observed and simulated data for different subsets K
# K={1},{1,2},{1,2,3},{1,2,3,4},{1,2,3,4,5}
observed2 <- as.data.frame(Gen_orig %>% filter(sim=="data") %>% dplyr::select(1:2) %>% apply(c(2),FUN=row_number)/(nrow(sims)*(1-v)+1)) %>% apply(MARGIN=1,FUN = max)
simulated2 <- as.data.frame(Gen_orig %>% filter(sim=="model") %>% dplyr::select(1:2) %>% apply(c(2),FUN=row_number)/(nrow(sims)*(1-v)+1)) %>% apply(MARGIN=1,FUN = max)
observed3 <- as.data.frame(Gen_orig %>% filter(sim=="data") %>% dplyr::select(1:3) %>% apply(c(2),FUN=row_number)/(nrow(sims)*(1-v)+1)) %>% apply(MARGIN=1,FUN = max)
simulated3 <- as.data.frame(Gen_orig %>% filter(sim=="model") %>% dplyr::select(1:3) %>% apply(c(2),FUN=row_number)/(nrow(sims)*(1-v)+1)) %>% apply(MARGIN=1,FUN = max)
observed4 <- as.data.frame(Gen_orig %>% filter(sim=="data") %>% dplyr::select(1:4) %>% apply(c(2),FUN=row_number)/(nrow(sims)*(1-v)+1)) %>% apply(MARGIN=1,FUN = max)
simulated4 <- as.data.frame(Gen_orig %>% filter(sim=="model") %>% dplyr::select(1:4) %>% apply(c(2),FUN=row_number)/(nrow(sims)*(1-v)+1)) %>% apply(MARGIN=1,FUN = max)
observed5 <- as.data.frame(Gen_orig %>% filter(sim=="data") %>% dplyr::select(1:5) %>% apply(c(2),FUN=row_number)/(nrow(sims)*(1-v)+1)) %>% apply(MARGIN=1,FUN = max)
simulated5 <- as.data.frame(Gen_orig %>% filter(sim=="model") %>% dplyr::select(1:5) %>% apply(c(2),FUN=row_number)/(nrow(sims)*(1-v)+1)) %>% apply(MARGIN=1,FUN = max)

datmax <- grid.arrange(PP_plot(observed = observed2,simulated = simulated2, title = TeX("$K=\\{1,2\\}$")),
                       PP_plot(observed = observed3,simulated = simulated3, title = TeX("$K=\\{1,2,3\\}$")),
                       PP_plot(observed = observed4,simulated = simulated4, title = TeX("$K=\\{1,2,3,4\\}$")),
                       PP_plot(observed = observed5,simulated = simulated5, title = TeX("$K=\\{1,2,3,4,5\\}$")),ncol=2)

ggsave(filename = "plots/datmax.png",resmax,width=10,height=10)
# for loop to condition on each variable ----
v <- 0.99 # threshold for conditioning
v_sim <- 0.99 # threshold for simulation
N_sim <- 500 # number of observations to simulate
p_est <- c() # numeric of estimated probabilities of all extreme
for (l in 1:5) {
fit3 <- RVineStructureSelect((observed_residuals(df = sims,given = l,v = v) %>% apply(c(2),FUN=row_number))/(nrow(sims)*(1-v)+1),
                             trunclevel = 3, indeptest = FALSE)
Zsim <- RVineSim(N=N_sim,RVM=fit3)
res <- c(1:d)[-l]
# transform back residuals to original margins
# can use kernel smoothed distribution as initial step
obs_res <- observed_residuals(df = sims,given = l,v = v) 
to_opt <- function(z) {
  return( (mean(pnorm((z-obs_res[,k] %>% pull())/density(obs_res[,k] %>% pull())$bw)) - Zsim[i,k])^2)
}
Z <- Zsim
for (i in 1:nrow(Zsim)) {
  for (k in 1:ncol(Zsim)) {
    Z[i,k] <- optim(fn=to_opt,par=1)$par
  }
}
# plot observed residuals
p <- rbind(
      Z %>% as.data.frame() %>%
        mutate(res=rep("model",N_sim)),
      obs_res %>% as.data.frame() %>%
        mutate(res=rep("data",nrow(sims)*(1-v)))) %>% 
ggpairs(columns = 1:4,ggplot2::aes(color=res,alpha=0.5), upper = list(continuous = wrap("cor", size = 2.5))) +
  scale_color_manual(values = c("data"="black","model" = "#C11432")) + scale_fill_manual(values = c("data"="black","model" = "#C11432"))
#ggsave(p,filename=(paste0("plots/giv",l,"obsmodelres.pdf")),device="pdf")
# simulate
Z_star <- as.data.frame(Z)
Y_1_gen <- -log(2*(1-v_sim)) + rexp(N_sim)
Gen_Y_1 <- data.frame(Y1=Y_1_gen)
# for each Y, generate a residual and calculate Y_2
Y1 <- Gen_Y_1$Y1
Y2 <- a_hat*Y1 + Y1^b_hat *Z_star[,1]
Y3 <-  a_hat*Y1 + Y1^b_hat *Z_star[,2]
Y4 <- a_hat*Y1 + Y1^b_hat *Z_star[,3]
Y5 <-  a_hat*Y1 + Y1^b_hat *Z_star[,4]

Gen_Y_1 <- Gen_Y_1 %>% mutate(Y2=Y2,Y3=Y3,Y4=Y4,Y5=Y5) %>% mutate(sim=rep("model",nrow(Z_star)))
names(Gen_Y_1) <- c(paste0("Y",l),paste0("Y",res[1]),paste0("Y",res[2]),paste0("Y",res[3]),paste0("Y",res[4]),"sim")
# generate Y_1 (extrapolate so above largest observed value)

#plot
Y_given_1_extreme <- sims %>% filter(sims[,l]>quantile(sims[,l],v))
Y1 <- Y_given_1_extreme[,l]
Y2 <- Y_given_1_extreme[,res[1]]
Y3 <- Y_given_1_extreme[,res[2]]
Y4 <- Y_given_1_extreme[,res[3]]
Y5 <- Y_given_1_extreme[,res[4]]
tmp <- data.frame(Y1,Y2,Y3,Y4,Y5) %>% mutate(sim="data")
names(tmp) <- c(paste0("Y",l),paste0("Y",res[1]),paste0("Y",res[2]),paste0("Y",res[3]),paste0("Y",res[4]),"sim")
# l <- min(Gen_Y_1 %>% dplyr::select(-sim),Y1,Y2,Y3,Y4,Y5)
# u <- max(Gen_Y_1 %>% dplyr::select(-sim),Y1,Y2,Y3,Y4,Y5)
Gen_orig <- rbind(Gen_Y_1,tmp)
p <- ggpairs(Gen_orig,columns = 1:5,ggplot2::aes(color=sim,alpha=0.5), upper = list(continuous = wrap("cor", size = 2.5))) +
  scale_color_manual(values = c("data"="black","model" = "#C11432")) + scale_fill_manual(values = c("data"="black","model" = "#C11432"))
#ggsave(p,filename=(paste0("plots/cond",l,".pdf")),device="pdf")
# calculate probability of all extreme
vL_sim <- frechet_laplace_pit(qfrechet(v_sim))
p_est[l] <- Gen_Y_1 %>% filter(Y1>vL_sim,Y2>vL_sim,Y3>vL_sim,Y4>vL_sim,Y5>vL_sim) %>% 
  nrow()/N_sim*(1-v_sim)
tmpz <- rbind(obs_res %>% as.data.frame() %>% mutate(res=rep("data",N*(1-v))),
              Z %>% as.data.frame() %>%
                mutate(res=rep("model",N*(1-v))))
observed1 <- tmpz %>% filter(res=="data") %>% dplyr::select(1) %>% pull()
simulated1 <- tmpz %>% filter(res=="model") %>% dplyr::select(1) %>% pull()
observed2 <- tmpz %>% filter(res=="data") %>% dplyr::select(2) %>% pull()
simulated2 <- tmpz %>% filter(res=="model") %>% dplyr::select(2) %>% pull()
observed3 <- tmpz %>% filter(res=="data") %>% dplyr::select(3) %>% pull()
simulated3 <- tmpz %>% filter(res=="model") %>% dplyr::select(3) %>% pull()
observed4 <- tmpz %>% filter(res=="data") %>% dplyr::select(4) %>% pull()
simulated4 <- tmpz %>% filter(res=="model") %>% dplyr::select(4) %>% pull()

reseach <- grid.arrange(PP_plot(observed = observed1,simulated = simulated1, title = TeX("$Z_2$")),
                        PP_plot(observed = observed2,simulated = simulated2, title = TeX("$Z_3$")),
                        PP_plot(observed = observed3,simulated = simulated3, title = TeX("$Z_4$")),
                        PP_plot(observed = observed4,simulated = simulated4, title = TeX("$Z_5$")),ncol=2)
ggsave(filename=paste0("plots/reseachcond",l,"Nsim",N_sim,".png"),reseach,width=10,height=10)
# PP_plot(observed=rnorm(500),simulated=rnorm(500))
observed <- Gen_orig %>% filter(sim=="data") %>% dplyr::select(1) %>% pull()
simulated <- Gen_orig %>% filter(sim=="model") %>% dplyr::select(1) %>% pull()
PP_plot(observed = observed,simulated = simulated)
# compare all observed and simulated residuals for different subsets K
# K={1},{1,2},{1,2,3},{1,2,3,4}
observed2 <- as.data.frame(obs_res %>% dplyr::select(1:2) %>% apply(c(2),FUN=row_number)/(nrow(sims)*(1-v)+1)) %>% apply(MARGIN=1,FUN = max)
simulated2 <- as.data.frame(Zsim) %>% dplyr::select(1:2) %>% apply(MARGIN=1,FUN=max) 
observed3 <- as.data.frame(obs_res %>% dplyr::select(1:3) %>% apply(c(2),FUN=row_number)/(nrow(sims)*(1-v)+1)) %>% apply(MARGIN=1,FUN = max)
simulated3 <- as.data.frame(Zsim) %>% dplyr::select(1:3) %>% apply(MARGIN=1,FUN=max) 
observed4 <- as.data.frame(obs_res %>% dplyr::select(1:4) %>% apply(c(2),FUN=row_number)/(nrow(sims)*(1-v)+1)) %>% apply(MARGIN=1,FUN = max)
simulated4 <- as.data.frame(Zsim) %>% dplyr::select(1:4) %>% apply(MARGIN=1,FUN=max) 

resmax <- grid.arrange(PP_plot(observed = observed2,simulated = simulated2, title = TeX("$K=\\{1,2\\}$")),
                       PP_plot(observed = observed3,simulated = simulated3, title = TeX("$K=\\{1,2,3\\}$")),
                       PP_plot(observed = observed4,simulated = simulated4, title = TeX("$K=\\{1,2,3,4\\}$")),ncol=3)

ggsave(filename = paste0("plots/resmaxcond",l,"Nsim",N_sim,".png"),resmax,width=10,height=3)

# compare observed and simulated data for each variable
observed1 <- Gen_orig %>% filter(sim=="data") %>% dplyr::select(1) %>% pull()
simulated1 <- Gen_orig %>% filter(sim=="model") %>% dplyr::select(1) %>% pull()
observed2 <- Gen_orig %>% filter(sim=="data") %>% dplyr::select(2) %>% pull()
simulated2 <- Gen_orig %>% filter(sim=="model") %>% dplyr::select(2) %>% pull()
observed3 <- Gen_orig %>% filter(sim=="data") %>% dplyr::select(3) %>% pull()
simulated3 <- Gen_orig %>% filter(sim=="model") %>% dplyr::select(3) %>% pull()
observed4 <- Gen_orig %>% filter(sim=="data") %>% dplyr::select(4) %>% pull()
simulated4 <- Gen_orig %>% filter(sim=="model") %>% dplyr::select(4) %>% pull()
observed5 <- Gen_orig %>% filter(sim=="data") %>% dplyr::select(5) %>% pull()
simulated5 <- Gen_orig %>% filter(sim=="model") %>% dplyr::select(5) %>% pull()
dateach <- grid.arrange(PP_plot(observed = observed1,simulated = simulated1, title = TeX("$Y_1$")),
                        PP_plot(observed = observed2,simulated = simulated2, title = TeX("$Y_2$")),
                        PP_plot(observed = observed3,simulated = simulated3, title = TeX("$Y_3$")),
                        PP_plot(observed = observed4,simulated = simulated4, title = TeX("$Y_4$")),
                        PP_plot(observed = observed5,simulated = simulated5, title = TeX("$Y_5$")),ncol=3)
ggsave(filename = paste0("plots/dateachcond",l,"Nsim",N_sim,".png"),dateach,width=10,height=7)
# compare all observed and simulated data for different subsets K
# K={1},{1,2},{1,2,3},{1,2,3,4},{1,2,3,4,5}
observed2 <- as.data.frame(Gen_orig %>% filter(sim=="data") %>% dplyr::select(1:2) %>% apply(c(2),FUN=row_number)/(nrow(sims)*(1-v)+1)) %>% apply(MARGIN=1,FUN = max)
simulated2 <- as.data.frame(Gen_orig %>% filter(sim=="model") %>% dplyr::select(1:2) %>% apply(c(2),FUN=row_number)/(nrow(sims)*(1-v)+1)) %>% apply(MARGIN=1,FUN = max)
observed3 <- as.data.frame(Gen_orig %>% filter(sim=="data") %>% dplyr::select(1:3) %>% apply(c(2),FUN=row_number)/(nrow(sims)*(1-v)+1)) %>% apply(MARGIN=1,FUN = max)
simulated3 <- as.data.frame(Gen_orig %>% filter(sim=="model") %>% dplyr::select(1:3) %>% apply(c(2),FUN=row_number)/(nrow(sims)*(1-v)+1)) %>% apply(MARGIN=1,FUN = max)
observed4 <- as.data.frame(Gen_orig %>% filter(sim=="data") %>% dplyr::select(1:4) %>% apply(c(2),FUN=row_number)/(nrow(sims)*(1-v)+1)) %>% apply(MARGIN=1,FUN = max)
simulated4 <- as.data.frame(Gen_orig %>% filter(sim=="model") %>% dplyr::select(1:4) %>% apply(c(2),FUN=row_number)/(nrow(sims)*(1-v)+1)) %>% apply(MARGIN=1,FUN = max)
observed5 <- as.data.frame(Gen_orig %>% filter(sim=="data") %>% dplyr::select(1:5) %>% apply(c(2),FUN=row_number)/(nrow(sims)*(1-v)+1)) %>% apply(MARGIN=1,FUN = max)
simulated5 <- as.data.frame(Gen_orig %>% filter(sim=="model") %>% dplyr::select(1:5) %>% apply(c(2),FUN=row_number)/(nrow(sims)*(1-v)+1)) %>% apply(MARGIN=1,FUN = max)

datmax <- grid.arrange(PP_plot(observed = observed2,simulated = simulated2, title = TeX("$K=\\{1,2\\}$")),
                       PP_plot(observed = observed3,simulated = simulated3, title = TeX("$K=\\{1,2,3\\}$")),
                       PP_plot(observed = observed4,simulated = simulated4, title = TeX("$K=\\{1,2,3,4\\}$")),
                       PP_plot(observed = observed5,simulated = simulated5, title = TeX("$K=\\{1,2,3,4,5\\}$")),ncol=2)

ggsave(filename = paste0("plots/datmaxcond",l,"Nsim",N_sim,".png"),datmax,width=10,height=10)
}

p_lowest <- (p_est*1000 - 1.96 *(p_est*1000*(1-p_est*1000)/N_sim)^(1/2))
p_uppest <- (p_est*1000 + 1.96 *(p_est*1000*(1-p_est*1000)/N_sim)^(1/2))
p_est*1000
df <- data.frame("p_est"=p_est*1000,"CI"=paste0("(",p_lowest,", ",p_uppest,")"))
df
#/(1-N_sim)

# calculate probabilities
set.seed(11)
N <- 5000000
v <- 0.99
sims <- generate_Y(N=N) %>% link_log(dep=1/2) %>%
  link_log(dep=1/2) %>% link_log(dep=1/2) %>% link_log(dep=1/2) %>%
  apply(c(1,2),FUN=frechet_laplace_pit) %>% as.data.frame()

combn(x=c(1,2,3,4,5),m=2,simplify=TRUE)

# calculate true probability of trivariate extreme
# 0.00194
# simulate above threshold
y <- 0.995
Y2_gen <- -log(2*(1-y)) + rexp(1000)
# transform to frechet margins
df <- data.frame(Y2=Y2_gen)%>%
  apply(c(1,2),FUN=laplace_frechet_pit) %>% 
  as.data.frame() %>% link_log() %>% relocate(Y2) %>% 
  link_log() 
thres <- qfrechet(0.995)
p <- mean(df$Y1>thres & df$Y2>thres & df$Y3>thres)*(1-y)
0.00196-1/(2*(sqrt(1000)))

# calculate true probability of all extreme
v_sim <- 0.999
Y2_gen <- -log(2*(1-v_sim)) + rexp(500000)
df <- data.frame(Y2=Y2_gen)%>%
  apply(c(1,2),FUN=laplace_frechet_pit) %>% 
  as.data.frame() %>% link_log() %>% relocate(Y2) %>% 
  link_log() %>% link_log() %>% link_log()
vF_sim <- qfrechet(v_sim)
p_true <- df %>% filter(Y1>vF_sim,Y2>vF_sim,Y3>vF_sim,Y4>vF_sim,Y5>vF_sim) %>% 
  nrow()/500000*(1-v_sim)

# explore summer and winter air pollution data ----
as.data.frame(winter %>% apply(c(2),FUN=row_number)/(nrow(winter)+1)) %>% ggpairs()
as.data.frame(summer %>% apply(c(2),FUN=row_number)/(nrow(summer)+1)) %>% ggpairs()
# fit vine copula to all the data
fitw <- RVineStructureSelect(as.data.frame(winter %>% apply(c(2),FUN=row_number)/(nrow(winter)+1)),
                     trunclevel = 3, indeptest = FALSE)
# plot simulated and observed data on unifrom margins
winter_sim <- RVineSim(N=nrow(winter),RVM=fitw)
sims <- as.data.frame(winter %>% apply(c(2),FUN=row_number)/(nrow(winter)+1))
sl <- winter_sim
grid.arrange(PP_plot(observed=as.numeric(sims[,1]),simulated = as.numeric(sl[,1]),title = names(sims)[1],CIcol = "#009ADA"),
             PP_plot(observed=as.numeric(sims[,2]),simulated = as.numeric(sl[,2]),title = names(sims)[2],CIcol = "#009ADA"),
             PP_plot(observed=as.numeric(sims[,3]),simulated = as.numeric(sl[,3]),title = names(sims)[3],CIcol = "#009ADA"),
             PP_plot(observed=as.numeric(sims[,4]),simulated = as.numeric(sl[,4]),title = names(sims)[4],CIcol = "#009ADA"),
             PP_plot(observed=as.numeric(sims[,5]),simulated = as.numeric(sl[,5]),title = names(sims)[5],CIcol = "#009ADA"),ncol=2)

rbind(as.data.frame(winter %>% apply(c(2),FUN=row_number)/(nrow(winter)+1)) %>% mutate(res=rep("data",nrow(winter))),
      winter_sim %>% as.data.frame() %>%
              mutate(res=rep("model",nrow(winter)))) %>% 
  ggpairs(columns = 1:5,ggplot2::aes(color=res,alpha=0.5), upper = list(continuous = wrap("cor", size = 2.5))) +
  scale_color_manual(values = c("data"="black","model" = "#009ADA")) + scale_fill_manual(values = c("data"="black","model" = "#009ADA"))

fits <- RVineStructureSelect(as.data.frame(summer %>% apply(c(2),FUN=row_number)/(nrow(summer)+1)),
                     trunclevel = 3, indeptest = FALSE)
summer_sim <- RVineSim(N=nrow(summer),RVM=fits)

rbind(as.data.frame(summer %>% apply(c(2),FUN=row_number)/(nrow(summer)+1)) %>% mutate(res=rep("data",nrow(summer))),
      summer_sim %>% as.data.frame() %>%
        mutate(res=rep("model",nrow(summer)))) %>% 
  ggpairs(columns = 1:5,ggplot2::aes(color=res,alpha=0.5), upper = list(continuous = wrap("cor", size = 2.5))) +
  scale_color_manual(values = c("data"="black","model" = "#C11432")) + scale_fill_manual(values = c("data"="black","model" = "#C11432"))

sims <- as.data.frame(summer %>% apply(c(2),FUN=row_number)/(nrow(summer)+1))
sl <- summer_sim
grid.arrange(PP_plot(observed=as.numeric(sims[,1]),simulated = as.numeric(sl[,1]),title = names(sims)[1]),
             PP_plot(observed=as.numeric(sims[,2]),simulated = as.numeric(sl[,2]),title = names(sims)[2]),
             PP_plot(observed=as.numeric(sims[,3]),simulated = as.numeric(sl[,3]),title = names(sims)[3]),
             PP_plot(observed=as.numeric(sims[,4]),simulated = as.numeric(sl[,4]),title = names(sims)[4]),
             PP_plot(observed=as.numeric(sims[,5]),simulated = as.numeric(sl[,5]),title = names(sims)[5]),ncol=2)


# transform to Laplace margins
winter_lap <- as.data.frame((winter %>% apply(c(2),FUN=row_number))/(nrow(winter)+1)) %>% apply(c(1,2),FUN=unif_laplace_pit) %>% as.data.frame()
summer_lap <- as.data.frame((summer %>% apply(c(2),FUN=row_number))/(nrow(summer)+1)) %>% apply(c(1,2),FUN=unif_laplace_pit) %>% as.data.frame()

# fit vine copula to the observed residuals
colnames(winter_lap) <- paste0("Y",1:5)
wintercol <- "#009ADA"
obsz <- observed_residuals(df = winter_lap,v=0.7,given = 1)
v <- 0.7
fit <- RVineStructureSelect((observed_residuals(df = winter_lap,given = 1,v = v) %>% apply(c(2),FUN=row_number))/(nrow(winter_lap)*(1-v)+1),
                     trunclevel = 3, indeptest = FALSE)
N_sim <- nrow(winter)*(1-v)
Zsim <- RVineSim(N=N_sim,RVM=fit)
# for loop to make PP plot for each variable
p <- list()
for (i in sequence(ncol(obsz))) {
 observed <- as.data.frame(obsz %>% dplyr::select(all_of(i)) %>% apply(c(2),FUN=row_number)/(nrow(winter)*(1-v)+1)) %>% pull()
 simulated <- as.data.frame(Zsim) %>% dplyr::select(all_of(i)) %>% pull()
 p[[i]] <- PP_plot(observed = observed,simulated = simulated, title = TeX(paste0("Z_",i+1)), CIcol = wintercol)
 } 
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],ncol=2)  
# for loop to make PP plot for maxima over different subsets K
p <- list()
for (i in sequence(ncol(obsz)-1)) {
  observed <- as.data.frame(obsz %>% dplyr::select(all_of(1:(i+1))) %>% apply(c(2),FUN=row_number)/(nrow(winter)*(1-v)+1)) %>% apply(MARGIN=1,FUN = max)
  simulated <- as.data.frame(Zsim) %>% dplyr::select(all_of(1:(i+1))) %>% apply(MARGIN=1,FUN=max)
  Ktitle <- TeX(paste0("$K=\\{",mapply(function(x) paste(1:(x+1), collapse = ","),i),"\\}$"))
  p[[i]] <- PP_plot(observed = observed,simulated = simulated, title = Ktitle, CIcol = wintercol)
} 
grid.arrange(p[[1]],p[[2]],p[[3]],ncol=3) 

# observed residuals conditioning on each variable
for (j in 1:5) {
  obsz <- observed_residuals(df = winter_lap,v=0.7,given = j)
  ggsave(ggpairs(obsz),filename = paste0("plots/pollution_winter_obs_z",j,".png"))
  print(RVineStructureSelect((observed_residuals(df = winter_lap,given = j,v = v) %>% apply(c(2),FUN=row_number))/(nrow(winter_lap)*(1-v)+1),
                       trunclevel = 3, indeptest = FALSE))
}
# plot observed and simulated residuals for each variable
for (j in 1:5) {
  obs_res <- observed_residuals(df = summer_lap,v=0.7,given = j)
  fit3 <- RVineStructureSelect((observed_residuals(df = summer_lap,given = j,v = v) %>% apply(c(2),FUN=row_number))/(nrow(summer_lap)*(1-v)+1),
                       trunclevel = 3, indeptest = FALSE)
  print(fit3)
  N_sim <- nrow(obs_res)
  Zsim <- RVineSim(N=N_sim,RVM=fit3)
  # transform back residuals to original margins
  # can use kernel smoothed distribution as initial step
  to_opt <- function(z) {
    return( (mean(pnorm((z-obs_res[,k] %>% pull())/density(obs_res[,k] %>% pull())$bw)) - Zsim[i,k])^2)
  }
  Z <- Zsim
  for (i in 1:nrow(Zsim)) {
    for (k in 1:ncol(Zsim)) {
      Z[i,k] <- optim(fn=to_opt,par=1)$par
    }
  }
 p <-  rbind(obs_res %>% as.data.frame() %>% mutate(res=rep("data",N_sim)),
        Z %>% as.data.frame() %>%
          mutate(res=rep("model",N_sim))) %>% 
    ggpairs(columns = 1:4,ggplot2::aes(color=res,alpha=0.5), upper = list(continuous = wrap("cor", size = 2.5))) +
    scale_color_manual(values = c("data"="black","model" = "#C11432")) + scale_fill_manual(values = c("data"="black","model" = "#C11432"))
 ggsave(ggpairs(obs_res),filename = paste0("plots/pollution_summer_obs_z",j,".png"))
  ggsave(p,filename = paste0("plots/pollution_summer_obs_sim_z",j,".png"))
}
# plot observed and simulated data on Laplace margins cond. on each variable
for (l in 1:5) {
fit3 <- RVineStructureSelect((observed_residuals(df = winter_lap,given = l,v = v) %>% apply(c(2),FUN=row_number))/(nrow(winter_lap)*(1-v)+1),
                             trunclevel = 3, indeptest = FALSE)
N_sim <- ceiling(nrow(winter_lap)*(1-v))
Zsim <- RVineSim(N=N_sim,RVM=fit3)
# transform back residuals to original margins
# can use kernel smoothed distribution as initial step
obs_res <- observed_residuals(df = winter_lap,given = l,v = v) 
to_opt <- function(z) {
  return( (mean(pnorm((z-obs_res[,k] %>% pull())/density(obs_res[,k] %>% pull())$bw)) - Zsim[i,k])^2)
}
Z <- Zsim
for (i in 1:nrow(Zsim)) {
  for (k in 1:ncol(Zsim)) {
    Z[i,k] <- optim(fn=to_opt,par=1)$par
  }
}

Z_star <- as.data.frame(Z)
v_sim <- 0.7
X1_gpd <- ismev::gpd.fit(winter[,l], threshold = quantile(probs=v_sim,winter[,l]))
X1_gen <- rgpd(n=N_sim,sigma = X1_gpd$mle[1],xi=X1_gpd$mle[2],u=X1_gpd$threshold)
Y1_gen <- as.data.frame((data.frame(X1_gen) %>% apply(c(2),FUN=row_number))/(N_sim+1)*(1-v_sim)+v_sim) %>% apply(c(1,2),FUN=unif_laplace_pit) %>% as.data.frame() 
Gen_Y1 <- data.frame("Y1"=Y1_gen)
names(Gen_Y1) <- c("Y1")
# for each Y, generate a residual and calculate Y_2
a_hat <- par_est(df = winter_lap, v=v,given = l,margin = "Normal", method = "one_step")$a
b_hat <- par_est(df = winter_lap, v=v,given = l,margin = "Normal", method = "one_step")$b
res <- c(1:5)[-l]
Y1 <- Gen_Y1$Y1
Y2 <- a_hat[1]*Y1 + Y1^b_hat[1] *Z_star[,1]
Y3 <-  a_hat[2]*Y1 + Y1^b_hat[2] *Z_star[,2]
Y4 <- a_hat[3]*Y1 + Y1^b_hat[3] *Z_star[,3]
Y5 <-  a_hat[4]*Y1 + Y1^b_hat[4] *Z_star[,4]

Gen_Y1 <- Gen_Y1 %>% mutate(Y2=Y2,Y3=Y3,Y4=Y4,Y5=Y5) %>% mutate(sim=rep("model",nrow(Z_star)))
names(Gen_Y1) <- c(paste0("Y",l),paste0("Y",res[1]),paste0("Y",res[2]),paste0("Y",res[3]),paste0("Y",res[4]),"sim")
# generate Y1 (extrapolate so above largest observed value)

#plot
Y_given_1_extreme <- winter_lap %>% filter(winter_lap[,l]>quantile(winter_lap[,l],v))
Y1 <- Y_given_1_extreme[,l]
Y2 <- Y_given_1_extreme[,res[1]]
Y3 <- Y_given_1_extreme[,res[2]]
Y4 <- Y_given_1_extreme[,res[3]]
Y5 <- Y_given_1_extreme[,res[4]]
tmp <- data.frame(Y1,Y2,Y3,Y4,Y5) %>% mutate(sim=c("data"))
names(tmp) <- c(paste0("Y",l),paste0("Y",res[1]),paste0("Y",res[2]),paste0("Y",res[3]),paste0("Y",res[4]),"sim")
# l <- min(Gen_Y_1 %>% dplyr::select(-sim),Y1,Y2,Y3,Y4,Y5)
# u <- max(Gen_Y_1 %>% dplyr::select(-sim),Y1,Y2,Y3,Y4,Y5)
Gen_orig <- rbind(Gen_Y1,tmp)
Gen_orig <- Gen_orig %>% relocate(1,.before = l+1) # relocate cond variable from first to original position
names(Gen_orig) <- c(names(winter_lap),"sim")
p <- ggpairs(Gen_orig,columns = 1:5,ggplot2::aes(color=sim,alpha=0.5), upper = list(continuous = wrap("cor", size = 2.5))) +
  scale_color_manual(values = c("data"="black","model" = "#009ADA")) + scale_fill_manual(values = c("data"="black","model" = "#009ADA"))
ggsave(p,filename = paste0("plots/pollution_winter_sim",l,".png"))
}
v <- 0.7
tbl_winter <- par_est(df = winter_lap, v=v,given = 1:5,margin = "Normal", method = "one_step") %>% 
  dplyr::select(lik,a,b,mu,sig,given,res) %>% 
  mutate(cond_pollutant = recode(given, `1` = names(winter_lap)[1], `2` = names(winter_lap)[2], `3` = names(winter_lap)[3],
                                 `4` = names(winter_lap)[4], `5` = names(winter_lap)[5])) %>% 
  mutate(res_pollutant = recode(res, `1` = names(winter_lap)[1], `2` = names(winter_lap)[2], `3` = names(winter_lap)[3],
                                 `4` = names(winter_lap)[4], `5` = names(winter_lap)[5])) %>% 
  relocate(6,.after=9) %>% relocate(6,.after = 9) %>% dplyr::select(c(1:7))
xtable(tbl_winter[c(1:4,5,9,13,17),])
xtable(tbl_winter[c(5:8,1,10,14,18),])
xtable(tbl_winter[c(9:12,2,6,15,19),])
xtable(tbl_winter[c(13:16,3,7,11,20),])
xtable(tbl_winter[c(17:20,4,8,12,16),])
 
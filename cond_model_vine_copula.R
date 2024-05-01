# libraries and functions ----
library(VineCopula)
library(evd)
library(tidyverse)
library(latex2exp)
library(gridExtra)
library(viridis)
library(MASS) #use dplyr::select to avoid function conflict
source("cond_model_helpers.R")

# set theme defaults to be black rectangle
theme_set(theme_bw())
theme_replace(
  panel.spacing = unit(2, "lines"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA) )

# generate from the model
set.seed(1)
N <- 50000
v <- 0.99
sims_tmp <- generate_Y(N=N) %>% link_log(dep=1/2) %>% link_log(dep=1/2) %>% link_log(dep=1/2) %>% link_log(dep=1/2)
sims <- sims_tmp %>% apply(c(1,2),FUN=frechet_laplace_pit) %>% as.data.frame()

# calculate the observed residuals
df <- sims %>% dplyr::select(starts_with("Y"))
j <- 1
a_hat <- c()
b_hat <- c()
res_var <- c()
tmp_z <- c()
d <- ncol(df)
  Y_given_1_extreme <- df %>% filter(df[,j]>quantile(df[,j],v))
  n_v <- nrow(Y_given_1_extreme)
  res <- c(1:d)[-j]
  init_par <- c()
  for (i in 2:d) {
    # optimise using the initial parameters
    init_opt <- optim(par=c(0.5,0,1), fn=Y_likelihood_initial,df=Y_given_1_extreme,given=j,sim=res[i-1],control = list(fnscale=-1))
    init_par <- c(init_opt$par[1],0.2,init_opt$par[2],init_opt$par[3])
    opt <- optim(par=init_par,fn = Y_likelihood,df=Y_given_1_extreme,given=j,sim=res[i-1],control = list(fnscale=-1))
   # a_hat <- opt$par[1]
    a_hat <- 1
  #  b_hat <- opt$par[2]
    b_hat <- 0
    res_var <- append(res_var,rep(paste0("Z",res[i-1]),n_v))
    Y1 <- Y_given_1_extreme[,j]
    Y2 <- Y_given_1_extreme[,res[i-1]]
    tmp_z <- append(tmp_z,(Y2-a_hat*Y1/(Y1^b_hat)))
  }
  
obs_res <- data.frame(res_var,tmp_z) %>% mutate(res_var=factor(res_var,levels=paste0("Z",res))) %>% group_by(res_var) %>% 
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = res_var, values_from = tmp_z) %>% 
  dplyr::select(-row)
pairs(obs_res)

# transform to Uniform margins
obs_res1 <- (obs_res %>% apply(c(2),FUN=row_number))/(N+1)
#sims1 <- (sims %>% apply(c(2),FUN=row_number))/(N+1) #transform data to uniform margins
# model using VineCopula package
fit <- RVineStructureSelect(obs_res1)
fit
plot(fit,edge.labels = "family")
<- # simulate from the copula

# transform back residuals to original margins

# simulate

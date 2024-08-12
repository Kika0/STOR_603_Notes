# libraries and functions ----
library(VineCopula)
library(evd)
library(ismev)
library(tidyverse)
library(latex2exp)
library(gridExtra)
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
# calculate the observed residuals
observed_residuals <- function(df=sims,given=1,v=0.99) {
  j <- given
  a_hat <- b_hat <- res_var <- c()
  tmp_z <- tmp_z1 <- c()
  df_orig <- df
  names(df) <- paste0("Y",1:ncol(df))
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
    a_hat <- opt$par[1]
    # a_hat <- 1
    b_hat <- opt$par[2]
    # b_hat <- 0
    res_var <- append(res_var,rep(paste0("Z",res[i-1]),n_v))
    Y1 <- Y_given_1_extreme[,j]
    Y2 <- Y_given_1_extreme[,res[i-1]]
    tmp_z <- append(tmp_z,(Y2-a_hat*Y1/(Y1^b_hat)))
  }
  Z <- data.frame(res_var,tmp_z) %>% mutate(res_var=factor(res_var,levels=paste0("Z",res))) %>% group_by(res_var) %>% 
    mutate(row = row_number()) %>%
    tidyr::pivot_wider(names_from = res_var, values_from = tmp_z) %>% 
    dplyr::select(-row)
  if (sum(names(df_orig)!=paste0("Y",1:ncol(df_orig)))==ncol(df_orig)) {
    names(Z) <- names(df_orig)[-j]
  }
  return(Z)
}

# generate from the model ----
set.seed(11)
N <- 5000
v <- 0.99
sims <- generate_Y(N=N) %>% link_log(dep=1/2) %>%
  link_log(dep=1/2) %>% link_log(dep=1/2) %>% link_log(dep=1/2) %>%
  apply(c(1,2),FUN=frechet_laplace_pit) %>% as.data.frame()
# explore residuals transformed to uniform margins
# ggpairs((observed_residuals(df = sims,given = 1,v = 0.99) %>% apply(c(2),FUN=row_number))/(n_v+1))
# transform to Uniform margins and fit a vine
fit3 <- RVineStructureSelect((observed_residuals(df = sims,given = 1,v = v) %>% apply(c(2),FUN=row_number))/(nrow(sims)*(1-v)+1),
                        trunclevel = 3, indeptest = TRUE)
fit3
fit2 <- RVineStructureSelect((observed_residuals(df = sims,given = 1,v = v) %>% apply(c(2),FUN=row_number))/(nrow(sims)*(1-v)+1),
                             trunclevel = 2, indeptest = FALSE)
fit2
fit1 <- RVineStructureSelect((observed_residuals(df = sims,given = 1,v = v) %>% apply(c(2),FUN=row_number))/(nrow(sims)*(1-v)+1),
                             trunclevel = 1, indeptest = FALSE)
fit1
RVineClarkeTest(data=(observed_residuals(df = sims,given = 2,v = 0.99) %>% apply(c(2),FUN=row_number))/(nrow(sims)*(1-v)+1),
                RVM1 = fit1,RVM2 = fit2)

plot(fit3,edge.labels = "family")
# simulate from the copula
N_sim <- 50
Zsim <- RVineSim(N=N_sim,RVM=fit3)
# transform back residuals to original margins
# can use kernel smoothed distribution as initial step
obs_res <- observed_residuals(df = sims,given = 1,v = 0.99) 
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
  for (j in 1:ncol(obsr)) {
    obsr[i,j] <- obsr[i,j]/(nrow(sims)*(1-v)+1)
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
Z_star <- as.data.frame(Z)
U <- runif(500)
Y_1_gen <- -log(2*(1-0.999)) + rexp(500)
Gen_Y_1 <- data.frame(Y1=Y_1_gen)

# for each Y, generate a residual and calculate Y_2
Y1 <- Gen_Y_1$Y1
Y2 <- a_hat*Y1 + Y1^b_hat *Z_star[,1]
Y3 <-  a_hat*Y1 + Y1^b_hat *Z_star[,2]
Y4 <- a_hat*Y1 + Y1^b_hat *Z_star[,3]
Y5 <-  a_hat*Y1 + Y1^b_hat *Z_star[,4]
res <- c(1:d)[-j]
Gen_Y_1 <- Gen_Y_1 %>% mutate(Y2=Y2,Y3=Y3,Y4=Y4,Y5=Y5) %>% mutate(sim=rep("model",nrow(Z_star)))
names(Gen_Y_1) <- c(paste0("Y",j),paste0("Y",res[1]),paste0("Y",res[2]),paste0("Y",res[3]),paste0("Y",res[4]),"sim")
# generate Y_1 (extrapolate so above largest observed value)

#plot
d <- ncol(df)
Y_given_1_extreme <- sims %>% filter(sims[,j]>quantile(sims[,j],v))

tmp <- data.frame(Y_given_1_extreme[,j],
                Y_given_1_extreme[,res[1]],
                Y_given_1_extreme[,res[2]],
                Y_given_1_extreme[,res[3]],
                Y_given_1_extreme[,res[4]]) %>%
  mutate(sim=rep("data",nrow(Y_given_1_extreme)))
names(tmp) <- c(paste0("Y",j),paste0("Y",res[1]),paste0("Y",res[2]),paste0("Y",res[3]),paste0("Y",res[4]),"sim")
Gen_orig <- rbind(Gen_Y_1,tmp)
ggpairs(Gen_orig,columns = 1:5,ggplot2::aes(color=sim,alpha=0.5), upper = list(continuous = wrap("cor", size = 2.5))) +
  scale_color_manual(values = c("data"="black","model" = "#C11432")) + scale_fill_manual(values = c("data"="black","model" = "#C11432"))

# for loop to condition on each variable ----
v <- 0.99 # threshold for conditioning
v_sim <- 0.999 # threshold for simulation
N_sim <- 100 # number of observations to simulate
d <- 5
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
tmp <- data.frame(Y1,Y2,Y3,Y4,Y5) %>% mutate(sim=rep("data",nrow(sims)*(1-v)))
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
RVineStructureSelect(as.data.frame(winter %>% apply(c(2),FUN=row_number)/(nrow(winter)+1)),
                     trunclevel = 3, indeptest = FALSE)
RVineStructureSelect(as.data.frame(summer %>% apply(c(2),FUN=row_number)/(nrow(summer)+1)),
                     trunclevel = 3, indeptest = FALSE)

# transform to Laplace margins
winter_lap <- as.data.frame((winter %>% apply(c(2),FUN=row_number))/(nrow(winter)+1)) %>% apply(c(1,2),FUN=unif_laplace_pit) %>% as.data.frame()
summer_lap <- as.data.frame((summer %>% apply(c(2),FUN=row_number))/(nrow(summer)+1)) %>% apply(c(1,2),FUN=unif_laplace_pit) %>% as.data.frame()

# fit vine copula to the observed residuals
colnames(winter_lap) <- paste0("Y",1:5)
obsz <- observed_residuals(df = winter_lap,v=0.7,given = 1)
v <- 0.7
RVineStructureSelect((observed_residuals(df = winter_lap,given = 1,v = v) %>% apply(c(2),FUN=row_number))/(nrow(winter_lap)*(1-v)+1),
                     trunclevel = 3, indeptest = FALSE)
for (j in 1:5) {
  obsz <- observed_residuals(df = winter_lap,v=0.7,given = j)
  ggsave(ggpairs(obsz),filename = paste0("plots/pollution_winter_obs_z",j,".png"))
  print(RVineStructureSelect((observed_residuals(df = winter_lap,given = j,v = v) %>% apply(c(2),FUN=row_number))/(nrow(winter_lap)*(1-v)+1),
                       trunclevel = 3, indeptest = FALSE))
}

for (j in 1:5) {
  obsz <- observed_residuals(df = winter_lap,v=0.7,given = j)
  ggsave(ggpairs(obsz),filename = paste0("plots/pollution_winter_obs_z",j,".png"))
  fit3 <- RVineStructureSelect((observed_residuals(df = winter_lap,given = j,v = v) %>% apply(c(2),FUN=row_number))/(nrow(winter_lap)*(1-v)+1),
                       trunclevel = 3, indeptest = FALSE)
  print(fit3)
  N_sim <- nrow(obsz)
  Zsim <- RVineSim(N=N_sim,RVM=fit3)
  # transform back residuals to original margins
  # can use kernel smoothed distribution as initial step
  obs_res <- observed_residuals(df = sims,given = j,v = 0.99) 
  to_opt <- function(z) {
    return( (mean(pnorm((z-obs_res[,k] %>% pull())/density(obs_res[,k] %>% pull())$bw)) - Zsim[i,k])^2)
  }
  Z <- Zsim
  for (i in 1:nrow(Zsim)) {
    for (k in 1:ncol(Zsim)) {
      Z[i,k] <- optim(fn=to_opt,par=1)$par
    }
  }
 p <-  rbind(obsz %>% as.data.frame() %>% mutate(res=rep("data",N_sim)),
        Z %>% as.data.frame() %>%
          mutate(res=rep("model",N_sim))) %>% 
    ggpairs(columns = 1:4,ggplot2::aes(color=res,alpha=0.5), upper = list(continuous = wrap("cor", size = 2.5))) +
    scale_color_manual(values = c("data"="black","model" = "#009ADA")) + scale_fill_manual(values = c("data"="black","model" = "#009ADA"))
  ggsave(p,filename = paste0("plots/pollution_winter_obs_sim_z",j,".png"))
}

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

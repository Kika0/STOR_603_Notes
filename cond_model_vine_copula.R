# libraries and functions ----
library(VineCopula)
library(evd)
library(tidyverse)
library(latex2exp)
library(gridExtra)
library(GGally) # for ggpairs function
library(MASS) # use dplyr::select to avoid function conflict
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
    a_hat1 <- opt$par[1]
    a_hat <- 1
    b_hat1 <- opt$par[2]
    b_hat <- 0
    res_var <- append(res_var,rep(paste0("Z",res[i-1]),n_v))
    Y1 <- Y_given_1_extreme[,j]
    Y2 <- Y_given_1_extreme[,res[i-1]]
    tmp_z <- append(tmp_z,(Y2-a_hat*Y1/(Y1^b_hat)))
    tmp_z1 <- append(tmp_z1,(Y2-a_hat1*Y1/(Y1^b_hat1)))
  }
  Z <- data.frame(res_var,tmp_z) %>% mutate(res_var=factor(res_var,levels=paste0("Z",res))) %>% group_by(res_var) %>% 
    mutate(row = row_number()) %>%
    tidyr::pivot_wider(names_from = res_var, values_from = tmp_z) %>% 
    dplyr::select(-row)
  return(Z)
}

# generate from the model
set.seed(11)
N <- 50000
v <- 0.99
sims <- generate_Y(N=N) %>% link_log(dep=1/2) %>%
  link_log(dep=1/2) %>% link_log(dep=1/2) %>% link_log(dep=1/2) %>%
  apply(c(1,2),FUN=frechet_laplace_pit) %>% as.data.frame()
# explore residuals transformed to uniform margins
# ggpairs((observed_residuals(df = sims,given = 1,v = 0.99) %>% apply(c(2),FUN=row_number))/(n_v+1))
# transform to Uniform margins and fit a vine
fit3 <- RVineStructureSelect((observed_residuals(df = sims,given = 3,v = 0.99) %>% apply(c(2),FUN=row_number))/(nrow(sims)*(1-v)+1),
                        trunclevel = 3, indeptest = FALSE)
fit3
fit2 <- RVineStructureSelect((observed_residuals(df = sims,given = 1,v = 0.99) %>% apply(c(2),FUN=row_number))/(nrow(sims)*(1-v)+1),
                             trunclevel = 2, indeptest = FALSE)
fit2
fit1 <- RVineStructureSelect((observed_residuals(df = sims,given = 1,v = 0.99) %>% apply(c(2),FUN=row_number))/(nrow(sims)*(1-v)+1),
                             trunclevel = 1, indeptest = FALSE)
fit1
RVineClarkeTest(data=(observed_residuals(df = sims,given = 1,v = 0.99) %>% apply(c(2),FUN=row_number))/(nrow(sims)*(1-v)+1),
                RVM1 = fit1,RVM2 = fit2)

plot(fit3,edge.labels = "family")
# simulate from the copula
Zsim <- RVineSim(N=500,RVM=fit3)
# transform back residuals to original margins
# can use kernel smoothed distribution as initial step
obs_res <- observed_residuals(df = sims,given = 3,v = 0.99) 
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

rbind(obs_res %>% as.data.frame() %>% mutate(res=rep("data",500)),
      Z %>% as.data.frame() %>%
        mutate(res=rep("model",500))) %>% 
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

# for loop to condition on each variable
for (l in 1:5) {
fit3 <- RVineStructureSelect((observed_residuals(df = sims,given = l,v = 0.99) %>% apply(c(2),FUN=row_number))/(nrow(sims)*(1-v)+1),
                             trunclevel = 3, indeptest = FALSE)
Zsim <- RVineSim(N=500,RVM=fit3)
d <- 5
res <- c(1:d)[-l]
# transform back residuals to original margins
# can use kernel smoothed distribution as initial step
obs_res <- observed_residuals(df = sims,given = l,v = 0.99) 
to_opt <- function(z) {
  return( (mean(pnorm((z-obs_res[,k] %>% pull())/density(obs_res[,k] %>% pull())$bw)) - Zsim[i,k])^2)
}
Z <- Zsim
for (i in 1:nrow(Zsim)) {
  for (k in 1:ncol(Zsim)) {
    Z[i,k] <- optim(fn=to_opt,par=1)$par
  }
}
# simulate
Z_star <- as.data.frame(Z)
Y_1_gen <- -log(2*(1-0.99)) + rexp(500)
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
tmp <- data.frame(Y1,Y2,Y3,Y4,Y5) %>% mutate(sim=rep("data",500))
names(tmp) <- c(paste0("Y",l),paste0("Y",res[1]),paste0("Y",res[2]),paste0("Y",res[3]),paste0("Y",res[4]),"sim")
thres <- frechet_laplace_pit( qfrechet(0.99))
v <- 0.99
# l <- min(Gen_Y_1 %>% dplyr::select(-sim),Y1,Y2,Y3,Y4,Y5)
# u <- max(Gen_Y_1 %>% dplyr::select(-sim),Y1,Y2,Y3,Y4,Y5)
Gen_orig <- rbind(Gen_Y_1,tmp)
p <- ggpairs(Gen_orig,columns = 1:5,ggplot2::aes(color=sim,alpha=0.5), upper = list(continuous = wrap("cor", size = 2.5))) +
  scale_color_manual(values = c("data"="black","model" = "#C11432")) + scale_fill_manual(values = c("data"="black","model" = "#C11432"))
ggsave(p,filename=(paste0("cond",l,".pdf")),device="pdf")
}

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

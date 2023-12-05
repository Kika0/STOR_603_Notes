# libraries and functions ----
library(evd)
library(tidyverse)
library(latex2exp)
library(gridExtra)
source("cond_model_helpers.R")

# set theme defaults to be black rectangle
theme_set(theme_bw())
theme_replace(
  panel.spacing = unit(2, "lines"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA) )


# generate trivariate sample ----
N <- 5000
sims <- generate_dep_X_Y_Y_Z(N=N)

# PIT to Laplace
sims <- sims %>% mutate(Y_1=as.numeric(map(.x=X_1,.f=frechet_laplace_pit))) %>% 
  mutate(Y_2=as.numeric(map(.x=X_2,.f=frechet_laplace_pit))) %>%
  mutate(Y_3=as.numeric(map(.x=X_3,.f=frechet_laplace_pit)))

ggplot(sims %>% select(Y_1,Y_2,Y_3) %>% pivot_longer(everything())) + geom_density(aes(x=value),,stat="density") + facet_wrap(~name)

grid.arrange(ggplot(sims) + geom_point(aes(x=Y_1,y=Y_2),alpha=0.5),
             ggplot(sims) + geom_point(aes(x=Y_2,y=Y_3),alpha=0.5),
             ggplot(sims) + geom_point(aes(x=Y_1,y=Y_3),alpha=0.5),ncol=3)

# filter for Y_1 being extreme -----
v <- 0.99
Y_given_1_extreme <- sims %>% filter(Y_1>quantile(Y_1,v))
Y_not_1_extreme <- sims %>% filter(Y_1<quantile(Y_1,v))


opt <- optim(par=c(1,0,0,1),fn = Y_2_likelihood,df=Y_given_1_extreme,given=1,sim=2,control = list(fnscale=-1))
a_hat <- opt$par[1]
b_hat <- opt$par[2]
# extrapolate using kernel smoothed residuals
N <- 5000
Y_1 <- Y_given_1_extreme[,4]
Y_2 <- Y_given_1_extreme[,5]
Z <- (Y_2-a_hat*Y_1)/(Y_1^b_hat)


ggplot(Y_given_1_extreme %>% select(Y_1,Y_2_sim,Y_2,Y_3) %>% pivot_longer(everything())) + geom_density(aes(x=value),,stat="density") + facet_wrap(~name)

Y_2_sim <- bind_rows(Y_not_1_extreme %>% mutate(Y_2_sim=Y_2),Y_given_1_extreme) %>% 
  mutate(X_2_sim=rep(0,N))

X_2_simdf <- Y_2_sim %>% 
  mutate(X_2=as.numeric(map(.x=Y_2,.f=laplace_frechet_pit))) %>%
  mutate(X_2_sim=as.numeric(map(.x=Y_2_sim,.f=laplace_frechet_pit)))

X_2_simdf <- X_2_simdf %>% mutate(v=c(rep("below_threshold",N*v),rep("above_threshold",N*(1-v))))

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
v <- 0.99
Y_given_1_extreme <- sims %>% filter(Y_1>quantile(Y_1,v))
Y_not_1_extreme <- sims %>% filter(Y_1<quantile(Y_1,v))

opt <- optim(par=c(1,0,0,1),fn = Y_2_likelihood,df=Y_given_1_extreme,given=1,sim=3,control = list(fnscale=-1))
a_hat <- opt$par[1]
b_hat <- opt$par[2]

# plot the values inferenced on ----


Y_given_1_extreme <- Y_given_1_extreme %>% mutate(Y_3_sim=Y_3_sim)
ggplot(Y_given_1_extreme %>% select(Y_1,Y_3_sim,Y_2,Y_3) %>% pivot_longer(everything())) + geom_density(aes(x=value),,stat="density") + facet_wrap(~name)

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


# generate residual Z ----
Y_1 <- Y_given_1_extreme[,4]
Y_2 <- Y_given_1_extreme[,5]
Z <- (Y_2-a_hat*Y_1)/(Y_1^b_hat)
plot(Y_1,Z)

# generate X_1 from Frechet distribution above 0.9 quantile
U <- runif(50000)
X_1_gen <- sort( -1/(log(0.99) )  -1/(log(U) ) )
set.seed(12)
N <- 50000
U <- runif(min=0.99,max=1,N)
X_1_gen <- sort( -1/(log(U) ) )

# transform to Laplace margins
Gen_Y_1 <- data.frame(X_1=X_1_gen) %>% mutate(Y_1=as.numeric(map(.x=X_1,.f=frechet_laplace_pit)))

# for each Y, generate a residual and calculate Y_2
Y_1 <- Gen_Y_1$Y_1
Z_gen <- sample(Z,N,replace=TRUE) +rnorm(N,mean=0,sd=density(Z)$bw) # plus noise
Y_2 <- a_hat*Y_1 + Y_1^b_hat *Z_gen
Gen_Y_1 <- Gen_Y_1 %>% mutate(Y_2=Y_2) %>% mutate(sim=rep("conditional_model",N))
# generate Y_1 (extrapolate so above largest observed value)

#plot
Gen_orig <- rbind(Gen_Y_1,Y_given_1_extreme %>% select(X_1,Y_1,Y_2) %>% mutate(sim=rep("original_laplace",50)))
ggplot(Gen_orig) + geom_point(aes(x=Y_1,y=Y_2,col=sim),alpha=0.5) + scale_color_manual(values = c("original_laplace"="black",
                                                                                                 "conditional_model" = "#C11432")) 
# specify threshold for Laplace margin
v_l <- c(5,12,5,12)
# transform to logistic margins
v <- log(sapply(X = v_l,FUN = laplace_frechet_pit))
laplace_frechet_pit(5)

# calculate empirical probability by simulating Y_2 from the model
((Gen_Y_1 %>% filter(Y_1>v[1],Y_1<v[2],Y_2>v[3],Y_2<v[4]) %>% dim())[1]/50000)*(1-0.99)

# need to log Fréchet margins to get the logistic model margin
evd::pbvevd(c(v[2],v[4]),dep=dep[1],model="log") -
  evd::pbvevd(c(v[1],v[4]),dep=dep[1],model="log") -
  evd::pbvevd(c(v[2],v[3]),dep=dep[1],model="log") +
  evd::pbvevd(c(v[1],v[3]),dep=dep[1],model="log")

# suppose we wish to simulate 1/10000 year event probability
p10_4 <- frechet_laplace_pit(-1/(log(0.9999)))
v_l <- c(p10_4,100,p10_4,100)
# because of dependence, the probability of two variables being large together is larger than p^2


# quite good estimation of probability but perhaps generate more large samples to verify
simulation_prob <- function(Z=Z,a_hat=a_hat,b_hat=b_hat,v_l=c(5,12,5,12))


# different methods of generating from Fréchet lead to different results, which should be looked at

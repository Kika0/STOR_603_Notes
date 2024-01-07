# libraries and functions ----
library(evd)
library(tidyverse)
library(latex2exp)
library(gridExtra)
library(MASS) #use dplyr::select to avoid function conflict
library(xtable)
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

ggplot(sims %>% select(Y_1,Y_2,Y_3) %>% pivot_longer(everything())) + geom_density(aes(x=value),stat="density") + facet_wrap(~name)

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
# extrapolate using kernel smoothed residuals ----
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
 # U <- runif(50000)
 # X_1_gen <- sort( -1/log(0.99) + 1 - exp(-U) )
set.seed(12)
N <- 50000
U <- runif(min=0.99,max=1,N)
X_1_gen <- sort( -1/(log(U) ) )

U <- runif(50000)
Y_1_gen <- -log(2*(1-0.99)) + rexp(50000)
Gen_Y_1 <- data.frame(Y_1=Y_1_gen,X_1=as.numeric(map(.x=X_1,.f=laplace_frechet_pit)))

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
ggplot(Gen_orig) + geom_point(aes(x=Y_1,y=Y_2,col=sim),alpha=0.5) + 
  scale_color_manual(values = c("original_laplace"="black","conditional_model" = "#C11432")) 

# specify threshold for Laplace margin
v_l <- c(5,12,5,12)
# transform to logistic margins
v <- log(sapply(X = v_l,FUN = laplace_frechet_pit))

# calculate empirical probability by simulating Y_2 from the model
((Gen_Y_1 %>% filter(Y_1>v[1],Y_1<v[2],Y_2>v[3],Y_2<v[4]) %>% dim())[1]/50000)*(1-0.99)

# need to log Fréchet margins to get the logistic model margin
p <- evd::pbvevd(c(v[2],v[4]),dep=dep[1],model="log") -
  evd::pbvevd(c(v[1],v[4]),dep=dep[1],model="log") -
  evd::pbvevd(c(v[2],v[3]),dep=dep[1],model="log") +
  evd::pbvevd(c(v[1],v[3]),dep=dep[1],model="log")

# calculate CI
CI <- c(p-(1.96*(p*(1-p)/50000)^(0.5)),p+(1.96*(p*(1-p)/50000)^(0.5)))

ran_bern <- rbinom(n=1000,size = 50000,p=p)/50000
# density(ran_bern) %>% plot()
ggplot(data.frame(x=ran_bern)) + geom_density(aes(x=x),stat="density") 


# suppose we wish to simulate 1/10000 year event probability
p10_4 <- frechet_laplace_pit(-1/(log(0.9999)))
v_l <- c(p10_4,100,p10_4,100)
# because of dependence, the probability of two variables being large together is larger than p^2


d <- data.frame(x=ran_bern)
# quite good estimation of probability but perhaps generate more large samples to verify
p1 <- ggplot(data = d) + theme_bw() + 
  geom_density(aes(x=x, y = ..density..), color = 'black')
x <- ran_bern
# new code is below
q25 <- quantile(x,.025)
q975 <- quantile(x,.975)
medx <- median(x)
x.dens <- density(x)
df.dens <- data.frame(x = x.dens$x, y = x.dens$y)
p1 + geom_area(data = subset(df.dens, x >= q25 & x <= q975), 
              aes(x=x,y=y), fill = 'lightblue') +
  geom_vline(xintercept = p) 
 # geom_vline(xintercept = medx)
library(HDInterval)
library(ggridges)
ggplot(d, aes(x = x, y = 0, fill = stat(quantile))) + 
  geom_density_ridges_gradient(quantile_lines = TRUE, quantile_fun = hdi, vline_linetype = 2) +
  scale_fill_manual(values = c("transparent", "lightblue", "transparent"), guide = "none")

# different methods of generating from Fréchet lead to different results, which should be looked at

# practice transforming from uniform margins to normal and back
U <- runif(10000)
plot(density(U))
N <- qnorm(U)
plot(density(N))
plot(density(rnorm(10000)))

# transform back to uniform
pnorm(1.96)

# start from the beginning
# generate trivariate sample ----
N <- 5000
sims <- generate_dep_X_Y_Y_Z(N=N)

# PIT to Laplace
sims <- sims %>% mutate(Y_1=as.numeric(map(.x=X_3,.f=frechet_laplace_pit))) %>% 
  mutate(Y_2=as.numeric(map(.x=X_1,.f=frechet_laplace_pit))) %>%
  mutate(Y_3=as.numeric(map(.x=X_2,.f=frechet_laplace_pit)))

ggplot(sims %>% dplyr::select(Y_1,Y_2,Y_3) %>% pivot_longer(everything())) + geom_density(aes(x=value),stat="density") + facet_wrap(~name)

grid.arrange(ggplot(sims) + geom_point(aes(x=Y_1,y=Y_2),alpha=0.5),
             ggplot(sims) + geom_point(aes(x=Y_2,y=Y_3),alpha=0.5),
             ggplot(sims) + geom_point(aes(x=Y_1,y=Y_3),alpha=0.5),ncol=3)

# filter for Y_1 being extreme -----
v <- 0.99
Y_given_1_extreme <- sims %>% filter(Y_1>quantile(Y_1,v))
Y_not_1_extreme <- sims %>% filter(Y_1<quantile(Y_1,v))


d <- 3
for (i in c(2,3)) {
opt[[i-1]] <- optim(par=c(1,0,0,1),fn = Y_2_likelihood,df=Y_given_1_extreme,given=1,sim=i,control = list(fnscale=-1))
}
a_hat <- c(opt[[1]]$par[1],opt[[2]]$par[1])
b_hat <- c(opt[[1]]$par[2],opt[[2]]$par[2])

# generate residual Z ----
Y_1 <- Y_given_1_extreme[,4]
Y_2 <- Y_given_1_extreme[,5]
Y_3 <- Y_given_1_extreme[,6]

Z_2 <- (Y_2-a_hat[1]*Y_1)/(Y_1^b_hat[1])
Z_3 <- (Y_3-a_hat[2]*Y_1)/(Y_1^b_hat[2])
plot(Y_1,Z_2)
plot(Y_1,Z_3)



# calculate the normal using the PIT
Z_N_2 <- qnorm(F_smooth_Z(Z_2))
Z_N_3 <- qnorm(F_smooth_Z(Z_3))

rho_hat <- cor(Z_N_2,Z_N_3)

Z_N <- mvrnorm(n=1000,mu=c(0,0),Sigma=matrix(c(1,rho_hat,rho_hat,1),2,2))
Z <- data.frame(Z_2,Z_3)

# transform back to original margins
Z_star <- norm_to_orig(Z_N=Z_N,emp_res = Z)

U <- runif(1000)
Y_1_gen <- -log(2*(1-0.999)) + rexp(1000)
Gen_Y_1 <- data.frame(Y_1=Y_1_gen,X_1=as.numeric(map(.x=X_1,.f=laplace_frechet_pit)))

# for each Y, generate a residual and calculate Y_2
Y_1 <- Gen_Y_1$Y_1
# Y_2 <- a_hat*Y_1 + Y_1^b_hat *x
# Y_3 <-  a_hat*Y_1 + Y_1^b_hat *y
Y_2 <- a_hat*Y_1 + Y_1^b_hat *Z_star[,1]
Y_3 <-  a_hat*Y_1 + Y_1^b_hat *Z_star[,2]
Gen_Y_1 <- Gen_Y_1 %>% mutate(Y_2=Y_2,Y_3=Y_3) %>% mutate(sim=rep("conditional_model",100))
# generate Y_1 (extrapolate so above largest observed value)

#plot
Gen_orig <- rbind(Gen_Y_1,Y_given_1_extreme %>% dplyr::select(X_1,Y_1,Y_2,Y_3) %>% mutate(sim=rep("original_laplace",50)))
p1 <- ggplot(Gen_orig) + geom_point(aes(x=Y_1,y=Y_2,col=sim),alpha=0.5) + 
  scale_color_manual(values = c("original_laplace"="black","conditional_model" = "#C11432")) 
p2 <- ggplot(Gen_orig) + geom_point(aes(x=Y_2,y=Y_3,col=sim),alpha=0.5) + 
  scale_color_manual(values = c("original_laplace"="black","conditional_model" = "#C11432")) 
p3 <- ggplot(Gen_orig) + geom_point(aes(x=Y_1,y=Y_3,col=sim),alpha=0.5) + 
  scale_color_manual(values = c("original_laplace"="black","conditional_model" = "#C11432")) 
grid.arrange(p1,p2,p3,ncol=3)

Z_comp <- Z_star %>% mutate(Compare=rep("Optimise_all",100))
Z_comp <- rbind(Z_comp,Z_star%>% mutate(Compare=rep("Linear_segments",100)))
ggplot() +
  geom_point(Z_comp,mapping = aes(x=X1,y=X2)) +
  facet_wrap(~Compare) +
  xlab(TeX("$Z^*_{2|1}$")) +
  ylab(TeX("$Z^*_{3|1}$"))

# show that linear approximation is reasonable to only optimise for values near 0 or 1
u <- seq(0.02,0.98,length.out=49)
u1 <- seq(0.0001,0.9999,length.out=998)
ggplot() + 
  geom_line(data.frame(x=qnorm(u1),u=u1),mapping=aes(x=x,y=u),alpha=0.5) +
  geom_point(data.frame(x=qnorm(u),u=u),mapping = aes(x=x,y=u),col="#C11432") +
  geom_line(data.frame(x=qnorm(u),u=u),mapping=aes(x=x,y=u),alpha=0.5,col="#C11432") +
  ylab(TeX("$\\Phi(x)$")) +
  xlab(TeX("$x$"))

# this is more useful shown on the actual distribution
u <- seq(0.02,0.98,length.out=49)
u1 <- seq(0.0001,0.9999,length.out=998)
Zu <- c()
Zu1 <- c()
to_opt <- function(z) {
  return( (mean(pnorm((z-Z_2)/density(Z_2)$bw)) - u[j])^2)
}
for (j in 1:length(u)) {
  Zu[j] <- optim(fn=to_opt,par=1)$par
}
to_opt <- function(z) {
  return( (mean(pnorm((z-Z_2)/density(Z_2)$bw)) - u1[j])^2)
}
for (j in 1:length(u1)) {
  Zu1[j] <- optim(fn=to_opt,par=1)$par
}

ggplot() + 
  geom_line(data.frame(x=Zu1,u=u1),mapping=aes(x=x,y=u),alpha=0.5) +
  geom_point(data.frame(x=Zu,u=u),mapping = aes(x=x,y=u),col="#C11432") +
  geom_line(data.frame(x=Zu,u=u),mapping=aes(x=x,y=u),alpha=0.5,col="#C11432") +
  ylab(TeX("$s$")) +
  xlab(TeX("$\\tilde{F}^{-1}_{2|1}\\left(s\\right)$"))

# print summary of the parameters ----
N <- 5000
set.seed(12)
sims <- generate_dep_X_Y_Y_Z(N=N)

# PIT to Laplace
sims <- sims %>% mutate(Y_1=as.numeric(map(.x=X_1,.f=frechet_laplace_pit))) %>% 
  mutate(Y_2=as.numeric(map(.x=X_2,.f=frechet_laplace_pit))) %>%
  mutate(Y_3=as.numeric(map(.x=X_3,.f=frechet_laplace_pit)))
# print summary
print(xtable(par_summary(sims=sims),digits=3,include.rownames=FALSE))

# do 1000 simulations to get CI for the estimates
sumar <- list()
for (i in 1:1000) {
  set.seed(123*i)
  sims <- generate_dep_X_Y_Y_Z(N=N)
  # PIT to Laplace
  sims <- sims %>% mutate(Y_1=as.numeric(map(.x=X_1,.f=frechet_laplace_pit))) %>% 
    mutate(Y_2=as.numeric(map(.x=X_2,.f=frechet_laplace_pit))) %>%
    mutate(Y_3=as.numeric(map(.x=X_3,.f=frechet_laplace_pit)))
  sumar[[i]] <- par_summary(sims=sims)
}


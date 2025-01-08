# libraries and functions ----
library(evd)
library(tidyverse)
library(latex2exp)
library(gridExtra)
library(viridis)
library(MASS) #use dplyr::select to avoid function conflict
library(xtable)
library(gnorm)
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


# generate trivariate sample ----
N <- 50000
sims <- generate_Y(N = N) %>% link_log(dep=1/2) %>% link_log(dep=1/2) %>% 
  apply(c(1,2),FUN=frechet_laplace_pit) %>% as.data.frame()

# check it looks Laplace
# ggplot(sims %>% dplyr::select(Y1,Y2,Y3) %>% pivot_longer(everything())) + geom_density(aes(x=value),stat="density") + facet_wrap(~name)
# simulate a laplace sample
# tmp <- data.frame(Y1=qfrechet(p=seq(0.001,0.999,length.out=5000000))) %>% mutate(Y1=as.numeric(map(.x=Y1,.f=frechet_laplace_pit)))
# p2 <- ggplot(tmp %>% dplyr::select(Y1)) +  geom_density(aes(x=Y1),stat="density") + xlim(c(-6,6)) + xlab("Laplace density function") + ylab("") + ylim(c(0,1))
# # plot uniform
# p1 <- ggplot(data.frame(x=seq(-6,6,length.out=100),y=dunif(seq(-6,6,length.out=100))))+ geom_line(aes(x=x,y=y))+ xlim(c(-6,6)) + xlab("Uniform density function") + ylab("")
# grid.arrange(p1,p2,ncol=2)

lim_min <- min(sims)
lim_max <- max(sims)
vL <- frechet_laplace_pit(qfrechet(0.99))
tmp <- sims %>% mutate(above_thres= as.character(sims$Y1>vL))
grid.arrange(ggplot(tmp) + geom_point(aes(x=Y1,y=Y2,col=above_thres),size=0.5,alpha=0.5) +
               scale_color_manual(values = c("FALSE"="black","TRUE" = "#009ADA")) +
               geom_vline(xintercept=vL,color="#009ADA",linetype="dashed") +
              xlab(TeX("$Y_1$")) + ylab(TeX("$Y_2$")) + xlim(c(lim_min,lim_max)) + ylim(c(lim_min,lim_max)) +
               theme(legend.position = "none") + coord_fixed(),
             ggplot(tmp) + geom_point(aes(x=Y1,y=Y3,col=above_thres),size=0.5,alpha=0.5) +
             scale_color_manual(values = c("FALSE"="black","TRUE" = "#009ADA")) +
               geom_vline(xintercept=vL,color="#009ADA",linetype="dashed") +
             xlab(TeX("$Y_1$")) + ylab(TeX("$Y_3$")) + xlim(c(lim_min,lim_max)) + ylim(c(lim_min,lim_max))+
               theme(legend.position="none") + coord_fixed(),ncol=2)

# evaluate the Gaussian copula dependence model -----
Gaus_cop_p_est <- function(df=sims,given=1,v=0.99,v_sim=0.999,Nsim=50000) {
l <- given
res <- c(1:3)[-l]
v <- 0.99
Y_given_1_extreme <- sims %>% filter(Y3>quantile(Y3,v))
# assume a_hat and b_hat are AD case
pe <- par_est(sims,given=l)
a_hat <- pe$a
b_hat <- pe$b
# extrapolate using kernel smoothed residuals ----
Y1 <- Y_given_1_extreme[l]
Y2 <- Y_given_1_extreme[res[1]]
Y3 <- Y_given_1_extreme[res[2]]

Z2 <- Z3 <- c()
for (i in 1:length(Y1)) {
Z2[i] <-   (Y2[i]-a_hat[1]*Y1[i])/(Y1[i]^b_hat[1]) 
Z3[i] <-   (Y3[i]-a_hat[2]*Y1[i])/(Y1[i]^b_hat[2]) 
}

# simulate from observed residuals and Gaussian copula with kernel smoothed density
# calculate the normal using the PIT
Z2 <- as.numeric(unlist(Z2))
Z3 <- as.numeric(unlist(Z3))
ZN2 <- qnorm(F_smooth_Z(Z2))
ZN3 <- qnorm(F_smooth_Z(Z3))
rho_hat <- cor(ZN2,ZN3)
Nsim <- 500
ZN <- mvrnorm(n=Nsim,mu=c(0,0),Sigma=matrix(c(1,rho_hat,rho_hat,1),2,2))
# transform back to original margins
Z_star <- norm_to_orig(ZN=ZN,emp_res = Z)
# ggplot(data.frame(Z2=Z_star$X1,Z3=Z_star$X2))+geom_point(aes(x=Z2,y=Z3),size=0.9,alpha=0.5,col="#C11432")+ 
  # xlab(TeX("$Z_2$")) + ylab(TeX("$Z_3$")) + xlim(min(Z_star),max(Z_star)) + ylim(min(Z_star),max(Z_star)) + coord_fixed()

# generate X_1 from Laplace distribution above 0.9 quantile
set.seed(12)
Nsim <- 50000
Y1_gen <- -log(2*(1-0.999)) + rexp(Nsim)
Gen_Y1 <- data.frame(Y1=Y1_gen)

# for each Y, generate a residual and calculate Y_2
Y1 <- Gen_Y1$Y1
#Z_gen <- sample(Z2,N,replace=TRUE) +rnorm(N,mean=0,sd=density(Z)$bw) # plus noise
Y2 <- pe$a[1]*Y1 + Y1^pe$b[1] *Z_star[,1]
Y3 <- pe$a[2]*Y1 + Y1^pe$b[2] *Z_star[,2]
Gen_Y1 <- Gen_Y1 %>% mutate(Y2=Y2, Y3=Y3) %>% mutate(sim=rep("model",Nsim))
# generate Y_1 (extrapolate so above largest observed value)

# #plot
# Gen_orig <- rbind(Gen_Y1,Y_given_1_extreme %>% dplyr::select(Y1,Y2,Y3) %>% mutate(sim=rep("original_laplace",500)))
# ggplot(Gen_orig) + geom_point(aes(x=Y1,y=Y2,col=sim),alpha=0.5) +
#   scale_color_manual(values = c("original_laplace"="black","model" = "#C11432"))

# specify threshold for Laplace margin
#v_l <- c(5,12,5,12)
vL <- frechet_laplace_pit(qfrechet(0.999))

# calculate empirical probability by simulating Y_2 from the model
p <- ((Gen_Y1 %>% filter(Y1>vL,Y2>vL,Y3>vL) %>% dim())[1]/Nsim)*(1-0.999)
# calculate CI
CI <- c(p-(1.96*(p*(1-p)/Nsim)^(0.5)),p+(1.96*(p*(1-p)/Nsim)^(0.5)))
}
ran_bern <- rbinom(n=1000,size = 50000,p=p)/50000
# density(ran_bern) %>% plot()
ggplot(data.frame(x=ran_bern)) + geom_density(aes(x=x),stat="density") 

# need to log Fréchet margins to get the logistic model margin
p <- evd::pbvevd(c(v[2],v[4]),dep=dep[1],model="log") -
  evd::pbvevd(c(v[1],v[4]),dep=dep[1],model="log") -
  evd::pbvevd(c(v[2],v[3]),dep=dep[1],model="log") +
  evd::pbvevd(c(v[1],v[3]),dep=dep[1],model="log")
x <- seq(min(sims)-1,max(sims)+1,length.out=10000)
yl <- cond_quantile(x,Z=Z2,q=0.025,a_hat=a_hat[1],b_hat=b_hat[1])
ym <- cond_quantile(x,Z=Z2,q=0.5,a_hat=a_hat[1],b_hat=b_hat[1])
yp <- cond_quantile(x,Z=Z2,q=0.975,a_hat=a_hat[1],b_hat=b_hat[1])

yl3 <- cond_quantile(x,Z3,q=0.025,a_hat=a_hat[2],b_hat=b_hat[2])
ym3 <- cond_quantile(x,Z3,q=0.5,a_hat=a_hat[2],b_hat=b_hat[2])
yp3 <- cond_quantile(x,Z3,q=0.975,a_hat=a_hat[2],b_hat=b_hat[2])

# plot again with data
grid.arrange(ggplot(tmp %>% filter(above_thres=="TRUE")) + geom_point(aes(x=Y1,y=Y2,col=above_thres),size=0.8,alpha=0.5) +
               scale_color_manual(values = c("FALSE"="black","TRUE" = "#009ADA")) +
               geom_vline(xintercept=vL,color="#009ADA",linetype="dashed") +
               geom_line(data=data.frame(x=x,y=yl),aes(x=x,y=y),col="#C11432",linetype="dashed")+
               geom_line(data=data.frame(x=x,y=ym),aes(x=x,y=y),col="#C11432")+
               # geom_ribbon(data=data.frame(x=x,yl=yl,yp=yp),aes(x=x,ymin=yl,ymax=yp), fill="#009ADA", alpha=0.2) +
               geom_line(data=data.frame(x=x,y=yp),aes(x=x,y=y),col="#C11432",linetype="dashed")+
               xlab(TeX("$Y_1$")) + ylab(TeX("$Y_2$")) + xlim(c(vL,lim_max)) + ylim(c(lim_min+8,lim_max)) +
               theme(legend.position = "none"),
             ggplot(tmp %>% filter(above_thres=="TRUE")) + geom_point(aes(x=Y1,y=Y3,col=above_thres),size=0.8,alpha=0.5) +
               geom_line(data=data.frame(x=x,y=yl3),aes(x=x,y=y),col="#C11432",linetype="dashed")+
               geom_line(data=data.frame(x=x,y=ym3),aes(x=x,y=y),col="#C11432")+
               geom_line(data=data.frame(x=x,y=yp3),aes(x=x,y=y),col="#C11432",linetype="dashed")+
               scale_color_manual(values = c("FALSE"="black","TRUE" = "#009ADA")) +
               geom_vline(xintercept=vL,color="#009ADA",linetype="dashed") +
               xlab(TeX("$Y_1$")) + ylab(TeX("$Y_3$")) + xlim(c(vL,lim_max)) + ylim(c(lim_min+8,lim_max))+
               theme(legend.position="none"),ncol=2)

# simulate residuals Zs from the empirical distribution
Z <- data.frame(Z2,Z3) # dataframe of observed residuals
Zs <- as.data.frame(matrix(ncol=2,nrow=0))
names(Zs) <- c("Z2","Z3")
for (i in 1:1000) {
  Zs[i,] <- Z[sample(1:nrow(Z),1,replace=TRUE),]
}
# plot simulated residuals
ggplot(Zs)+geom_point(aes(x=Z2,y=Z3),size=0.9,alpha=0.5,col="#C11432") + xlab(TeX("$Z_2$")) + ylab(TeX("$Z_3$")) +
  xlim(min(Zs),max(Zs)) + ylim(min(Zs),max(Zs)) + coord_fixed()

# create a loop to count points in each region ----
q <- seq(0.05,0.95,by=0.05)
tmp <- rep(NA,length(Y1))
# bottom edge case
y1 <- cond_quantile(x,Z,q=min(q),a_hat=a_hat,b_hat=b_hat)
A <- c(min(x),y1[which.min(x)])
B <- c(max(x),y1[which.max(x)])
D <- c(min(x),-10)
C <- c(max(x),-10)
df_tmp <- data.frame(matrix(c(A,B,C,D),ncol=2,byrow=TRUE))
no <-   data.frame(pointsInPolygon(Y_given_1_extreme[,4:5],df_tmp))
tmp[pointsInPolygon(Y_given_1_extreme[,4:5],df_tmp)] <- paste0(1)
for (i in 1:(length(q)-1)) {
  y1 <- cond_quantile(x,Z,q=q[i],a_hat=a_hat,b_hat=b_hat)
  y2 <- cond_quantile(x,Z,q=q[i+1],a_hat=a_hat,b_hat=b_hat)
  A <- c(min(x),y1[which.min(x)])
  B <- c(max(x),y1[which.max(x)])
  D <- c(min(x),y2[which.min(x)])
  C <- c(max(x),y2[which.max(x)])
  df_tmp <- data.frame(matrix(c(A,B,C,D),ncol=2,byrow=TRUE))
  no <-   cbind(no,pointsInPolygon(Y_given_1_extreme[,4:5],df_tmp))
  tmp[pointsInPolygon(Y_given_1_extreme[,4:5],df_tmp)] <- paste0(i+1)
}
# top edge case
y1 <- cond_quantile(x,Z,q=max(q),a_hat=a_hat,b_hat=b_hat)
A <- c(min(x),y1[which.min(x)])
B <- c(max(x),y1[which.max(x)])
D <- c(min(x),20)
C <- c(max(x),20)
df_tmp <- data.frame(matrix(c(A,B,C,D),ncol=2,byrow=TRUE))
no <- cbind(no,pointsInPolygon(Y_given_1_extreme[,4:5],df_tmp))
rowSums(no) # check it is all ones
colSums(no) %>% as.vector()
tmp[pointsInPolygon(Y_given_1_extreme[,4:5],df_tmp)] <- paste0(length(q)+1)

ggplot() + geom_point(data=Y_given_1_extreme %>% mutate(check=tmp),aes(x=Y1,y=Y2,col=check),alpha=0.5)+
  geom_polygon(data=df_tmp,aes(x=X1,y=X2),fill="#009ada",alpha=0.2)

tmp_df %>% filter(Lik>quantile(Lik,0)) %>% 
  ggplot(aes(x = Var1, y = Var2, z = Lik)) + 
  geom_tile( aes(fill = Lik)) + 
  #scale_fill_viridis() +  
  geom_contour(color = "black", alpha=0.5) +
  xlab(TeX("$\\alpha$")) +
  ylab(TeX("$\\beta$"))

ggplot(Y_given_1_extreme %>% select(Y1,Y2_sim,Y2,Y3) %>% pivot_longer(everything())) + geom_density(aes(x=value),,stat="density") +
  facet_wrap(~name)

Y2_sim <- bind_rows(Y_not_1_extreme %>% mutate(Y_2_sim=Y_2),Y_given_1_extreme) %>% 
  mutate(X_2_sim=rep(0,N))

X2_simdf <- Y2_sim %>% 
  mutate(X_2=as.numeric(map(.x=Y_2,.f=laplace_frechet_pit))) %>%
  mutate(X_2_sim=as.numeric(map(.x=Y_2_sim,.f=laplace_frechet_pit)))

X2_simdf <- X2_simdf %>% mutate(v=c(rep("below_threshold",N*v),rep("above_threshold",N*(1-v))))

grid.arrange(ggplot(X2_simdf) + geom_point(aes(x=X1,y=X2,col=v),alpha=0.5),
             ggplot(X2_simdf) + geom_point(aes(x=X2,y=X3,col=v),alpha=0.5),
             ggplot(X2_simdf) + geom_point(aes(x=X1,y=X2_sim,col=v),alpha=0.5),
             ggplot(X2_simdf) + geom_point(aes(x=X2_sim,y=X3,col=v),alpha=0.5),ncol=2)

grid.arrange(ggplot(X2_simdf) + geom_point(aes(x=Y1,y=Y2,col=v),alpha=0.5),
             ggplot(X2_simdf) + geom_point(aes(x=Y2,y=Y3,col=v),alpha=0.5),
             ggplot(X2_simdf) + geom_point(aes(x=Y1,y=Y2_sim,col=v),alpha=0.5),
             ggplot(X2_simdf) + geom_point(aes(x=Y2_sim,y=Y3,col=v),alpha=0.5),ncol=2)
# plot only extremes
grid.arrange(ggplot(X_2_simdf %>% filter(v=="above_threshold")) + geom_point(aes(x=X_1,y=X_2),alpha=0.5),
             ggplot(X_2_simdf%>% filter(v=="above_threshold")) + geom_point(aes(x=X_2,y=X_3),alpha=0.5),
             ggplot(X_2_simdf%>% filter(v=="above_threshold")) + geom_point(aes(x=X_1,y=X_2_sim),alpha=0.5),
             ggplot(X_2_simdf%>% filter(v=="above_threshold")) + geom_point(aes(x=X_2_sim,y=X_3),alpha=0.5),ncol=2)

grid.arrange(ggplot(X_2_simdf %>% filter(v=="above_threshold")) + geom_point(aes(x=Y1,y=Y_2),alpha=0.5),
             ggplot(X_2_simdf%>% filter(v=="above_threshold")) + geom_point(aes(x=Y2,y=Y_3),alpha=0.5),
             ggplot(X_2_simdf%>% filter(v=="above_threshold")) + geom_point(aes(x=Y1,y=Y_2_sim),alpha=0.5),
             ggplot(X_2_simdf%>% filter(v=="above_threshold")) + geom_point(aes(x=Y2_sim,y=Y_3),alpha=0.5),ncol=2)


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

# start from the beginning ----
# generate trivariate sample
N <- 5000000
sims <- generate_Y(N = N) %>% link_log(dep=1/2) %>% link_log(dep=1/2) %>% 
  apply(c(1,2),FUN=frechet_laplace_pit) %>% as.data.frame()
# filter for Y_1 being extreme -----
v <- 0.99
Y_given1extreme <- sims %>% filter(Y1>quantile(Y1,v))
Y_not1extreme <- sims %>% filter(Y1<quantile(Y1,v))

d <- 3
opt <- list()
for (i in c(2,3)) {
opt[[i-1]] <- optim(par=c(0.6,0.2,0,1),fn = Y_likelihood,df=Y_given1extreme,given=1,sim=i,control = list(fnscale=-1))
}
a_hat <- c(opt[[1]]$par[1],opt[[2]]$par[1])
b_hat <- c(opt[[1]]$par[2],opt[[2]]$par[2])

# generate residual Z ----
Y1 <- Y_given1extreme[,1]
Y2 <- Y_given1extreme[,2]
Y3 <- Y_given1extreme[,3]

Z2 <- (Y2-a_hat[1]*Y1)/(Y1^b_hat[1])
Z3 <- (Y3-a_hat[2]*Y1)/(Y1^b_hat[2])
plot(Y1,Z2)
plot(Y1,Z3)

# calculate the normal using the PIT
Z_N2 <- qnorm(F_smooth_Z(Z2))
Z_N3 <- qnorm(F_smooth_Z(Z3))
rho_hat <- cor(Z_N2,Z_N3)
Nsim <- 1000
Z_N <- mvrnorm(n=Nsim,mu=c(0,0),Sigma=matrix(c(1,rho_hat,rho_hat,1),2,2))
Z <- data.frame(Z2,Z3)

# transform back to original margins
Z_star <- norm_to_orig(ZN=Z_N,emp_res = Z)

U <- runif(Nsim)
Y1_gen <- -log(2*(1-0.9999)) + rexp(Nsim)
Gen_Y1 <- data.frame(Y1=Y1_gen)

# calculate empirical probability
vL <- frechet_laplace_pit(qfrechet(0.9999))
(sims %>% filter(Y1>vL,Y2>vL,Y3>vL) %>% nrow())/N

# for each Y, generate a residual and calculate Y_2
Y1 <- Gen_Y1$Y1
Y2 <- a_hat[1]*Y1 + Y1^b_hat[1] *Z_star[,1]
Y3 <-  a_hat[2]*Y1 + Y1^b_hat[2] *Z_star[,2]
Gen_Y1 <- Gen_Y1 %>% mutate(Y2=Y2,Y3=Y3) %>% mutate(sim=rep("model",Nsim))
# generate Y_1 (extrapolate so above largest observed value)

#plot
Gen_orig <- rbind(Gen_Y1,Y_given_1_extreme %>% dplyr::select(Y1,Y2,Y3) %>% mutate(sim=rep("original_laplace",500)))
p1 <- ggplot(Gen_orig) + geom_point(aes(x=Y1,y=Y2,col=sim),alpha=0.5) + 
  scale_color_manual(values = c("original_laplace"="black","model" = "#C11432")) 
p2 <- ggplot(Gen_orig) + geom_point(aes(x=Y2,y=Y3,col=sim),alpha=0.5) + 
  scale_color_manual(values = c("original_laplace"="black","model" = "#C11432")) 
p3 <- ggplot(Gen_orig) + geom_point(aes(x=Y1,y=Y3,col=sim),alpha=0.5) + 
  scale_color_manual(values = c("original_laplace"="black","model" = "#C11432")) 
grid.arrange(p1,p2,p3,ncol=3)

Z_comp <- Z_star %>% mutate(Compare=rep("Optimise_all",Nsim))
Z_comp <- rbind(Z_comp,Z_star%>% mutate(Compare=rep("Linear_segments",Nsim)))
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
  return( (mean(pnorm((z-Z2)/density(Z2)$bw)) - u[j])^2)
}
for (j in 1:length(u)) {
  Zu[j] <- optim(fn=to_opt,par=1)$par
}
to_opt <- function(z) {
  return( (mean(pnorm((z-Z2)/density(Z2)$bw)) - u1[j])^2)
}
for (j in 1:length(u1)) {
  Zu1[j] <- optim(fn=to_opt,par=1)$par
}

ggplot() + 
  geom_line(data.frame(x=Zu1,u=u1),mapping=aes(x=x,y=u),alpha=0.5) +
  geom_point(data.frame(x=Zu,u=u),mapping = aes(x=x,y=u),col="#C11432") +
  geom_line(data.frame(x=Zu,u=u),mapping=aes(x=x,y=u),alpha=0.5,col="#C11432") +
  ylab(TeX("$s$")) +
  xlab(TeX("$\\tilde{G}^{-1}_{2}\\left(s\\right)$"))

# print summary of the parameters ----
print(xtable(par_summary(sims=sims),digits=3,include.rownames=FALSE))

# do 1000 simulations to get CI for the estimates
sumar <- list()
for (i in 1:100) {
  set.seed(12*i)
  sims <- generate_Y(N=N) %>% link_log(dep=1/2) %>%
    link_log(dep=1/2) %>%
    apply(c(1,2),FUN=frechet_laplace_pit) %>% as.data.frame()
  sumar[[i]] <- par_summary(sims=sims)
}

# plot boxplots for a,b,mu,sig estimated values
muboot <- c()
sigboot <- c()
aboot <- c()
bboot <- c()
rhoboot <- c()
Ys <- rep(c("Y1","Y2","Y3"),length(sumar))
res <- rep(c("Z_21","Z_31","Z_12","Z_32","Z_13","Z_23"),length(sumar))
for (i in 1:length(sumar)) {
  aboot <- append(aboot,unname(unlist(sumar[[i]][7,])))
  bboot <- append(bboot,unname(unlist(sumar[[i]][8,])))
  muboot <- append(muboot,unname(unlist(sumar[[i]][9,])))
  sigboot <- append(sigboot,unname(unlist(sumar[[i]][10,])))
  rhoboot <- append(rhoboot,unname(unlist(sumar[[i]][11,c(2,4,6)])))
}

p1 <- ggplot(data.frame(x=factor(res,levels=c("Z_21","Z_31","Z_12","Z_32","Z_13","Z_23")),y=muboot),aes(x=x,y=y))+
  geom_boxplot(fill= c(rep("#C11432",2),rep("#66A64F",2),rep("#009ADA",2))) + xlab("Observed residuals") + ylab(TeX("$\\hat{\\mu}$")) 
p2 <- ggplot(data.frame(x=factor(res,levels=c("Z_21","Z_31","Z_12","Z_32","Z_13","Z_23")),y=sigboot),aes(x=x,y=y))+
  geom_boxplot(fill= c(rep("#C11432",2),rep("#66A64F",2),rep("#009ADA",2))) + xlab("Observed residuals") + ylab(TeX("$\\hat{\\sigma}$")) 
p3 <- ggplot(data.frame(x=factor(res,levels=c("Z_21","Z_31","Z_12","Z_32","Z_13","Z_23")),y=aboot),aes(x=x,y=y))+
  geom_boxplot(fill= c(rep("#C11432",2),rep("#66A64F",2),rep("#009ADA",2))) + xlab("Observed residuals") + ylab(TeX("$\\hat{\\alpha}$")) 
p4 <- ggplot(data.frame(x=factor(res,levels=c("Z_21","Z_31","Z_12","Z_32","Z_13","Z_23")),y=bboot),aes(x=x,y=y))+
  geom_boxplot(fill= c(rep("#C11432",2),rep("#66A64F",2),rep("#009ADA",2))) + xlab("Observed residuals") + ylab(TeX("$\\hat{\\beta}$")) 
grid.arrange(p3,p4,p1,p2,ncol=2)
p1 <- ggplot(data.frame(x=factor(res,levels=c("Z_21","Z_31","Z_12","Z_32","Z_13","Z_23")),y=aboot),aes(x=y,fill=x))+
  geom_density(alpha=0.5) + xlab("Observed residuals") + ylab(TeX("$\\hat{\\alpha}$")) +
  scale_fill_manual(values=c(rep("#C11432",1),"#009ADA","#C11432",rep("#66A64F",1),rep("#009ADA",1),"#66A64F"))
p2 <- ggplot(data.frame(x=factor(res,levels=c("Z_21","Z_31","Z_12","Z_32","Z_13","Z_23")),y=bboot),aes(x=y,fill=x))+
  geom_density(alpha=0.5) + xlab("Observed residuals") + ylab(TeX("$\\hat{\\beta}$")) +
  scale_fill_manual(values=c(rep("#C11432",1),"#009ADA","#C11432",rep("#66A64F",1),rep("#009ADA",1),"#66A64F"))
p3 <- ggplot(data.frame(x=factor(res,levels=c("Z_21","Z_31","Z_12","Z_32","Z_13","Z_23")),y=muboot),aes(x=y,fill=x))+
  geom_density(alpha=0.5) + xlab("Observed residuals") + ylab(TeX("$\\hat{\\mu}$")) +
  scale_fill_manual(values=c(rep("#C11432",1),"#009ADA","#C11432",rep("#66A64F",1),rep("#009ADA",1),"#66A64F"))
p4 <- ggplot(data.frame(x=factor(res,levels=c("Z_21","Z_31","Z_12","Z_32","Z_13","Z_23")),y=sigboot),aes(x=y,fill=x))+
  geom_density(alpha=0.5) + xlab("Observed residuals") + ylab(TeX("$\\hat{\\sigma}$")) +
  scale_fill_manual(values=c(rep("#C11432",1),"#009ADA","#C11432",rep("#66A64F",1),rep("#009ADA",1),"#66A64F"))
grid.arrange(p1,p2,p3,p4,ncol=2)

ggplot(data.frame(x=factor(Ys,levels=c("Y_1","Y_2","Y_3")),y=rhoboot),aes(x=x,y=y))+
  geom_boxplot(fill= c(rep("#C11432",1),rep("#66A64F",1),rep("#009ADA",1))) + xlab("Conditional variable") + ylab(TeX("$\\hat{\\rho}$"))

# plot the residual pairs
tmp_df <- plot_residual(sims=sims)
l <- c(min(tmp_df$Z_2,tmp_df$Z_3),max(tmp_df$Z_2,tmp_df$Z_3))
ggplot(tmp_df) + geom_point(aes(x=Z_2,y=Z_3)) + facet_wrap(~given) + xlim(l) + ylim(l)
rm(l)

# plot the normalised residuals
l <- c(min(tmp_df$Z_N_2,tmp_df$Z_N_3),max(tmp_df$Z_N_2,tmp_df$Z_N_3))
ggplot(tmp_df) + geom_point(aes(x=Z_N_2,y=Z_N_3)) + facet_wrap(~given)  + xlim(l) + ylim(l)
rm(tmp_df,l)
# plot the residuals after the transformation


# plot residuals
tmp_df <- rbind(plot_residual(sims=sims) %>% mutate(Estimate="ab_estimated"),plot_residual_trueab(sims=sims) %>% mutate(Estimate="ab_true"))
p1 <- tmp_df %>% filter(given=="1") %>% ggplot() + geom_point(mapping = aes(x=Z_2,y=Z_3,color=Estimate)) +
ggtitle("Given 1")+ scale_color_manual(values=c("#C11432","black"))
p2 <- tmp_df %>% filter(given=="2") %>% ggplot() + geom_point(mapping = aes(x=Z_2,y=Z_3,color=Estimate)) +
  ggtitle("Given 2")+ scale_color_manual(values=c("#66A64F","black"))
p3 <- tmp_df %>% filter(given=="3") %>% ggplot() + geom_point(mapping = aes(x=Z_2,y=Z_3,color=Estimate)) +
  ggtitle("Given 3")+ scale_color_manual(values=c("#009ADA","black"))
grid.arrange(p1,p2,p3,ncol=3)

# tmp_df <- rbind(plot_residual(sims=sims) %>% mutate(Estimate="ab_estimated"),plot_residual_trueab(sims=sims) %>% mutate(Estimate="ab_true"))
p1 <- tmp_df %>% filter(given=="1") %>% ggplot() + geom_point(mapping = aes(x=Z_N_2,y=Z_N_3,color=Estimate)) +
  ggtitle("Given 1")+ scale_color_manual(values=c("#C11432","black"))
p2 <- tmp_df %>% filter(given=="2") %>% ggplot() + geom_point(mapping = aes(x=Z_N_2,y=Z_N_3,color=Estimate)) +
  ggtitle("Given 2")+ scale_color_manual(values=c("#66A64F","black"))
p3 <- tmp_df %>% filter(given=="3") %>% ggplot() + geom_point(mapping = aes(x=Z_N_2,y=Z_N_3,color=Estimate)) +
  ggtitle("Given 3")+ scale_color_manual(values=c("#009ADA","black"))
grid.arrange(p1,p2,p3,ncol=3)

p1 <- plot_simulated(sims=sims,given=1)
p2 <- plot_simulated(sims=sims,given=2)
p3 <- plot_simulated(sims=sims,given=3)
grid.arrange(p1,p2,p3,ncol=1)

v_l <- rep(frechet_laplace_pit( qfrechet(0.999)),2)
# transform to frechet margins
# v <- log(sapply(X = v_l,FUN = laplace_frechet_pit))
annotate('rect', xmin=3, xmax=5, ymin=3, ymax=7, alpha=.2, fill='red')

# calculate empirical probability by simulating Y_2 from the model
set.seed(12)
sim_val <- plot_simulated(sims=sims,given=3)
giv1 <- c(((sim_val %>% filter(sim=="model") %>% filter(Y_1>v_l[1],Y_2>v_l[2],Y_3>v_l[1]) %>% dim())[1]/1000)*(1-0.999))
giv1 <- c(((sim_val %>% filter(sim=="model") %>% filter(Y_1>v_l[1],Y_2>v_l[2],Y_3>v_l[1]) %>% dim())[1]/1000)*(1-0.999),
           ((sim_val %>% filter(sim=="model") %>% filter(Y_2>v_l[1],Y_3>v_l[2]) %>% dim())[1]/1000)*(1-0.999),
           ((sim_val %>% filter(sim=="model") %>% filter(Y_1>v_l[1],Y_3>v_l[2]) %>% dim())[1]/1000)*(1-0.999)
)
sim_val <- plot_simulated(sims=sims,given=2)
giv2 <-  c(((sim_val %>% filter(sim=="model") %>% filter(Y_1>v_l[1],Y_2>v_l[2],Y_3>v_l[1]) %>% dim())[1]/1000)*(1-0.999))

giv2 <- c(((sim_val %>% filter(sim=="model") %>% filter(Y_1>v_l[1],Y_2>v_l[2]) %>% dim())[1]/1000)*(1-0.999),
           ((sim_val %>% filter(sim=="model") %>% filter(Y_2>v_l[1],Y_3>v_l[2]) %>% dim())[1]/1000)*(1-0.999),
           ((sim_val %>% filter(sim=="model") %>% filter(Y_1>v_l[1],Y_3>v_l[2]) %>% dim())[1]/1000)*(1-0.999)
)
sim_val <- plot_simulated(sims=sims,given=3)
giv3 <-  c(((sim_val %>% filter(sim=="model") %>% filter(Y_1>v_l[1],Y_2>v_l[2],Y_3>v_l[1]) %>% dim())[1]/1000)*(1-0.999))

giv3 <- c(((sim_val %>% filter(sim=="model") %>% filter(Y_1>v_l[1],Y_2>v_l[2]) %>% dim())[1]/1000)*(1-0.999),
           ((sim_val %>% filter(sim=="model") %>% filter(Y_2>v_l[1],Y_3>v_l[2]) %>% dim())[1]/1000)*(1-0.999),
           ((sim_val %>% filter(sim=="model") %>% filter(Y_1>v_l[1],Y_3>v_l[2]) %>% dim())[1]/1000)*(1-0.999)
)

tmp_df <- data.frame(giv1,giv2,giv3)
# need to log Fréchet margins to get the logistic model margin

p <- 1 - 
  0.999 -
  0.999 +
  evd::pbvevd(c(log(qfrechet(0.999)),log(qfrechet(0.999))),dep=dep[1],model="log")

# calculate CI
p <- giv_1[1]
CI <- c(p-(1.96*(p*(1-p)/100000)^(0.5)),p+(1.96*(p*(1-p)/100000)^(0.5)))

tmp_CI <- data.frame(giv_1=rep(NA,1), giv_2=rep(NA,1),giv_3=rep(NA,1))
for (i in 1:1) {
  for (j in 1:3) {
    p <- tmp_df[i,j]
    tmp_CI[i,j] <- paste0("(",p-(1.96*(p*(1-p)/100000)^(0.5)),", ",p+(1.96*(p*(1-p)/100000)^(0.5)),")")
  }
}

sims[(sims$Y1>quantile(sims$Y1,v1) & sims$Y2>quantile(sims$Y2,v) & sims$Y3>quantile(sims$Y_3,v1)),] %>% glimpse()
sims[(sims$Y1>quantile(sims$Y1,v) & sims$Y2>quantile(sims$Y2,v)),] %>% glimpse()

# calculate the exact probability
to_opt <- function(z) {
  (  (  y^(-(1/a)+1)*(y^(-1/a)+z^(-1/a))^(a-1)*exp(-(y^(-1/a)+z^(-1/a))^a)*exp(1/y)  )-Unif)^2
}

# integrand <- function(y) {
#   (1-  y^(-(1/a)+1)*(y^(-1/a)+x^(-1/a))^(a-1)*exp(-(y^(-1/a)+x^(-1/a))^a)*exp(1/y)  )*
#     y^(-2)*exp(-1/y)*
#     (1-  y^(-(1/a)+1)*(y^(-1/a)+z^(-1/a))^(a-1)*exp(-(y^(-1/a)+z^(-1/a))^a)*exp(1/y)  )
# }

integrand <- function(y) {
  (1-  (1+(x/y)^(-1/a))^(a-1)*exp( y^(-1)*( 1- ( 1+(x/y)^(-1/a) )^a) )  )*
    y^(-2)*exp(-1/y)*
    (1-  (1+(z/y)^(-1/a))^(a-1)*exp( y^(-1)*( 1- ( 1+(z/y)^(-1/a) )^a) ) )
}

x <- z <- qfrechet(0.995)
integrate(integrand,lower=qfrechet(0.995),upper=Inf)

s <- seq(0.9,0.997,length.out=50)
cdf <- c()
cdf_emp <- c()
for (i in 1:length(s)) {
  x <- z <- qfrechet(s[i])
  v <- s[i]
  cdf[i] <- integrate(integrand,lower=qfrechet(s[i]),upper=Inf)$value
 cdf_emp[i] <- nrow(sims[(sims$Y_1>quantile(sims$Y_1,v) & sims$Y_2>quantile(sims$Y_2,v) & sims$Y_3>quantile(sims$Y_3,v)),])/500000
  
}

nrow(sims[(sims$Y_1>quantile(sims$Y_1,v) & sims$Y_2>quantile(sims$Y_2,v) & sims$Y_3>quantile(sims$Y_3,v)),])/500000

ggplot(data.frame(x=qfrechet(s),cdf,cdf_emp)) + geom_line(aes(x=x,y=cdf),col="#C11432",alpha=0.6) + 
  geom_line(aes(x=x,y=cdf_emp),col="#009ADA",alpha=0.6) + geom_point(aes(x=x,y=cdf),col="#C11432",alpha=0.6) + 
  geom_point(aes(x=x,y=cdf_emp),col="#009ADA",alpha=0.6)

generate_Y(N=N) %>% link_log(dep=1/2) %>% link_norm(dep=1/2) %>% plot_pairs()
set.seed(12)
sims_tmp <- generate_Y(N=N) %>% link_log(dep=1/2) %>% link_log(dep=1/2) 
sims <- sims_tmp %>% apply(c(1,2),FUN=frechet_laplace_pit) %>% as.data.frame()
names(sims) <- c("Y2","Y1","Y3")
grid.arrange(ggplot(sims) + geom_point(aes(x=Y1,y=Y2),alpha=0.5,size=0.1,col="#009ADA"),
             ggplot(sims) + geom_point(aes(x=Y2,y=Y3),alpha=0.5,size=0.1,col="#009ADA"),
             ggplot(sims) + geom_point(aes(x=Y1,y=Y3),alpha=0.5,size=0.1),ncol=3)

U <- runif(N)
X <- (-log(U))^(-1)
U <- runif(N)
Z <- (-log(U))^(-1)
X <- sims_tmp[,2]
Z <- sims_tmp[,3]

to_opt <- function(y) {
  (  (  exp(-(x^(-1/a) +y^(-1/a)+z^(-1/a))^a + (x^(-1/a)+z^(-1/a))^a )*
          (x^(-1/a) +y^(-1/a)+z^(-1/a))^(a-2)*((1-a)/a+(x^(-1/a) +y^(-1/a)+z^(-1/a))^a)/
          ((x^(-1/a) +z^(-1/a))^(a-2)*((1-a)/a+(x^(-1/a) +z^(-1/a))^a))
          )-Unif)^2
}
y <- c()
for (i in 1:N){
  Unif <- runif(1) # generate U
  x <- X[i]
  z <- Z[i]
  y[i] <- optim(par=1,fn=to_opt,lower=0,upper=10^6,method="Brent")$par
}

sims1 <- data.frame(Y1=X,Y2=y,Y3=Z)%>% apply(c(1,2),FUN=frechet_laplace_pit) %>% as.data.frame()
l <- min(sims,sims1)
u <- max(sims,sims1)
grid.arrange(ggplot(sims) + geom_point(aes(x=Y1,y=Y2),alpha=0.5,size=0.1,col="#009ADA") + xlim(c(l,u))+ ylim(c(l,u)) + xlab(TeX("$Y_1$"))+ ylab(TeX("$Y_2$")),
             ggplot(sims) + geom_point(aes(x=Y2,y=Y3),alpha=0.5,size=0.1,col="#009ADA")+ xlim(c(l,u))+ ylim(c(l,u))+ xlab(TeX("$Y_2$"))+ ylab(TeX("$Y_3$")),
             ggplot(sims) + geom_point(aes(x=Y1,y=Y3),alpha=0.5,size=0.1)+ xlim(c(l,u))+ ylim(c(l,u))+ xlab(TeX("$Y_1$"))+ ylab(TeX("$Y_3$")),
             ggplot(sims1) + geom_point(aes(x=Y1,y=Y2),alpha=0.5,size=0.1,col="#009ADA")+ xlim(c(l,u))+ ylim(c(l,u))+ xlab(TeX("$Y_1$"))+ ylab(TeX("$Y_2^*$")),
             ggplot(sims1) + geom_point(aes(x=Y2,y=Y3),alpha=0.5,size=0.1,col="#009ADA")+ xlim(c(l,u))+ ylim(c(l,u))+ xlab(TeX("$Y_2^*$"))+ ylab(TeX("$Y_3$")),
             ggplot(sims1) + geom_point(aes(x=Y1,y=Y3),alpha=0.5,size=0.1)+ xlim(c(l,u))+ ylim(c(l,u))+ xlab(TeX("$Y_1$"))+ ylab(TeX("$Y_3$")),ncol=3)
#(x2*x3)^{-(1/a)-1}
integrand_subs <- function(t) {
t^(-2/a)* (1+(x2*t)^(-1/a))^(a-2) *(1+(x3*t)^(-1/a))^(a-2)*
    (((1-a)/a) + t*(1+(x2*t)^(-1/a))^a )*(((1-a)/a) + t*(1+(x3*t)^(-1/a))^a )*
    exp(t*( 1-(1+(x2*t)^(-1/a))^a-(1+(x3*t)^(-1/a))^a  ))
}
\texttt{x2*10^1.4}
integral <- function(s) {
  s^(2/a-2)* (1+(x2/s)^(-1/a))^(a-2) *(1+(x3/s)^(-1/a))^(a-2)*
    (((1-a)/a) + s^(-1)*(1+(x2/s)^(-1/a))^a )*(((1-a)/a) + s^(-1)*(1+(x3/s)^(-1/a))^a )*
    exp(s^(-1)*( 1-(1+(x2/s)^(-1/a))^a-(1+(x3/s)^(-1/a))^a  ))
}
U <- runif(1)
X2 <- sims_tmp[,2]
X3 <- sims_tmp[,3]
x2 <- X2[1]
x3 <- X3[1]
u <- integrate(integral,lower=0,upper=qfrechet(0.9999))$value
l <- integrate(integral,lower=0,upper=u)$value
U
u/l

U <- seq(from=0.0001,to=0.99999,length.out=1998)
int <- c()
for (i in 1:1998){
  Unif <- U[i]
  x2 <- max(X2)
  x3 <- X3[which.max(X2)]
  # x2 <- X2[4]
  # x3 <- X3[4]
  # u[i] <- integrate(integral,lower=0,upper=qfrechet(Unif))$value
  # l[i] <- integrate(integral,lower=0,upper=10^4)$value
  d[i] <- integral(qfrechet(Unif))
  # x <- c()
  # for (j in 1:1998) {
  #   x[j] <- integral(s=qfrechet(U[j]))
  # }
}
  y <- c()
  for (j in 1:1997) {
    y[j] <- abs(d[j+1]+d[j])
  }

  u <-  qfrechet(U[which(y==y[y<0.0000001][1])])
  if (length(u)==0) {
    u <- qfrechet(max(U))
  }
 # int[i] <-  u[i]/l[i]

plot(qfrechet(U),d,xlim = c(0,u))
plot(qfrechet(U),d)
plot(y)

to_opt <- function(x) {
  (  (  integrate(integral,lower=0,upper=x)$value/l[i] )-Unif)^2
}
x4 <-c()
for (i in 1:N){
  Unif <- runif(1) # generate U
  x2 <- X2[i]
  x3 <- X3[i]
  # x <- c()
  # for (j in 1:1998) {
  #   x[j] <- integral(s=1/qfrechet(U[j]))
  # }
  # y <- c()
  # for (j in 1:1997) {
  #   y[j] <- abs(x[j+1]+x[j])
  # }
  # u <- c()
  # u <-  qfrechet(U[which(y==y[y<(10^(-12))][1])])
  # if (length(u)==0) {
  #   u <- qfrechet(max(U))
  # }
   u <- (x2+x3)/2*10^(1.2)
  l[i] <- integrate(integral,lower=0,upper=u)$value 
  x4[i] <- optim(par=1,fn=to_opt,lower=0,upper=u,method="Brent")$par
}

# col="#009ADA"

sims1 <- data.frame(Y1=X2,Y2=x4,Y3=X3)%>% apply(c(1,2),FUN=frechet_laplace_pit) %>% as.data.frame()
sims <- sims %>% mutate(out=((X2+X3)<(1)))
sims1 <- sims1 %>% mutate(out=((X2+X3)<(1)))
l <- min(sims,sims1)
u <- max(sims,sims1)
grid.arrange(ggplot(sims) + geom_point(aes(x=Y1,y=Y2),alpha=0.5,size=0.1,col="#009ADA") + xlim(c(l,u))+ ylim(c(l,u)) + xlab(TeX("$Y_1$"))+ ylab(TeX("$Y_2$")),
             ggplot(sims) + geom_point(aes(x=Y2,y=Y3),alpha=0.5,size=0.1,col="#009ADA")+ xlim(c(l,u))+ ylim(c(l,u))+ xlab(TeX("$Y_2$"))+ ylab(TeX("$Y_3$")),
             ggplot(sims) + geom_point(aes(x=Y1,y=Y3),alpha=0.5,size=0.1)+ xlim(c(l,u))+ ylim(c(l,u))+ xlab(TeX("$Y_1$"))+ ylab(TeX("$Y_3$")),
             ggplot(sims1) + geom_point(aes(x=Y1,y=Y2),alpha=0.5,size=0.1,col="#009ADA")+ xlim(c(l,u))+ ylim(c(l,u))+ xlab(TeX("$Y_1$"))+ ylab(TeX("$Y_2^*$")),
             ggplot(sims1) + geom_point(aes(x=Y2,y=Y3),alpha=0.5,size=0.1,col="#009ADA")+ xlim(c(l,u))+ ylim(c(l,u))+ xlab(TeX("$Y_2^*$"))+ ylab(TeX("$Y_3$")),
             ggplot(sims1) + geom_point(aes(x=Y1,y=Y3),alpha=0.5,size=0.1)+ xlim(c(l,u))+ ylim(c(l,u))+ xlab(TeX("$Y_1$"))+ ylab(TeX("$Y_3$")),ncol=3)

x <- c()
U <- seq(from=0.01,to=0.99,length.out=198)
x2 <- X2[1]
x3 <- X3[1]
l <- integrate(integral,lower=0,upper=10^4)$value
for (i in 1:198) {
x[i] <- integral(s=1/qfrechet(U[i]))/l 
}
plot(qfrechet(U),x)
integrate(integral,lower=0,upper=10^4)

# plot several different densities
tmp <- data.frame(z=as.numeric(),cdf=as.numeric())
tmp1 <- data.frame(z=as.numeric(),cdf=as.numeric(),ite=as.character())
set.seed(12)
for (j in 1:20) {
  # x2 <- max(X2)
  # x3 <- X3[X2==max(X2)]
  x2 <- X2[j]
  x3 <- X3[j]
  u <- (x2+x3)/2*10^(1.2)
  l <- integrate(integral,lower=0,upper=10^5)$value
  for (i in 1:1998) {
    x[i] <- integral(s=1/qfrechet(U[i]))/l 
  }
  from <- (j-1)*1998+1
  to <- 1998*j
  tmp[from:to,] <- cbind(data.frame(z=qfrechet(U),cdf=x))
  tmp1[from:to,] <- cbind(tmp[from:to,],data.frame(ite=rep(paste0("x2=",round(x2,2),", x3=",round(x3,2)),1998)))
}
ggplot(tmp1[from:to,]) + geom_line(aes(x=z,y=cdf,col=ite),alpha=1,linewidth=1)+ylab(TeX(paste0("$h($","$x_4$","$)")))+ xlim(c(0,u))+
  xlab(TeX("$x_4$"))

for (i in 1:198) {
  x[i] <- integral(s=1/qfrechet(U[i]))/l 
}
y <- c()
for (i in 1:197) {
 y[i] <- abs(x[i+1]-x[i]) 
}

# generate from 4 variables with logistic links ----

x_y <- as.data.frame(exp(evd::rbvevd(N,dep=1/2,model="log")))
names(x_y) <- c("x1","y1")
x_y <- x_y %>% mutate(x=as.numeric(map(.x=x1,.f=frechet_laplace_pit))) %>% 
  mutate(y=as.numeric(map(.x=y1,.f=frechet_laplace_pit)))
x_y$tf <- (x_y$x>frechet_laplace_pit(qfrechet(0.99)) & x_y$y>frechet_laplace_pit(qfrechet(0.99)))
xy <- as.data.frame(mvrnorm(N,mu=c(0,0),Sigma = matrix(c(1,1/2,1/2,1),ncol=2)))
names(xy) <- c("x","y")
xy$tf <- (xy$x>qnorm(0.99) &  xy$y>qnorm(0.99))
xy <- xy %>% mutate(x=as.numeric(map(.x=qfrechet(pnorm(x)),.f=frechet_laplace_pit))) %>% 
  mutate(y=as.numeric(map(.x=qfrechet(pnorm(y)),.f=frechet_laplace_pit)))
vL <- frechet_laplace_pit(qfrechet(0.99))
p1 <- ggplot(x_y) + geom_point(aes(x=x,y=y,col=tf),size=0.8,alpha=0.5) + xlab("") + ylab("") + scale_color_manual(values = c("black", "#009ADA")) + theme(legend.position="none") + coord_fixed() + xlim(-12,12) + ylim(-12,12)+
               geom_vline(xintercept=vL,color="#009ADA",linetype="dashed") +
              geom_hline(yintercept=vL,color="#009ADA",linetype="dashed") 
                 
p2 <- ggplot(xy) + geom_point(aes(x=x,y=y,col=tf),size=0.8,alpha=0.5) + xlab("") + ylab("")+ scale_color_manual(values = c("FALSE"="black","TRUE" = "#009ADA")) + theme(legend.position="none") + coord_fixed()+ xlim(-12,12) + ylim(-12,12) +
  geom_vline(xintercept=vL,color="#009ADA",linetype="dashed") +
  geom_hline(yintercept=vL,color="#009ADA",linetype="dashed") 
# p3 <- ggplot(x_y %>% filter(x> frechet_laplace_pit(qfrechet(0.99)))%>% filter(y> frechet_laplace_pit(qfrechet(0.99)))) + geom_point(aes(x=x,y=y),alpha=0.5,col="#009ADA") + xlab("") + ylab("")
# p4 <- ggplot(xy %>% filter(x>qnorm(0.99))%>% filter(y>qnorm(0.99))) + geom_point(aes(x=x,y=y),alpha=0.5,col="#009ADA") + xlab("") + ylab("")
p <- grid.arrange(p1,p2,ncol=2)
ggsave(p,filename = "AD_AI_Laplacemargin.png",width=10,height=5)
ggsave(p,filename = "AD_AI_Laplacemargin.pdf",width=10,height=5)

# simulation study for different margin methods ----
set.seed(12)
N <- 50000
v <- 0.99
sims <- generate_Y(N=N) %>% link_log(dep=1/2) %>%
  apply(c(1,2),FUN=frechet_laplace_pit) %>% as.data.frame()
tmp_est <- par_est(sims,v=v,given=c(1),method="Normal")
tmp_est1 <- par_est(sims,v=v,given=c(1),method="AGG")
 
for (i in 1:100) {
  set.seed(12*i)
  N <- 50000
  sims <- generate_Y(N=N) %>% link_log(dep=1/2) %>%
    apply(c(1,2),FUN=frechet_laplace_pit) %>% as.data.frame()
  tmp_est <- rbind(tmp_est,par_est(sims,v=v,given=c(1),method="Normal"))
  tmp_est1 <- rbind(tmp_est1,par_est(sims,v=v,given=c(1),method="AGG"))
}

# calculate root mean square error
mse_an <- (mean((tmp_est$a - 1)^2))
mse_aa <- (mean((tmp_est1$a - 1)^2))
mse_bn <- (mean((tmp_est$b - 0)^2))
mse_ba <- (mean((tmp_est1$b - 0)^2))
bias_an <- mean(tmp_est$a)-1
bias_aa <- mean(tmp_est1$a)-1
bias_bn <- mean(tmp_est$b)-0
bias_ba <- mean(tmp_est1$b)-0
# calculate variance
var_an <- sum((tmp_est$a-mean(tmp_est$a))^2)/nrow(tmp_est)
var_aa <- sum((tmp_est1$a-mean(tmp_est1$a))^2)/nrow(tmp_est1)
var_bn <- sum((tmp_est$b-mean(tmp_est$b))^2)/nrow(tmp_est)
var_ba <- sum((tmp_est1$b-mean(tmp_est1$b))^2)/nrow(tmp_est1)
var_an+bias_an^2
var_aa+bias_aa^2
var_bn+bias_bn^2
var_ba+bias_ba^2

summary(tmp_est)
summary(tmp_est1)
ggpairs(tmp_est,col=2:3)
ggpairs(tmp_est1,col=2:3)
names(tmp_est1) <- paste0(names(tmp_est1),"1")
tmp <- cbind(tmp_est,tmp_est1)
# exploratory plots to try identify source of the rmse
ggplot(tmp[1:100,]) + geom_segment(aes(x=a,y=b,xend=a1,yend=b1)) +
  geom_point(aes(x=a,y=b))+
  geom_point(aes(x=a1,y=b1),col="#C11432") +
  xlab(TeX("$\\hat{\\alpha}$")) +
  ylab(TeX("$\\hat{\\beta}$"))
ggpairs(tmp_est1 %>% filter(a1<0.9 | b1>0.3),col=2:7)

ggplot(tmp[1:100,]) + geom_segment(aes(x=mu,y=sig,xend=mu1,yend=sig1)) +
  geom_point(aes(x=mu,y=sig))+
  geom_point(aes(x=mu1,y=sig1),col="#C11432") +
  xlab(TeX("$\\hat{\\alpha}$")) +
  ylab(TeX("$\\hat{\\beta}$"))

ggplot(tmp[1:100,]) + geom_segment(aes(x=deltal1,y=deltau1,xend=deltal1,yend=deltau1)) +
  geom_point(aes(x=deltal1,y=deltau1))+
  geom_point(aes(x=deltal1,y=deltau1),col="#C11432")

#ggplot(tmp) + geom_point(aes(x=a,y=deltal1))+ geom_point(aes(x=a,y=deltau1),col="red")
  
  ggplot(tmp_est) + geom_boxplot(aes(y=a))

# incomplete gamma function
library(pracma)
# sample from AGG ----
AGG_cdf  <- function(x,theta) {
    sigl <- theta[1]
    sigu <- theta[2]
    deltal <- theta[3]
    deltau <- theta[4]
    C <- (sigl/deltal*gamma(1/deltal) + sigu/deltau*gamma(1/deltau)  )^(-1)
if (x<0) {
  C*sigl/deltal*pracma::gammainc((-x/sigl)^deltal,1/deltal)[2]
}
    else {
      C*sigl/deltal*gamma(1/deltal) + C*sigu/deltau*pracma::gammainc((x/sigu)^deltau,1/deltau)[1]
    }
  }

s <- seq(-3,3,length.out=500)
cdf <- c()
s <- seq(0.01,0.99,length.out=1000)
  cdf <- c()
  for (i in 1:length(s)) {
    x <- s[i]
    cdf[i] <- AGG_cdf(x,theta = c(1,1,1,2))
  }
# generate from CDF
theta <- c(1,1,2,2)
to_opt <- function(x) {
  (AGG_cdf(x=x,theta=theta)-s[i])^2
}
x <- c()
for (i in 1:1000) {
  x[i] <- optim(par=0,fn=to_opt,lower=-10,upper=10,method="Brent")$par
}
plot(density(x)) 
  plot(x)

rAGG <- function(theta) {
  s <- seq(0.0001,0.9999,length.out=10000)
  sigl <- theta[1]
  sigu <- theta[2]
  deltal <- theta[3]
  deltau <- theta[4]
  C <-  (sigl/deltal*gamma(1/deltal) + sigu/deltau*gamma(1/deltau)  )^(-1)
  AGG_sample <- c()
  for (i in 1:length(s)) {
  U <- runif(1)
  to_opt <- function(x) {
  if (x<0) {
   ( C*sigl/deltal*pracma::gammainc((-x/sigl)^deltal,1/deltal)[2] - U)^2
  }
  else {
   ( C*sigl/deltal*gamma(1/deltal) + C*sigu/deltau*pracma::gammainc((x/sigu)^deltau,1/deltau)[1] - U)^2
  }
  }
  AGG_sample[i] <- optim(par=0,fn=to_opt)$par
  }
  return(AGG_sample)
} 

AGG_density <- function(x,theta) {
  mu <- theta[1]
  sigl <- theta[2]
  sigu <- theta[3]
  deltal <- theta[4]
  deltau <- theta[5]
  C_AGG <-  (sigl/deltal*gamma(1/deltal) + sigu/deltau*gamma(1/deltau)  )^(-1)
  y <- c()
  for (i in seq(x)) {
  if (x[i]<mu) {
  y[i] <- C_AGG*exp(-abs((mu-x[i])/sigl)^deltal)
  }
  else {y[i] <- C_AGG*exp(-((x[i]-mu)/sigu)^deltau)}
  }
  return(y)
}

df <- data.frame(matrix(nrow=0,ncol=3))
names(df) <- c("AGG_sample", "param","sim")
thetas <- data.frame("sigl"=c(1/2,1,1/2,1),"sigu"=c(1,1,1,1),
           "deltal"=c(2,2,1,1),"deltau"=c(2,2,2,2))
x <- seq(-5,5,0.0005)
N <- length(x)
for (i in 1:nrow(thetas)) {
  df <-   rbind(df,data.frame("AGG_sample"=AGG_density(x=x,theta=as.numeric(thetas[i,])),
                              "param"=rep(paste0(thetas[i,1],",",thetas[i,2],",",thetas[i,3],",",thetas[i,4]),N),
                              "sim"=rep(i,N)))
}
df <- df %>% mutate(sim=as.factor(sim))
ggplot(df,aes(x=rep(seq(-5,5,0.0005),4),y=AGG_sample,color=sim)) + geom_line(linewidth=1,alpha=0.7) + xlim(c(-5,5)) +
  scale_color_manual(name="Left tail\nparameters",
                     labels=c(TeX(paste0("$\\sigma_l=",thetas[1,1],", \\delta_l=$",thetas[1,3])),
                              TeX(paste0("$\\sigma_l=",thetas[2,1],", \\delta_l=$",thetas[2,3])),
                              TeX(paste0("$\\sigma_l=",thetas[3,1],", \\delta_l=$",thetas[3,3])),
                              TeX(paste0("$\\sigma_l=",thetas[4,1],", \\delta_l=$",thetas[4,3]))),
                     breaks=c(1,2,3,4),
                     values=c("#C11432","#66A64F","#FDD10A","#009ADA")) +
  xlab(TeX("$z$")) + ylab(TeX("$f_{AGG} (z)$"))

# plot limited extrapolation
set.seed(123)
N <- 5000
v <- 0.99
sims <- generate_Y(N=N) %>% link_log(dep=1/2) %>%
  apply(c(1,2),FUN=frechet_laplace_pit) %>% as.data.frame()
pe <- par_est(df = sims, given = 1, v = v)
Y_given1extreme <- sims %>% filter(Y1>quantile(Y1,v))
Y1 <- Y_given1extreme$Y1
Y2 <- Y_given1extreme$Y2
# calculate observed residuals
Z2 <- (Y2-pe$a[1]*Y1)/(Y1^pe$b[1])
# sample from the observed residuals
Nsim <- 5000
# Zstar <- data.frame(Z2 = sample(Z2, size = Nsim, replace = TRUE) # without noise
Zstar <- data.frame(Z2 = sample(Z2, size = Nsim, replace = TRUE)+ rnorm(n=Nsim,mean=0,sd=density(Z2)$bw))
# simulate Y1
Y1_gen <- -log(2*(1-v)) + rexp(Nsim)
Gen_Y1 <- data.frame(Y1=Y1_gen)
Gen_Y1$Y2 <- pe$a[1]*Y1_gen + Y1_gen^pe$b[1] *Zstar[,1]
Gen_Y1$sim <- "model"
Gen_orig <- rbind(Gen_Y1,Y_given1extreme %>% dplyr::select(Y1,Y2) %>% mutate(sim=rep("original_laplace",N*(1-v))))
ggplot(Gen_orig) + geom_point(aes(x=Y1,y=Y2,col=sim),alpha=0.5) +
  scale_color_manual(values = c("original_laplace"="black","model" = "#C11432"))+ 
  xlab(TeX("$Y_1$")) + ylab(TeX("$Y_2$"))

ggplot() +
  geom_point(data=Gen_Y1,mapping = aes(x=Y1,y=Y2),alpha=0.5,col = "#C11432") + coord_fixed() + geom_point(data=Y_given1extreme,mapping=aes(x=Y1,Y2), alpha = 1) 


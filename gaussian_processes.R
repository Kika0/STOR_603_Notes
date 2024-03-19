library(MASS)
library(tidyverse)
library(latex2exp)
library(viridis)
library(plgp)
# set theme defaults to be black rectangle
theme_set(theme_bw())
theme_replace(
  panel.spacing = unit(2, "lines"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA) )

gaussprocess <- function(from = 0, to = 1, K = function(s, t) {exp(-(abs(s-t)/lambda)^alpha)},
                         start = NULL, m = 50,alpha=1,lambda=1) {
  
  t <- seq(from = from, to = to, length.out = m)
  Sigma <- sapply(t, function(s1) {
    sapply(t, function(s2) {
      K(s1, s2)
    })
  })
  
  path <- mvrnorm(mu = rep(0, times = m), Sigma = Sigma)
  # if (!is.null(start)) {
  #   path <- path - path[1] + start  # Must always start at "start"
  # }
  
  return(data.frame("t" = t, "xt" = path))
}

# run for different values of lambda ----
m <- 50
set.seed(1)
l1 <- gaussprocess(lambda=0.1)
set.seed(1)
l2 <- gaussprocess(lambda=1)
set.seed(1)
l3 <- gaussprocess(lambda=10)
df <- cbind(rbind(l1,l2,l3),
            data.frame(lambda=c(rep("0.1",m),rep("1",m),rep("10",m))))
ggplot(df) + geom_line(aes(x=t,y=xt,col=lambda))+ ylab(TeX(paste0("$X($","$s$","$)")))

# run for different values of alpha ----
m <- 50
set.seed(1)
l1 <- gaussprocess(alpha = 0.1,m=m)
set.seed(1)
l2 <- gaussprocess(alpha = 1,m=m)
set.seed(1)
l3 <- gaussprocess(alpha = 1.9,m=m)
df <- cbind(rbind(l1,l2,l3),
            data.frame(alpha=c(rep("0.1",m),rep("1",m),rep("1.9",m))))
ggplot(df) + geom_line(aes(x=t,y=xt,col=alpha))+ ylab(TeX(paste0("$X($","$s$","$)")))

# sample many times to illustrate the Gaussian density at each 1:m ----
tmp <- data.frame(t=as.numeric(),xt=as.numeric())
tmp1 <- data.frame(t=as.numeric(),xt=as.numeric(),ite=as.character())
set.seed(1234)
m <- 10
for (i in 1:1000) {
  from <- (i-1)*m+1
  to <- m*i
  tmp[from:to,] <- gaussprocess(m=m)
  tmp1[from:to,] <- cbind(tmp[from:to,],data.frame(ite=rep(as.character(i),m)))
}

ggplot(tmp1) + geom_line(aes(x=t,y=xt,col=ite),alpha=0.2,linewidth=0.05)+ylab(TeX(paste0("$X($","$s$","$)")))+   theme(legend.position="none") + 
  scale_color_manual(values=c(rep("#C11432",500),rep("#009ADA",500)))
ggplot() + geom_density(tmp1 %>% mutate(t=as.character(t)), mapping=aes(x = xt, col = t),alpha = 0.1,linewidth=0.3)+ theme(legend.position="none") +
  geom_density(data.frame(x=rnorm(1000)),mapping=aes(x=x))

# conditioning on given values use GP regression ----
gaussprocessadd <- function(from = 0, to = 1,df, K = function(s1,s2 ) {exp(-(abs(s1-s2)/lambda)^alpha)},
                         start = NULL, m = 50,alpha=1,lambda=1) {
  t <- df$t
  y <- df$xt
  ts <- seq(from = from, to = to, length.out = m)
  Kxx <- sapply(t, function(s1) {
    sapply(t, function(s2) {
      K(s1, s2)
    })
  })
  
  Kxxs <- sapply(ts, function(s1) {
    sapply(t, function(s2) {
      K(s1, s2)
    })
  })
  
  Kxsx <- sapply(t, function(s1) {
    sapply(ts, function(s2) {
      K(s1, s2)
    })
  })
  
  Kxsxs <- sapply(ts, function(s1) {
    sapply(ts, function(s2) {
      K(s1, s2)
    })
  })
  
  mu <- Kxsx%*%solve(Kxx)%*%matrix(y,ncol=1)
  Sigma <- Kxsxs-Kxsx%*%solve(Kxx)%*%Kxxs
  path <- mvrnorm(mu = mu, Sigma = Sigma)
  # if (!is.null(start)) {
  #   path <- path - path[1] + start  # Must always start at "start"
  # }
  
  return(data.frame("t" = ts, "xt" = path))
}
set.seed(123)
df <- gaussprocess(m=11)
tmp <- gaussprocessadd(df=df,m=51)
tmp1 <- gaussprocessadd(df=df,m=51)
tmp2 <- gaussprocessadd(df=df,m=51)
ggplot() + geom_line(data=tmp,aes(x=t,y=xt),col="#C11432") +geom_point(data=df,aes(x=t,y=xt),size=2) + geom_line(data=tmp1,aes(x=t,y=xt),col="#009ADA")+ geom_line(data=tmp2,aes(x=t,y=xt),col="#66A64F")

# sample in two dimensions ----
x <- seq(0,1)
gaussprocess2d <- function(from = 0, to = 1, K = function(s, t) {exp(-(sqrt((s[1]-t[1])^2+(s[2]-t[2])^2)/lambda)^alpha)},
                         start = NULL, rho=NULL, sig=NULL, m = 1000,alpha=1,lambda=1) {
  
  x <- seq(from = from, to = to, length.out = m)
  y <- seq(from = from, to = to, length.out = m)
  xy <- expand.grid(x,y) 
  xy.list <- split(xy, seq(nrow(xy)))
  Sigma <- matrix(ncol=m^2,nrow=m^2)
for (i in 1:nrow(xy)){
    for (j in 1:nrow(xy)) {
    Sigma[i,j] <-  as.numeric(K(s= as.numeric(xy[i,]), t=as.numeric(xy[j,])) )
    }
  }
  
  path <- mvrnorm(mu = rep(0, times = m^2), Sigma = Sigma)
  # if (!is.null(start)) {
  #   path <- path - path[1] + start  # Must always start at "start"
  # }
  
  return(data.frame("x" = xy$Var1, "y"=xy$Var2, "xt" = path))
}

K = function(s, t) {exp(-(mahalanobis(x=s,center = t,
                    cov = matrix(c(1,rho*sig,rho*sig,sig^2),ncol=2))/lambda)^alpha)}

set.seed(1)
tmp <- gaussprocess2d(m=10,K = function(s, t) {exp(-(mahalanobis(x=s,center = t,
                                                                 cov = matrix(c(1,rho*sig,rho*sig,sig^2),ncol=2))/lambda)^alpha)},
              rho=1/3,sig=1/2        )
ggplot(tmp,aes(x=x,y=y)) +
  geom_raster(aes(fill=xt), interpolate = TRUE) +
  geom_contour(aes(z=xt), bins = 12, color = "gray30", 
               size = 0.5, alpha = 0.5) +
  coord_equal() +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_viridis_c(option = "viridis")
# try mahalanobis distance for a range of parameters
# number of samples
mx <- 20
x <- seq(0,1,length=mx)
# grid of pairwise values
X <- expand.grid(x, x)
rho <- c(0.1,0.5,0.9)
rho <- c(0,0.1,0.2)
sig <- c(0.5,1,2)
pp <- data.frame(y=as.numeric(),x1=as.numeric(),x2=as.numeric(),ite=as.character())

for (i in 1:length(rho)) {
  for (j in 1:length(sig)) {
    # sample from multivariate normal with mean zero, sigma = sigma
    set.seed(1)
    Y <- gaussprocess2d(m=mx,rho = rho[i],sig=sig[j],K = function(s, t) {exp(-(mahalanobis(x=s,center = t,
                                                                                           cov = matrix(c(1,rho[i]*sig[j],rho[i]*sig[j],sig[j]^2),ncol=2))/lambda)^alpha)})$xt
    pp <- rbind(pp,data.frame(y=Y,x1=X[,1],x2=X[,2],ite=paste0("rho=",rho[i],", sigma=",sig[j])))
  }
}

# kernel function
rbf_D <- function(X,lambda=1, alpha=1, eps = sqrt(.Machine$double.eps) ){
  D <- plgp::distance(X)
  Sigma <- exp(-(D/lambda)^alpha) + diag(eps, nrow(X))
}
# number of samples
mx <- 30
x <- seq(0,1,length=mx)
# grid of pairwise values
X <- expand.grid(x, x)
alpha <- c(0.5,1,1.5)
lambda <- c(0.5,1,2)
pp <- data.frame(y=as.numeric(),x1=as.numeric(),x2=as.numeric(),ite=as.character())
for (i in 1:length(alpha)) {
  for (j in 1:length(lambda)) {
    # sample from multivariate normal with mean zero, sigma = sigma
    Y <- gaussprocess2d(m=mx,alpha = alpha[i],lambda=lambda[j])$xt
    pp <- rbind(pp,data.frame(y=Y,x1=X[,1],x2=X[,2],ite=paste0("lambda=",lambda[j],", alpha=",alpha[i])))
  }
} 

ggplot(pp,aes(x=x1,y=x2)) +
  geom_raster(aes(fill=y), interpolate = TRUE) +
  geom_contour(aes(z=y), bins = 12, color = "gray30", 
               size = 0.5, alpha = 0.5) +
  coord_equal() +
  facet_wrap(~ite) +
  # ylab(TeX("$\\Phi(\\boldsymbol(x))$"))+
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_viridis_c(option = "viridis")

# compute squared exponential kernel on pairwise values
Sigma <- rbf_D(X,lambda = 0.01,alpha = 1)

# sample from multivariate normal with mean zero, sigma = sigma
Y <- MASS::mvrnorm(1,rep(0,dim(Sigma)[1]), Sigma)

# plot results
pp <- data.frame(y=Y,x1=X[,1],x2=X[,2])

# combine GP with conditional models----
m <- 10
n <- 1000

# sample many times to illustrate the Gaussian density at each 1:m ----
tmp <- data.frame(t=as.numeric(),xt=as.numeric())
tmp1 <- data.frame(t=as.numeric(),xt=as.numeric(),ite=as.character())
set.seed(1234)
for (i in 1:n) {
  from <- (i-1)*m+1
  to <- m*i
  tmp[from:to,] <- gaussprocess(m=m)
  tmp1[from:to,] <- cbind(tmp[from:to,],data.frame(ite=rep(as.character(i),m)))
}

ggplot(tmp1) + geom_line(aes(x=t,y=xt,col=ite),alpha=0.5,linewidth=0.1)+ylab(TeX(paste0("$X($","$s$","$)")))+   theme(legend.position="none") 
 # scale_color_manual(values=c(rep("#C11432",500),rep("#009ADA",500)))
ggplot() + geom_density(tmp1 %>% mutate(t=as.character(t)), mapping=aes(x = xt, col = t),alpha = 0.1,linewidth=0.3)+ theme(legend.position="none") +
  geom_density(data.frame(x=rnorm(n)),mapping=aes(x=x))

paste0("Y",1:10)
# transform data to a correct format
tmp <- tmp1 %>% select(xt,t) %>%   mutate(id = row_number(), .by =t) %>%
  pivot_wider(names_from = t, values_from = xt, id_cols = id) %>% select(-id)
colnames(tmp) <- paste0("Y",1:(m+1))
# PIT to Laplace
sims <- tmp %>% mutate(Y_1=as.numeric(map(.x=Y1,.f=frechet_laplace_pit))) %>% 
  mutate(Y_2=as.numeric(map(.x=Y2,.f=frechet_laplace_pit))) %>%
  mutate(Y_3=as.numeric(map(.x=Y3,.f=frechet_laplace_pit)))

# check it looks Laplace 
ggplot(sims %>% dplyr::select(Y_1,Y_2,Y_3) %>% pivot_longer(everything())) + geom_density(aes(x=value),stat="density") + facet_wrap(~name)
grid.arrange(ggplot(sims) + geom_point(aes(x=Y_1,y=Y_2),alpha=0.5),
             ggplot(sims) + geom_point(aes(x=Y_2,y=Y_3),alpha=0.5),
             ggplot(sims) + geom_point(aes(x=Y_1,y=Y_3),alpha=0.5),ncol=3)



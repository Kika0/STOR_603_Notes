library(tidyverse)
set.seed(123456789)
library(MASS)
library(gridExtra)
library(latex2exp)
library(mev)
library(evd)
library(xtable)

# use this to block out the upper part
upper_corner = tibble(x = c(0, 1, 1),
                       y = c(1, 1, 0))

# dependence between r.v. for husler reiss dist ----
dep_x1_x2 = 1
dep_x2_x3 = 1

# dependence matrix 
sigma=matrix(c(0, dep_x2_x3, dep_x1_x2,
               dep_x2_x3, 0, dep_x2_x3,
               dep_x1_x2, dep_x2_x3, 0), 
             nrow = 3, byrow = T)
r <- 0.5
#sigma=matrix(c(1,r,r,r,1,r,r,r,1),nrow=3,byrow = T)

# simulate from extreme MVD
sims_low_dependence = mev::rmev(n=5000, sigma = sigma, d = 3, model = "hr")
#sims_low_dependence = mvrnorm(n=5000,mu=c(0,0,0),Sigma=sigma,empirical = TRUE)

# pseudo polar decomposition to find angular components
dat = tibble(x1 = sims_low_dependence[,1], x2 = sims_low_dependence[,2], x3 = sims_low_dependence[,3])%>%
  mutate(R = x1 + x2 + x3)  %>% 
  mutate(u = quantile(R, 0.9)) %>%
  filter(R>u) %>%
  mutate(w1 = x1/R, 
         w2 = x2/R,
         w3 = x3/R) 



# use this to block out the upper
upper_courner = tibble(x = c(0, 1, 1),
                       y = c(1, 1, 0))


dat %>%
  ggplot()+ 
  geom_density_2d_filled(aes(w1,w2,)) + 
  geom_polygon(data = upper_courner, aes(x,y),col = 'white', fill = "white")+
  theme_minimal()+
  scale_y_continuous(limits = c(0, 1), breaks = c())+
  scale_x_continuous(limits = c(0, 1), breaks = c())+
  theme(panel.grid.major = element_blank(), 
        legend.position = 'none',
        panel.grid.minor = element_blank())+
  labs(x = "", y = "")

angles <- function(rho,u) {
  r <- rho
  sigma=matrix(c(1,r,r,r,1,r,r,r,1),nrow=3,byrow = T)
  
  # simulate from normal
  sims_low_dependence = mvrnorm(n=500000,mu=c(0,0,0),Sigma=sigma,empirical=TRUE)
  
#   X <- dat$X1
#   Y <- dat$X2
#   X <- rnorm(n=1000,mean=20)
#   Y <- rnorm(n=1000,mean=20)
# data.frame(X=X,Y=Y,R=X+Y,W=X/(X+Y)) %>% view()
#   X <- x
#   Y <- y
#   X <- sims_low_dependence[,2]
#   Y <- sims_low_dependence[,3]
#   R <- X+Y
#  R <-  R
#   W <- X/R
#  plot(W,R)
  
  # pseudo polar decomposition to find angular components
  df <-  data.frame(sims_low_dependence) %>% rowid_to_column()
  #df <- df[(df$X1+df$X2+df$X3)>6,]
  df <- data.frame(X1=rnorm(100),X2=rnorm(100),X3=rnorm(100)) %>% rowid_to_column()
  
  # create x and y quantile variables
  U1 <- df %>% dplyr::select(rowid,X1) %>% arrange(X1) %>% mutate(U1=row_number()/(nrow(df)+1)) %>% dplyr::select(-X1)
  U2 <- df %>% dplyr::select(rowid,X2) %>% arrange(X2) %>% mutate(U2=row_number()/(nrow(df)+1))%>% dplyr::select(-X2)
  U3 <- df %>% dplyr::select(rowid,X3) %>% arrange(X3) %>% mutate(U3=row_number()/(nrow(df)+1))%>% dplyr::select(-X3)
  dat <- df %>% dplyr::select(rowid) %>%  left_join(U1,by="rowid") %>% left_join(U2,by="rowid") %>% left_join(U3,by="rowid")
  dat <- dat %>% mutate(G1=-log(-log(U1))) %>% mutate(G2=-log(-log(U2))) %>% mutate(G3=-log(-log(U3)))
  plot(dat$G1,dat$G2)
  
  dat <- dat %>%
    mutate(R = G1 + G2 + G3)  %>%
    mutate(u = quantile(R, u)) %>%
    filter(R>u) %>%
    mutate(w1 = G1/R,
           w2 = G2/R,
           w3 = G3/R)

  # use this to block out the upper
  upper_courner = tibble(x = c(0, 1, 1),
                         y = c(1, 1, 0))
  dat %>%
    ggplot()+ 
    geom_density_2d_filled(aes(w1,w2)) + 
    geom_polygon(data = upper_courner, aes(x,y),col = 'white', fill = "white")+
    theme_minimal()+
    scale_y_continuous(limits = c(0, 1), breaks = c())+
    scale_x_continuous(limits = c(0, 1), breaks = c(seq(0,200,20)))+
   # theme(panel.grid.major = element_blank(),
   #       legend.position = 'none',
   #      panel.grid.minor = element_blank())+
    labs(x = "", y = "") +ggtitle(TeX(paste0("$\\rho=.",r*10,",u=$",u,"")))
}
p1 <- angles(rho=0,u=0.9)
p2 <- angles(rho=0,u=0.9)
p3 <- angles(rho=0,u=0.9999)




p1 <- angles(rho=0.1,u=0.9)
p4 <- angles(rho=0.5,u=0.9)
p7 <- angles(rho=0.9,u=0.9)
p2 <- angles(rho=0.1,u=0.99)
p5 <- angles(rho=0.5,u=0.99)
p8 <- angles(rho=0.9,u=0.99)
p3 <- angles(rho=0.1,u=0.999)
p6 <- angles(rho=0.5,u=0.999)
p9 <- angles(rho=0.9,u=0.999)
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,ncol=3)

# dependence between r.v. for husler reiss dist
h <- function(l,u,N=500000){
  dep_x1_x2 = l
  dep_x2_x3 = l
  
  # dependence matrix 
  sigma=matrix(c(0, dep_x2_x3, dep_x1_x2,
                 dep_x2_x3, 0, dep_x2_x3,
                 dep_x1_x2, dep_x2_x3, 0), 
               nrow = 3, byrow = T)
  
  #sigma=matrix(c(1,r,r,r,1,r,r,r,1),nrow=3,byrow = T)
  
  # simulate from extreme MVD
  sims_low_dependence = mev::rmev(n=N, sigma = sigma, d = 3, model = "hr")
  #sims_low_dependence = mvrnorm(n=5000,mu=c(0,0,0),Sigma=sigma)
  
  # pseudo polar decomposition to find angular components
  dat = tibble(x1 = sims_low_dependence[,1], x2 = sims_low_dependence[,2], x3 = sims_low_dependence[,3])%>%
    mutate(R = x1 + x2 + x3)  %>% 
    mutate(u = quantile(R, u)) %>%
    filter(R>u) %>%
    mutate(w1 = x1/R, 
           w2 = x2/R,
           w3 = x3/R) 
  
  
  
  # use this to block out the upper
  upper_courner = tibble(x = c(0, 1, 1),
                         y = c(1, 1, 0))
  

  
  dat %>%
    ggplot()+ 
    geom_density_2d_filled(aes(w1,w2),bins=7) + 
    geom_polygon(data = upper_courner, aes(x,y),col = 'white', fill = "white")+
    theme_minimal()+
    scale_y_continuous(limits = c(0, 1))+
    scale_x_continuous(limits = c(0, 1))+
    theme(panel.grid.major = element_blank(), 
         # legend.position = 'none',
          panel.grid.minor = element_blank())+
    labs(x = "", y = "") + ggtitle(TeX(paste0("$\\lambda=$",l,", $u=\\hat{F}_R^{-1}($",u,"$)$"))) +
    guides(fill=guide_legend(title="Density estimate"))
  
}

grid.arrange(h(0.01,0.9),h(0.01,0.99),h(0.01,0.999),
             h(0.25,0.9),h(0.25,0.99),h(0.25,0.999),
             h(10,0.9),h(10,0.99),h(10,0.999),
             ncol=3)

# generate with coefficients with limiting point distributions ---
# adapt function form GEV.qmd
# try another
a1 <- c(17,3,0,1,4)
a <- a1/sum(a1)
b1 <- c(0,11,6,2,1)
b <- b1/sum(b1)
c1 <- c(0,0,3,3,4)
c <- c1/sum(c1)
abc <- data.frame(a,b,c)

abcw1w2 <- abc %>%    mutate(R = a+b+c)  %>% 
  mutate(w1 = a/R, 
         w2 = b/R)


generate_dependent_X_Y_Z <- function(N,abc=abc) {
  set.seed(12)
  d <- 5
  U <- runif(d*N)
  # a <- c(2/3,1/12,0,1/12,1/6)
  # b <- c(0,1/3,1/3,1/6,1/6)
  # c <- c(0,0,1/3,1/3,1/3)
a <- abc[,1]
b <- abc[,2]
c <- abc[,3]

  # generate Y
  Y <- c()
  Y <- -1/(log(U) ) 
  # generate X
  X_1 <- c()
  X_2 <- c()
  X_3 <- c()
  for (j in 1:N) {
    X_1[j] <- max(Y[(d*(j-1)+1):(d*j)]*a)/N
    X_2[j] <- max(Y[(d*(j-1)+1):(d*j)]*b)/N
    X_3[j] <- max(Y[(d*(j-1)+1):(d*j)]*c)/N
  }
  return(data.frame(X_1,X_2,X_3))
}

plot_clusters <- function(sims,u=0.9) {
  sims_low_dependence <- sims
  dat = tibble(x1 = sims_low_dependence[,1], x2 = sims_low_dependence[,2], x3 = sims_low_dependence[,3])%>%
    mutate(R = x1 + x2 + x3)  %>% 
    mutate(u = quantile(R, u)) %>%
    filter(R>u) %>%
    mutate(w1 = x1/R, 
           w2 = x2/R,
           w3 = x3/R) 
  
  dat %>%
    ggplot()+ 
    geom_density_2d_filled(aes(w1,w2),bins=7) + 
    geom_polygon(data = upper_corner, aes(x,y),col = 'white', fill = "white")+
    theme_minimal()+
    scale_y_continuous(limits = c(0, 1))+
    scale_x_continuous(limits = c(0, 1))+
    theme(panel.grid.major = element_blank(), 
          # legend.position = 'none',
          panel.grid.minor = element_blank())+
    labs(x = "", y = "") + ggtitle(TeX(paste0("$u=\\hat{F}_R^{-1}($",u,"$)$"))) +
    guides(fill=guide_legend(title="Density estimate")) 
}

# plot for different thresholds
p1 <- generate_dependent_X_Y_Z(N=50000,abc=abc) %>% plot_clusters(u=0.9)
p2 <- generate_dependent_X_Y_Z(N=50000,abc=abc) %>% plot_clusters(u=0.99)
p3 <- generate_dependent_X_Y_Z(N=50000,abc=abc) %>% plot_clusters(u=0.999)
grid.arrange(p1,p2,p3,ncol=3)

# keep common scale across plots
generate_dep_X_Y_Z <- function(N,abc=abc,U=c(0.9,0.99,0.999)) {
  set.seed(12)
  d <- 5
  unif <- runif(d*N)
  # a <- c(2/3,1/12,0,1/12,1/6)
  # b <- c(0,1/3,1/3,1/6,1/6)
  # c <- c(0,0,1/3,1/3,1/3)
  a <- abc[,1]
  b <- abc[,2]
  c <- abc[,3]
  
  # generate Y
  Y <- c()
  Y <- -1/(log(unif) ) 
  # generate X
  X_1 <- c()
  X_2 <- c()
  X_3 <- c()
  for (j in 1:N) {
    X_1[j] <- max(Y[(d*(j-1)+1):(d*j)]*a)/N
    X_2[j] <- max(Y[(d*(j-1)+1):(d*j)]*b)/N
    X_3[j] <- max(Y[(d*(j-1)+1):(d*j)]*c)/N
  }
 sims <- data.frame(X_1,X_2,X_3)

  sims_low_dependence <- sims
  
  dat <- tibble(x1=numeric(),x2=numeric(),x3=numeric(),R=numeric(),u=numeric(),q=numeric(),w1=numeric(),w2=numeric(),w3=numeric())
  for (i in 1:length(U)) {
    dat1 <- tibble(x1=numeric(),x2=numeric(),x3=numeric(),R=numeric(),u=numeric(),q=numeric(),w1=numeric(),w2=numeric(),w3=numeric())
  dat1 = tibble(x1 = sims_low_dependence[,1], x2 = sims_low_dependence[,2], x3 = sims_low_dependence[,3])%>%
    mutate(R = x1 + x2 + x3)  %>% 
    mutate(u = quantile(R, U[i])) %>%
    mutate(q=rep(U[i],N)) %>% 
    filter(R>u) %>%
    mutate(w1 = x1/R, 
           w2 = x2/R,
           w3 = x3/R) 
 dat <-  rbind(dat,dat1)
  
}
  
p <-   dat %>%
    ggplot()+ 
    geom_density_2d_filled(aes(w1,w2),bins=7) + 
    geom_polygon(data = upper_corner, aes(x,y),col = 'white', fill = "white")+
    facet_wrap(~q) +
    theme_minimal()+
    scale_y_continuous(limits = c(0, 1))+
    scale_x_continuous(limits = c(0, 1))+
    theme(panel.grid.major = element_blank(), 
          # legend.position = 'none',
          panel.grid.minor = element_blank())+
    labs(x = "", y = "") +
  #ggtitle(TeX(paste0("$u=\\hat{F}_R^{-1}($",U,"$)$"))) +
    guides(fill=guide_legend(title="Density estimate")) 
return(p)
}
generate_dep_X_Y_Z(N=50000,abc=abc,U=c(0.9,0.99,0.999))


# create density along edges instead of point mass ----
generate_dependent_X_Y_Z <- function(N,abc=abc,dep) {
  set.seed(12)
  d <- 5
  U <- runif(d*N)
  # a <- c(2/3,1/12,0,1/12,1/6)
  # b <- c(0,1/3,1/3,1/6,1/6)
  # c <- c(0,0,1/3,1/3,1/3)
  a <- abc[,1]
  b <- abc[,2]
  c <- abc[,3]
  
  # generate Y
  # asy<-list(.4,.1,.6,c(.3,.2),c(.1,.1),c(.4,.1),c(.2,.3,.2))
  # Y <- rmvevd(n=N*d,dep=c(0.6,0.5,0.2,0.9),d=3,asy=asy,model="alog",mar=c(1,1,1)) 
 Y <-  exp(rmvevd(n=N*d,dep=dep,model="log",d=3))
  # generate X
  X_1 <- c()
  X_2 <- c()
  X_3 <- c()
  for (j in 1:N) {
    X_1[j] <- max(Y[(d*(j-1)+1):(d*j),1]*a)
    X_2[j] <- max(Y[(d*(j-1)+1):(d*j),2]*b)
    X_3[j] <- max(Y[(d*(j-1)+1):(d*j),3]*c)
  }
  return(data.frame(X_1,X_2,X_3))
}

plot_clusters <- function(sims,u=0.9,dep) {
  sims_low_dependence <- sims
  dat = tibble(x1 = sims_low_dependence[,1], x2 = sims_low_dependence[,2], x3 = sims_low_dependence[,3])%>%
    mutate(R = x1 + x2 + x3)  %>% 
    mutate(u = quantile(R, u)) %>%
    filter(R>u) %>%
    mutate(w1 = x1/R, 
           w2 = x2/R,
           w3 = x3/R) 
  
  dat %>%
    ggplot()+ 
    geom_density_2d_filled(aes(w1,w2),bins=7) + 
    geom_polygon(data = upper_corner, aes(x,y),col = 'white', fill = "white")+
    theme_minimal()+
    scale_y_continuous(limits = c(0, 1))+
    scale_x_continuous(limits = c(0, 1))+
    theme(panel.grid.major = element_blank(), 
          # legend.position = 'none',
          panel.grid.minor = element_blank())+
    labs(x = "", y = "") + ggtitle(TeX(paste0("$u=\\hat{F}_R^{-1}($",u,"$)$",","," $\\alpha=$",dep))) +
    guides(fill=guide_legend(title="Density estimate")) 
}

# plot for different thresholds
p1 <- generate_dependent_X_Y_Z(N=50000,abc=abc,dep=0.99) %>% plot_clusters(u=0.9,dep=0.99)
p2 <- generate_dependent_X_Y_Z(N=50000,abc=abc,dep=0.99) %>% plot_clusters(u=0.99,dep=0.99)
p3 <- generate_dependent_X_Y_Z(N=50000,abc=abc,dep=0.99) %>% plot_clusters(u=0.999,dep=0.99)
p4 <- generate_dependent_X_Y_Z(N=50000,abc=abc,dep=1/2) %>% plot_clusters(u=0.9,dep=1/2)
p5 <- generate_dependent_X_Y_Z(N=50000,abc=abc,dep=1/2) %>% plot_clusters(u=0.99,dep=1/2)
p6 <- generate_dependent_X_Y_Z(N=50000,abc=abc,dep=1/2) %>% plot_clusters(u=0.999,dep=1/2)
p7 <- generate_dependent_X_Y_Z(N=50000,abc=abc,dep=0.1) %>% plot_clusters(u=0.9,dep=0.1)
p8 <- generate_dependent_X_Y_Z(N=50000,abc=abc,dep=0.1) %>% plot_clusters(u=0.99,dep=0.1)
p9 <- generate_dependent_X_Y_Z(N=50000,abc=abc,dep=0.1) %>% plot_clusters(u=0.999,dep=0.1)
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,ncol=3)

# make density common scale using facet_wrap ----
generate_dep_X_Y_Y_Z <- function(N,abc=abc,dep=1/2,U=c(0.9,0.99,0.999)) {
  set.seed(12)
  d <- 5
  unif <- runif(d*N)
  # a <- c(2/3,1/12,0,1/12,1/6)
  # b <- c(0,1/3,1/3,1/6,1/6)
  # c <- c(0,0,1/3,1/3,1/3)
  a <- abc[,1]
  b <- abc[,2]
  c <- abc[,3]
  dat <- tibble(x1=numeric(),x2=numeric(),x3=numeric(),R=numeric(),u=numeric(),q=numeric(),w1=numeric(),w2=numeric(),w3=numeric())
  
  for (l in 1:length(dep)) {
  # generate Y
  Y <- c()
  Y <-  exp(rmvevd(n=N*d,dep=dep,model="log",d=3))
  # generate X
  X_1 <- c()
  X_2 <- c()
  X_3 <- c()
  for (j in 1:N) {
    X_1[j] <- max(Y[(d*(j-1)+1):(d*j),1]*a)/N
    X_2[j] <- max(Y[(d*(j-1)+1):(d*j),2]*b)/N
    X_3[j] <- max(Y[(d*(j-1)+1):(d*j),3]*c)/N
  }
sims <- data.frame(X_1=X_1,X_2=X_2,X_3=X_3)
  sims %>% head()
  sims_low_dependence <- sims
  
  for (i in 1:length(U)) {
    dat1 <- tibble(x1=numeric(),x2=numeric(),x3=numeric(),R=numeric(),u=numeric(),q=numeric(),w1=numeric(),w2=numeric(),w3=numeric())
    dat1 = tibble(x1 = sims_low_dependence[,1], x2 = sims_low_dependence[,2], x3 = sims_low_dependence[,3])%>%
      mutate(R = x1 + x2 + x3)  %>% 
      mutate(u = quantile(R, U[i])) %>%
      mutate(q=rep(U[i],N)) %>% 
      filter(R>u) %>%
      mutate(w1 = x1/R, 
             w2 = x2/R,
             w3 = x3/R) 
    dat <-  rbind(dat,dat1)
    
  }
  
  }
  p <-   dat %>%
    ggplot()+ 
    geom_density_2d_filled(aes(w1,w2),bins=7) + 
    geom_polygon(data = upper_corner, aes(x,y),col = 'white', fill = "white")+
    facet_wrap(~q) +
    theme_minimal()+
    scale_y_continuous(limits = c(0, 1))+
    scale_x_continuous(limits = c(0, 1))+
    theme(panel.grid.major = element_blank(), 
          # legend.position = 'none',
          panel.grid.minor = element_blank())+
    labs(x = "", y = "") +
    #ggtitle(TeX(paste0("$u=\\hat{F}_R^{-1}($",U,"$)$"))) +
    guides(fill=guide_legend(title="Density estimate")) 
  return(p)
}
generate_dep_X_Y_Y_Z(N=50000,abc=abc,U=c(0.9,0.99,0.999))


# try asymmetric case ----
# generate x and y
a <- 1/2
x_y <- evd::rbvevd(100,dep=a,model="log")
x <- exp(x_y[,1])/N
Y <- exp(x_y[,2])/N
# generate z
to_opt <- function(z) {
  (  (  y^(-(1/a)+1)*(y^(-1/a)+z^(-1/a))^(a-1)*exp(-(y^(-1/a)+z^(-1/a))^a)*exp(1/y)  )-U)^2
}
z <- c()
for (i in 1:nrow(x_y)){
# generate U
U <- runif(1)
y <- Y[i]
# F_Y_Z <- function(z) {
#   -a *y^(-(1/a)+1)*(y^(-1/a)+z^(-1/a))^(a-1)*exp(-(y^(-1/a)+z^(-1/a))^a)
# }
z[i] <- optim(par=1,fn=to_opt)$par
}

plot(Y,z)
plot(x,Y)
plot(x,z)
df <- data.frame(x,Y,z)
p1 <- ggplot(df) + geom_point(aes(x,Y))
p2 <- ggplot(df) + geom_point(aes(z,Y))
p3 <- ggplot(df) + geom_point(aes(x,z))
grid.arrange(p1,p2,p3,ncol=3)

# plot the triangles



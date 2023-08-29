library(tidyverse)
set.seed(123456789)
library(MASS)
library(gridExtra)


# dependence between r.v. for husler reiss dist
dep_x1_x2 = 1
dep_x2_x3 = 1

# dependence matrix 
sigma=matrix(c(0, dep_x2_x3, dep_x1_x2,
               dep_x2_x3, 0, dep_x2_x3,
               dep_x1_x2, dep_x2_x3, 0), 
             nrow = 3, byrow = T)
r <- 0.5
sigma=matrix(c(1,r,r,r,1,r,r,r,1),nrow=3,byrow = T)

# simulate from extreme MVD
sims_low_dependence = mev::rmev(n=5000, sigma = sigma, d = 3, model = "hr")
sims_low_dependence = mvrnorm(n=5000,mu=c(0,0,0),Sigma=sigma)

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
  geom_density_2d_filled(aes(w1,w2)) + 
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
  sims_low_dependence = mvrnorm(n=50000,mu=c(0,0,0),Sigma=sigma)
  
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
    geom_density_2d_filled(aes(w1,w2)) + 
    geom_polygon(data = upper_courner, aes(x,y),col = 'white', fill = "white")+
    theme_minimal()+
    scale_y_continuous(limits = c(0, 1), breaks = c())+
    scale_x_continuous(limits = c(0, 1), breaks = c())+
    theme(panel.grid.major = element_blank(), 
          legend.position = 'none',
          panel.grid.minor = element_blank())+
    labs(x = "", y = "") +ggtitle(TeX(paste0("$\\rho=.",r*10,",u=$",u,"")))
}
p1 <- angles(rho=0,u=0.9)
p2 <- angles(rho=0,u=0.9)
p3 <- angles(rho=0,u=0.99)


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
h <- function(l){
  dep_x1_x2 = l
  dep_x2_x3 = l
  
  # dependence matrix 
  sigma=matrix(c(0, dep_x2_x3, dep_x1_x2,
                 dep_x2_x3, 0, dep_x2_x3,
                 dep_x1_x2, dep_x2_x3, 0), 
               nrow = 3, byrow = T)
  
  #sigma=matrix(c(1,r,r,r,1,r,r,r,1),nrow=3,byrow = T)
  
  # simulate from extreme MVD
  sims_low_dependence = mev::rmev(n=500000, sigma = sigma, d = 3, model = "hr")
  #sims_low_dependence = mvrnorm(n=5000,mu=c(0,0,0),Sigma=sigma)
  
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
    geom_density_2d_filled(aes(w1,w2)) + 
    geom_polygon(data = upper_courner, aes(x,y),col = 'white', fill = "white")+
    theme_minimal()+
    scale_y_continuous(limits = c(0, 1), breaks = c())+
    scale_x_continuous(limits = c(0, 1), breaks = c())+
    theme(panel.grid.major = element_blank(), 
          legend.position = 'none',
          panel.grid.minor = element_blank())+
    labs(x = "", y = "") + ggtitle(TeX(paste0("$\\lambda=$",l)))
  
}

grid.arrange(h(0.01),h(1/4),h(10),ncol=3)

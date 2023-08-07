library(tidyverse)
library(gridExtra)
library(evd)
library(latex2exp)
library(gridExtra)
uk_winter <- readRDS("data/uk_1999_2018_autumn.RDS")
# calculate dependence between X(London) and Y(some other location)
X <- uk_winter[uk_winter$is_location=="london",6:ncol(uk_winter)] %>% as_vector()
Y <- uk_winter[1,6:ncol(uk_winter)] %>% as_vector()
df <- data.frame(X=X,Y=Y) %>% remove_rownames()
# use PIT
# create x and y quantile variables
U <- df %>% select(X) %>% arrange(X) %>% mutate(u=row_number()/(nrow(df)+1))
V <- df %>% select(Y) %>% arrange(Y) %>% mutate(v=row_number()/(nrow(df)+1)) 
# examine top and bottom
U %>% head()
U %>% tail()
V %>% head()
V %>% tail()
df <- df %>% left_join(U,by="X") %>% left_join(V,by="Y")
df %>% head()

# plot
p1 <- ggplot(df) + geom_point(aes(x=X,y=Y),size=0.1) + xlab("X (London temperature)") +
  theme_bw() + ggtitle("Temperature time series") +
  theme(panel.spacing = unit(2, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
p2 <- ggplot(df) + geom_point(aes(x=u,y=v),size=0.1) + xlab("u") +
  theme_bw() + coord_fixed() + ggtitle("Uniform transform") +
  theme(panel.spacing = unit(2, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
grid.arrange(p1,p2,ncol=2)

# calculate the dependence between X and Y
threshold <- 0.95
2-log( mean(df$u < threshold & df$v < threshold))/log(threshold)
# create empirical chi(u) function
chi <- function(df,threshold) {
  2-log( mean(df$u < threshold & df$v < threshold))/log(mean(df$u < threshold))  
}
chi(df,threshold=0.95)

threshold <- sort(df$u)[5:(length(df$u)-5)]
threshold <- seq(0.1,0.995,length.out=1000)  
chi_u <- c()
var_chi_u <- c()
for (i in 1:length(threshold)) {
  chi_u[i] <- chi(df,threshold[i])
  cu <- mean(df$u < threshold[i] & df$v < threshold[i])
  cnst <- qnorm((1 + 0.95)/2)
  varchi <- ((1/log(threshold[i])^2 * 1)/cu^2 * cu * (1 - cu))/length(df$u)
  varchi <- cnst * sqrt(varchi)
var_chi_u[i] <- varchi
 
}

df_chi <- data.frame(u=threshold,chi_u=chi_u,
        CI_upper=chi_u+var_chi_u,CI_lower=chi_u-var_chi_u)
ggplot(df_chi) + ylim(c(0,1)) + geom_line(aes(x=threshold,y=chi_u))+ 
  geom_line(aes(x=threshold,y=CI_upper),lty=2,col="#C11432") +
  ylab(TeX("$\\chi(u)$")) + xlab("u") +
 geom_line(aes(x=threshold,y=CI_lower),lty=2,col="#C11432") + coord_fixed() +
  theme_bw() +
  theme(panel.spacing = unit(2, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

  
# calculate variance to give CI

# pick 0.95 as representative and calculate for every other location
threshold <- 0.95
chi_X_Y <- c()
var_X_Y <- c()

for (i in 1:nrow(uk_winter)) {
  X <- uk_winter[uk_winter$is_location=="london",6:ncol(uk_winter)] %>% as_vector()
  Y <- uk_winter[i,6:ncol(uk_winter)] %>% as_vector()
  df <- data.frame(X=X,Y=Y) %>% remove_rownames()
  # use PIT
  # create x and y quantile variables
  U <- df %>% select(X) %>% arrange(X) %>% mutate(u=row_number()/(nrow(df)+1))
  V <- df %>% select(Y) %>% arrange(Y) %>% mutate(v=row_number()/(nrow(df)+1)) 
  df <- df %>% left_join(U,by="X") %>% left_join(V,by="Y")
  chi_X_Y[i] <- chi(df,threshold)
  cu <- mean(df$u < threshold & df$v < threshold)
  cnst <- qnorm((1 + 0.95)/2)
  varchi <- ((1/log(threshold)^2 * 1)/cu^2 * cu * (1 - cu))/length(df$u)
  varchi <- cnst * sqrt(varchi)
  var_X_Y[i] <- varchi
}

df_X_Y <- data.frame(chi_X_Y,CI_u=chi_X_Y+var_X_Y,CI_l=chi_X_Y-var_X_Y)
df_X_Y <- cbind(uk_winter[,c(1:2,3)],df_X_Y)
p1 <- ggplot(df_X_Y) + ylim(c(-0.2,1)) + geom_point(aes(x=dist_london,y=chi_X_Y),size=0.5)+ 
  geom_point(aes(x=dist_london,y=CI_u),size=0.5,col="#C11432") +
  ylab(TeX("$\\chi_{i,j}(.95)$")) + xlab("Distance from central London") +
  geom_point(aes(x=dist_london,y=CI_l),size=0.5,col="#C11432") +
  geom_line(aes(x=dist_london,y=CI_l),lty=2,col="#C11432") +
  geom_line(aes(x=dist_london,y=CI_u),lty=2,col="#C11432") +
  theme_bw() +
  theme(panel.spacing = unit(2, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

# birmingham
chi_X_Y <- c()
var_X_Y <- c()

for (i in 1:nrow(uk_winter)) {
  X <- uk_winter[uk_winter$is_location=="birmingham",6:ncol(uk_winter)] %>% as_vector()
  Y <- uk_winter[i,6:ncol(uk_winter)] %>% as_vector()
  df <- data.frame(X=X,Y=Y) %>% remove_rownames()
  # use PIT
  # create x and y quantile variables
  U <- df %>% select(X) %>% arrange(X) %>% mutate(u=row_number()/(nrow(df)+1))
  V <- df %>% select(Y) %>% arrange(Y) %>% mutate(v=row_number()/(nrow(df)+1)) 
  df <- df %>% left_join(U,by="X") %>% left_join(V,by="Y")
  chi_X_Y[i] <- chi(df,threshold)
  cu <- mean(df$u < threshold & df$v < threshold)
  cnst <- qnorm((1 + 0.95)/2)
  varchi <- ((1/log(threshold)^2 * 1)/cu^2 * cu * (1 - cu))/length(df$u)
  varchi <- cnst * sqrt(varchi)
  var_X_Y[i] <- varchi
}

df_X_Y <- data.frame(chi_X_Y,CI_u=chi_X_Y+var_X_Y,CI_l=chi_X_Y-var_X_Y)
df_X_Y <- cbind(uk_winter[,c(1:2,4)],df_X_Y)
p2 <- ggplot(df_X_Y) + ylim(c(-0.2,1)) + geom_point(aes(x=dist_birmingham,y=chi_X_Y),size=0.5)+ 
  geom_point(aes(x=dist_birmingham,y=CI_u),size=0.5,col="#C11432") +
  ylab(TeX("$\\chi_{i,j}(.95)$")) + xlab("Distance from central Birmingham") +
  geom_point(aes(x=dist_birmingham,y=CI_l),size=0.5,col="#C11432") +
  geom_line(aes(x=dist_birmingham,y=CI_l),lty=2,col="#C11432") +
  geom_line(aes(x=dist_birmingham,y=CI_u),lty=2,col="#C11432") +
  theme_bw() +
  theme(panel.spacing = unit(2, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

# glasgow
chi_X_Y <- c()
var_X_Y <- c()

for (i in 1:nrow(uk_winter)) {
  X <- uk_winter[uk_winter$is_location=="glasgow",6:ncol(uk_winter)] %>% as_vector()
  Y <- uk_winter[i,6:ncol(uk_winter)] %>% as_vector()
  df <- data.frame(X=X,Y=Y) %>% remove_rownames()
  # use PIT
  # create x and y quantile variables
  U <- df %>% select(X) %>% arrange(X) %>% mutate(u=row_number()/(nrow(df)+1))
  V <- df %>% select(Y) %>% arrange(Y) %>% mutate(v=row_number()/(nrow(df)+1)) 
  df <- df %>% left_join(U,by="X") %>% left_join(V,by="Y")
  chi_X_Y[i] <- chi(df,threshold)
  cu <- mean(df$u < threshold & df$v < threshold)
  cnst <- qnorm((1 + 0.95)/2)
  varchi <- ((1/log(threshold)^2 * 1)/cu^2 * cu * (1 - cu))/length(df$u)
  varchi <- cnst * sqrt(varchi)
  var_X_Y[i] <- varchi
}

df_X_Y <- data.frame(chi_X_Y,CI_u=chi_X_Y+var_X_Y,CI_l=chi_X_Y-var_X_Y)
df_X_Y <- cbind(uk_winter[,c(1:2,5)],df_X_Y)
p3 <- ggplot(df_X_Y) + ylim(c(-0.2,1)) + geom_point(aes(x=dist_glasgow,y=chi_X_Y),size=0.5)+ 
  geom_point(aes(x=dist_glasgow,y=CI_u),size=0.5,col="#C11432") +
  ylab(TeX("$\\chi_{i,j}(.95)$")) + xlab("Distance from central Glasgow") +
  geom_point(aes(x=dist_glasgow,y=CI_l),size=0.5,col="#C11432") +
  geom_line(aes(x=dist_glasgow,y=CI_l),lty=2,col="#C11432") +
  geom_line(aes(x=dist_glasgow,y=CI_u),lty=2,col="#C11432") +
  theme_bw() +
  theme(panel.spacing = unit(2, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

grid.arrange(p1,p2,p3,ncol=3)

# add other seasons

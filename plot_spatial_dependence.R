library(tidyverse)
library(gridExtra)
library(evd)
library(latex2exp)
uk_winter <- readRDS("data/uk_1999_2018_winter.RDS")
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

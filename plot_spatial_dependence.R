library(tidyverse)
library(gridExtra)
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
# create empirical chi(u) function
chi <- function(df,threshold) {
  2-log( mean(df$u < threshold & df$v < threshold))/log(mean(df$u < threshold))  
}
chi(df,threshold=0.95)
    
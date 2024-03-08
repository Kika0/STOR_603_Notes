library(MASS)
library(tidyverse)
library(latex2exp)
# set theme defaults to be black rectangle
theme_set(theme_bw())
theme_replace(
  panel.spacing = unit(2, "lines"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA) )

gaussprocess <- function(from = 0, to = 1, K = function(s, t) {exp(-(abs(s-t)/lambda)^alpha)},
                         start = NULL, m = 1000,alpha=1,lambda=1) {
  
  t <- seq(from = from, to = to, length.out = m)
  Sigma <- sapply(t, function(s1) {
    sapply(t, function(s2) {
      K(s1, s2)
    })
  })
  
  path <- mvrnorm(mu = rep(0, times = m), Sigma = Sigma)
  if (!is.null(start)) {
    path <- path - path[1] + start  # Must always start at "start"
  }
  
  return(data.frame("t" = t, "xt" = path))
}

set.seed(1)
l1 <- gaussprocess(lambda=0.1)
set.seed(1)
l2 <- gaussprocess(lambda=1)
set.seed(1)
l3 <- gaussprocess(lambda=10)
df <- cbind(rbind(l1,l2,l3),
            data.frame(lambda=c(rep("0.1",1000),rep("1",1000),rep("10",1000))))

ggplot(df) + geom_line(aes(x=t,y=xt,col=lambda))+ ylab(TeX(paste0("$X($","$s$","$)")))


set.seed(1)
l1 <- gaussprocess(alpha = 0.1)
set.seed(1)
l2 <- gaussprocess(alpha = 1)
set.seed(1)
l3 <- gaussprocess(alpha = 1.9)
df <- cbind(rbind(l1,l2,l3),
            data.frame(alpha=c(rep("0.1",1000),rep("1",1000),rep("1.9",1000))))

ggplot(df) + geom_line(aes(x=t,y=xt,col=alpha))+ ylab(TeX(paste0("$X($","$s$","$)")))

tmp <- data.frame(t=as.numeric(),xt=as.numeric())
tmp1 <- data.frame(t=as.numeric(),xt=as.numeric(),ite=as.character())
set.seed(123)
for (i in 1:1000) {
  from <- (i-1)*100+1
  to <- 100*i
  tmp[from:to,] <- gaussprocess(m=100)
  tmp1[from:to,] <- cbind(tmp[from:to,],data.frame(ite=rep(as.character(i),100)))
}

ggplot(tmp1) + geom_line(aes(x=t,y=xt,col=ite),alpha=0.3,linewidth=0.1)+ylab(TeX(paste0("$X($","$s$","$)")))+   theme(legend.position="none")
plot(density(tmp1$xt[tmp1$t==3/99]))
plot(density(rnorm(1000)))
ggplot() + geom_density(tmp1 %>% mutate(t=as.character(t)), mapping=aes(x = xt, col = t),alpha = 0.1,linewidth=0.3)+ theme(legend.position="none") +
  geom_density(data.frame(x=rnorm(1000)),mapping=aes(x=x))

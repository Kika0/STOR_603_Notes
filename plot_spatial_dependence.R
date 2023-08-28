library(tidyverse)
library(gridExtra)
library(evd)
library(latex2exp)
library(gridExtra)
library(xtable)
library(units)
uk <- list()
uk[[1]] <- readRDS("data/uk_1999_2018_winter.RDS")
uk[[2]] <- readRDS("data/uk_1999_2018_spring.RDS")
uk[[3]] <- readRDS("data/uk_1999_2018_summer.RDS")
uk[[4]] <- readRDS("data/uk_1999_2018_autumn.RDS")
# set theme defaults to be black rectangle
theme_set(theme_bw())
theme_replace(
          panel.spacing = unit(2, "lines"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA) )

# for a nice table in latex ----
# for_latex <- uk_spring[c(1:3,76:78,160:162,347:349,445),c(1:7,1806)] %>%
#   rowid_to_column() 
#  print(xtable(for_latex), include.rownames=FALSE

# create empirical chi(u) function
chi <- function(df,threshold) {
   # 2-log( mean(df$u < threshold & df$v < threshold))/log(mean(df$u < threshold)) 
  ( mean(df$u > threshold & df$v > threshold))/(1-threshold)
    #mean(df$u > threshold))
}

plot_dependence <- function(city="london",Y=147,season=3) {
  #specify seasonal colour
  season_colour <- c("#009ADA","#66A64F","#C11432","#DF5D22")
  CI_col <- season_colour[season]
  tmp <- uk[[season]]
  X <- tmp[tmp$is_location==city,7:ncol(tmp)] %>% as_vector()
  Y <- tmp[Y,7:ncol(tmp)] %>% as_vector()
  df <- data.frame(X=X,Y=Y) %>% rowid_to_column()
  # use PIT
  # create x and y quantile variables
  U <- df %>% select(rowid,X) %>% arrange(X) %>% mutate(u=row_number()/(nrow(df)+1))
  V <- df %>% select(rowid,Y) %>% arrange(Y) %>% mutate(v=row_number()/(nrow(df)+1)) 
  df <- df %>% select(rowid) %>%  left_join(U,by="rowid") %>% left_join(V,by="rowid")
  
  # plot
  # p1 <- ggplot(df) + geom_point(aes(x=X,y=Y),size=0.1) + xlab(TeX("$Y(s_1)$ (London temperature)")) + ylab(TeX("$Y(s_2)$"))+
  #   ggtitle("Summer temperature time series")
  # p2 <- ggplot(df) + geom_point(aes(x=u,y=v),size=0.1) + xlab("u") +
  #   coord_fixed() + ggtitle("Uniform transform")
  # grid.arrange(p1,p2,ncol=2)

  threshold <- sort(df$u)[13:(length(df$u)-6)]
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
  df_chi$CI_lower[is.na( df_chi$CI_lower)] <- 0
  df_chi$CI_upper[is.na( df_chi$CI_upper)] <- 1
  df_chi$CI_lower[df_chi$CI_lower<(0)] <- rep(0,sum(df_chi$CI_lower<(0)))
  df_chi$CI_upper[df_chi$CI_upper>1] <- rep(1,sum(df_chi$CI_upper>1))
  p1 <- ggplot(df_chi)  + geom_line(aes(x=threshold,y=chi_u))+ 
    geom_line(aes(x=threshold,y=CI_upper),lty=2,col=CI_col) +
    ylab(TeX("$\\chi_u(s_1,s_2)$")) + xlab("u") +
    geom_line(aes(x=threshold,y=CI_lower),lty=2,col=CI_col) + coord_fixed() +
    geom_ribbon(aes(x=u,ymin=CI_lower,ymax=CI_upper), fill=CI_col, alpha=0.2)+ ylim(c(0,1))
  return(p1)
}
p1 <- plot_dependence(city="london",Y=147,season=1)
p2 <- plot_dependence(city="birmingham",Y=228,season=1)
p3 <- plot_dependence(city="glasgow",Y=310,season=1)
p4 <- plot_dependence(city="london",Y=147,season=2)
p5 <- plot_dependence(city="birmingham",Y=228,season=2)
p6 <- plot_dependence(city="glasgow",Y=310,season=2)
p7 <- plot_dependence(city="london",Y=147,season=3)
p8 <- plot_dependence(city="birmingham",Y=228,season=3)
p9 <- plot_dependence(city="glasgow",Y=310,season=3)
p10 <- plot_dependence(city="london",Y=147,season=4)
p11 <- plot_dependence(city="birmingham",Y=228,season=4)
p12 <- plot_dependence(city="glasgow",Y=310,season=4)
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,ncol=3,nrow=4)

# plot dependence decaying with distance ----
plot_X_dist_dependence <- function(city="London",season=3,threshold=0.95,y_low_limit=1) {
  #specify seasonal colour
  season_colour <- c("#009ADA","#66A64F","#C11432","#DF5D22")
  CI_col <- season_colour[season]
  tmp <- uk[[season]]
  chi_X_Y <- c()
  var_X_Y <- c()
  tmp <- uk[[season]]
  for (i in 1:nrow(tmp)) {
    X <- tmp[tmp$is_location==tolower(city),7:ncol(tmp)] %>% as_vector()
    Y <- tmp[i,7:ncol(tmp)] %>% as_vector()
    df <- data.frame(X=X,Y=Y) %>% rowid_to_column()
    # use PIT
    # create x and y quantile variables
    U <- df %>% select(rowid,X) %>% arrange(X) %>% mutate(u=row_number()/(nrow(df)+1))
    V <- df %>% select(rowid,Y) %>% arrange(Y) %>% mutate(v=row_number()/(nrow(df)+1)) 
    df <- df %>% select(rowid) %>%  left_join(U,by="rowid") %>% left_join(V,by="rowid")
    chi_X_Y[i] <- chi(df,threshold)
    cu <- mean(df$u < threshold & df$v < threshold)
    cnst <- qnorm((1 + 0.95)/2)
    varchi <- ((1/log(threshold)^2 * 1)/cu^2 * cu * (1 - cu))/length(df$u)
    varchi <- cnst * sqrt(varchi)
    var_X_Y[i] <- varchi
  }
  
  df_X_Y <- data.frame(chi_X_Y,CI_u=chi_X_Y+var_X_Y,CI_l=chi_X_Y-var_X_Y)
  df_X_Y <- cbind(tmp[,c(1:6)],df_X_Y)
  df_X_Y$dist_london <- drop_units(df_X_Y$dist_london)/1000
  df_X_Y$dist_birmingham <- drop_units(df_X_Y$dist_birmingham)/1000
  df_X_Y$dist_glasgow <- drop_units(df_X_Y$dist_glasgow)/1000
  df_X_Y$CI_l[is.na( df_X_Y$CI_l)] <- y_low_limit
  df_X_Y$CI_u[is.na( df_X_Y$CI_u)] <- 1
  df_X_Y$CI_l[df_X_Y$CI_l<(y_low_limit)] <- rep(y_low_limit,sum(df_X_Y$CI_l<(y_low_limit)))
  df_X_Y$CI_u[df_X_Y$CI_u>1] <- rep(1,sum(df_X_Y$CI_u>1))
  dist_min <- min(c(df_X_Y$dist_london,df_X_Y$dist_birmingham,df_X_Y$dist_glasgow))
  dist_max <- max(c(df_X_Y$dist_london,df_X_Y$dist_birmingham,df_X_Y$dist_glasgow))
  km <- TeX(" $[km]$")
  dist_city <- paste0("dist_",tolower(city))
  p1 <- ggplot(df_X_Y) + ylim(c(y_low_limit,1)) + xlim(c(dist_min,dist_max))+ geom_point(aes(x=.data[[dist_city]],y=chi_X_Y),size=0.5)+ 
    #geom_point(aes(x=dist_london,y=CI_u),size=0.5,col=CI_col) +
    ylab(TeX(paste0("$\\chi_{.",threshold*100,"}(s_1,s_2)$"))) + xlab(TeX(paste0("Distance from ", city ," $[km]$"))) +
    #geom_point(aes(x=dist_london,y=CI_l),size=0.5,col=CI_col) +
    geom_line(aes(x=.data[[dist_city]] ,y=CI_l),lty=2,col=CI_col) +
    geom_line(aes(x=.data[[dist_city]],y=CI_u),lty=2,col=CI_col) +
    geom_ribbon(aes(x=.data[[dist_city]],ymin=CI_l,ymax=CI_u), fill=CI_col, alpha=0.2)
  return(p1)
}

y_low_limit <- 0
threshold <- 0.99
p1 <-plot_X_dist_dependence(city="London",season=1,threshold=threshold,y_low_limit=y_low_limit)
p2 <-plot_X_dist_dependence(city="Birmingham",season=1,threshold=threshold,y_low_limit=y_low_limit)
p3 <-plot_X_dist_dependence(city="Glasgow",season=1,threshold=threshold,y_low_limit=y_low_limit)
p4 <-plot_X_dist_dependence(city="London",season=2,threshold=threshold,y_low_limit=y_low_limit)
p5 <-plot_X_dist_dependence(city="Birmingham",season=2,threshold=threshold,y_low_limit=y_low_limit)
p6 <-plot_X_dist_dependence(city="Glasgow",season=2,threshold=threshold,y_low_limit=y_low_limit)
p7 <-plot_X_dist_dependence(city="London",season=3,threshold=threshold,y_low_limit=y_low_limit)
p8 <-plot_X_dist_dependence(city="Birmingham",season=3,threshold=threshold,y_low_limit=y_low_limit)
p9 <-plot_X_dist_dependence(city="Glasgow",season=3,threshold=threshold,y_low_limit=y_low_limit)
p10 <-plot_X_dist_dependence(city="London",season=4,threshold=threshold,y_low_limit=y_low_limit)
p11 <-plot_X_dist_dependence(city="Birmingham",season=4,threshold=threshold,y_low_limit=y_low_limit)
p12 <-plot_X_dist_dependence(city="Glasgow",season=4,threshold=threshold,y_low_limit=y_low_limit)
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,ncol=3,nrow=4)

# plot winter minima
uk[[1]][,7:ncol(uk[[1]])] <- (-uk[[1]][,(7:ncol(uk[[1]]))] )
p1 <-plot_X_dist_dependence(city="London",season=1)
p2 <-plot_X_dist_dependence(city="Birmingham",season=1)
p3 <-plot_X_dist_dependence(city="Glasgow",season=1)
grid.arrange(p1,p2,p3,ncol=3)
p1 <- plot_dependence(city="london",Y=147,season=1)
p2 <- plot_dependence(city="birmingham",Y=228,season=1)
p3 <- plot_dependence(city="glasgow",Y=310,season=1)

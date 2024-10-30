library(MASS)
library(tidyverse)
library(latex2exp)
library(viridis)
library(plgp)
library(gridExtra)
library(here)
library(tmap)
library(units)
library(sf)

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

#' Gaussian process simulation in one dimension
#'
#' @param from lower bound
#' @param to upper bound
#' @param K kernel function, default is exponential family
#' @param start initial value
#' @param m number of points equally spaced between lower and upper bound
#' @param alpha parameter of kernel function
#' @param lambda parameter of kernel function
#'
#' @return
#' @export
#'
#' @examples
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

# sample in two dimensions ----
gaussprocess2d <- function(from = 0, to = 1,alpha=1,lambda=1,rho=NULL, sig=NULL, K = function(s, t) {exp(-(sqrt((s[1]-t[1])^2+(s[2]-t[2])^2)/lambda)^alpha)},
                         start = NULL,  m = 10) {
  
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

# UKCP 18 data (summer max daily temperatures 1999-2018) ----
ukcp18 <- readRDS("data/uk_1999_2018_summer.RDS") %>% relocate(dist_london,.after=dist_glasgow)
# remove last year of the data due to error (same data as first year)
ukcp18 <- ukcp18[,1:1716]
# remove mean from the data
# try shifting all the data with the corresponding yearly coefficients

y <- as.vector(unlist(ukcp18[sample(1:445,size=50),7:1716])) # try pooling all sites together
I <- seq_along(y)
gradient <- summary(lm(y~I))$coefficients[2,1]

shift_tmp <- c()
for (i in seq_along(y)) {
  shift_tmp[i] <- y[i] - (i-1)*gradient 
}

plot(I,y)
plot(I,shift_tmp)
ggplot() + geom_point(data.frame(I=I,y=y),mapping=aes(x=I,y=y),size=0.1,alpha=0.5)+ geom_point(data.frame(I=I,shift_tmp=shift_tmp),mapping = aes(x=I,y=shift_tmp),col="#C11432",size=0.1,alpha=0.5)  + geom_smooth(data.frame(I=I,y=y),mapping=aes(x=I,y=y))
conv <- CnvRttPol(latlon = data.frame(long=ukcp18$Longitude,lat=ukcp18$Latitude),spol_coor = c(gr_npole_lon, gr_npole_lat))
uk_sf_rot <- data.frame(lon=conv$lon,lat=conv$lat) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)
uk_sf_rot <- st_transform(uk_sf_rot,27700)
# coast_dist <- st_geometry(obj = uk) %>%
#   st_cast(to = 'LINESTRING') %>%
#   st_distance(y=uk_sf_rot)
# add buffer
uk_buffered <- st_buffer(uk, 50000)   
# polygon(s) of only the buffer
buffer_only <- st_difference(uk_buffered, uk) 
tm_shape(buffer_only) + tm_polygons()
coast_dist <- st_distance(uk_sf_rot,buffer_only)
uk_sf_rot <- cbind(uk_sf_rot,data.frame("coast_dist" = coast_dist))
tm_shape(uk_sf_rot) + tm_dots("coast_dist",size=1)
ukcp18 <- ukcp18 %>% mutate("Longitude" = conv$lon, "Latitude" = conv$lat, "coast_dist" = uk_sf_rot$coast_dist) %>% relocate(coast_dist,.before = dist_birmingham)
sims <- ukcp18 %>% dplyr::arrange(is_location)%>% dplyr::select(!contains("i")) %>% t() %>% as.data.frame()
# ordered alphabetically so Y1 Birmingham, Y2 Glasgow and Y3 is London
colnames(sims) <- paste0("Y",1:ncol(sims))
# transform to Laplace margins
sims <- as.data.frame((sims %>% apply(c(2),FUN=row_number))/(nrow(sims)+1)) %>% apply(c(1,2),FUN=unif_laplace_pit) %>% as.data.frame()
# calculate the residuals for Birmingham (1), then Glasgow (2) and London (3)
sites <- c("Birmingham","Glasgow","London")
uk_tmp3 <- data.frame("lik" = numeric(),"lika" = numeric() ,"likb" = numeric(),"lik2" = numeric(),
                      "a" = numeric() , "b" = numeric(),
                      "mu" = numeric() ,"mu_agg" = numeric(),
                      "sig" = numeric() ,"sig_agg" = numeric(),"sigl" = numeric(),"sigu" = numeric(),
                      "delta" = numeric(),"deltal" = numeric(),"deltau" = numeric(),
                      "given" = numeric(), "res" = numeric(),
                      "margin" = character(), "method" =character(), "cond_site" = character() )
v <- 0.98
for (cond_site in 1:3) {
cond_var <- cond_site
tmp_est <- par_est(sims,v=v,given=c(cond_var),margin = "AGG", method="two_step")
tmp_est$pair_dist <- ukcp18 %>% arrange(is_location) %>% filter(is_location != tolower(sites[cond_var])) %>%  dplyr::select(4+cond_var) %>% pull()
tmp_est$coast_dist <- ukcp18 %>% arrange(is_location) %>% filter(is_location != tolower(sites[cond_var])) %>%  dplyr::select(4) %>% pull()
tmp <- tmp_est %>% mutate(given=factor(given,levels = cond_var))
tmp1 <- tmp %>% add_row(.before=cond_var)
# match back to spatial locations and plot
uk_tmp <- uk_temp_sf %>% dplyr::select() %>% cbind(ukcp18[,1:8]) %>% 
   arrange(is_location) 
uk_tmp1 <- cbind(uk_tmp,tmp1) %>% mutate(margin=rep("AGG",nrow(uk_tmp)),method=rep("two_step",nrow(uk_tmp)))

tmp_est <- par_est(sims,v=v,given=c(cond_var),margin = "AGG", method="one_step")
tmp_est$pair_dist <- ukcp18 %>% arrange(is_location) %>% filter(is_location != tolower(sites[cond_var])) %>%  dplyr::select(4+cond_var) %>% pull()
tmp_est$coast_dist <- ukcp18 %>% arrange(is_location) %>% filter(is_location != tolower(sites[cond_var])) %>%  dplyr::select(4) %>% pull()
tmp <- tmp_est %>% mutate(given=factor(given,levels = cond_var))
tmp1 <- tmp %>% add_row(.before = cond_var)
# match back to spatial locations and plot
uk_tmp <- uk_temp_sf %>% dplyr::select() %>% cbind(ukcp18[,1:8]) %>% 
  arrange(is_location) 
uk_tmp2 <- rbind(cbind(uk_tmp,tmp1) %>% mutate(margin=rep("AGG",nrow(uk_tmp)),method=rep("one_step",nrow(uk_tmp))),uk_tmp1)

tmp_est <- par_est(sims,v=v,given=c(cond_var),margin = "AGG", method="sequential")
tmp_est$pair_dist <- NA
# tmp_est$pair_dist <- ukcp18 %>% arrange(is_location) %>% filter(is_location != tolower(sites[cond_var])) %>%  dplyr::select(4+cond_var) %>% pull()
tmp_est$coast_dist <- ukcp18 %>% arrange(is_location) %>% filter(is_location != tolower(sites[cond_var])) %>%  dplyr::select(4) %>% pull()
tmp <- tmp_est %>% mutate(given=factor(given,levels = cond_var))
tmp1 <- tmp %>% add_row(.before=cond_var)
# match back to spatial locations and plot
uk_tmp <- uk_temp_sf %>% dplyr::select() %>% cbind(ukcp18[,1:8]) %>% 
  arrange(is_location) 
uk_tmp3 <- rbind(uk_tmp3,rbind(cbind(uk_tmp,tmp1) %>% mutate(margin=rep("AGG",nrow(uk_tmp)),method=rep("sequential",nrow(uk_tmp))),uk_tmp2) %>% 
  mutate(cond_site = sites[cond_var]))  
}

pa <- tmap_arrange(tm_shape(uk_tmp3 %>% filter(given==cond_var & method=="two_step")) + tm_dots(col="a",style="cont",size=0.3,palette="viridis",title=TeX("$\\alpha$")) + tm_layout(main.title="AGG 2 step"),
             tm_shape(uk_tmp3 %>% filter(given==cond_var & method=="sequential")) + tm_dots(col="a",style="cont",size=0.3,palette="viridis",title=TeX("$\\alpha$")) + tm_layout(main.title=TeX("$\\beta=0 \\rightarrow \\hat{\\alpha} \\rightarrow \\hat{\\beta}$")),             
             tm_shape(uk_tmp3 %>% filter(given==cond_var & method=="one_step")) + tm_dots(col="a",style="cont",size=0.3,palette="viridis",title=TeX("$\\alpha$")) + tm_layout(main.title="AGG 1 step"),ncol=3)
tmap::tmap_save(tm = pa, filename = paste0("plots/map_a_agg3methods",sites[cond_var],".png"),height =6,width=9)
pb <- tmap_arrange(tm_shape(uk_tmp3 %>% filter(given==cond_var & method=="two_step")) + tm_dots(col="b",style="cont",size=0.3,palette="viridis",title=TeX("$\\beta$"))+ tm_layout(main.title="AGG 2 step"),
             tm_shape(uk_tmp3 %>% filter(given==cond_var & method=="sequential")) + tm_dots(col="b",style="cont",size=0.3,palette="viridis",title=TeX("$\\beta$")) + tm_layout(main.title=TeX("$\\beta=0 \\rightarrow \\hat{\\alpha} \\rightarrow \\hat{\\beta}$")),             
             tm_shape(uk_tmp3 %>% filter(given==cond_var & method=="one_step")) + tm_dots(col="b",style="cont",size=0.3,palette="viridis",title=TeX("$\\beta$"))+ tm_layout(main.title="AGG 1 step"),ncol=3)
tmap::tmap_save(tm = pb, filename = paste0("plots/map_b_agg3methods",sites[cond_var],".png"),height = 6,width=9)
pmu <- tmap_arrange(tm_shape(uk_tmp3 %>% filter(given==cond_var & method=="two_step")) + tm_dots(col="mu_agg",style="cont",size=0.3,palette="viridis",title=TeX("$\\mu_{AGG}$"))+ tm_layout(main.title="AGG 2 step"),
             tm_shape(uk_tmp3 %>% filter(given==cond_var & method=="sequential")) + tm_dots(col="mu_agg",style="cont",size=0.3,palette="viridis",title=TeX("$\\mu_{AGG}$")) + tm_layout(main.title=TeX("$\\beta=0 \\rightarrow \\hat{\\alpha} \\rightarrow \\hat{\\beta}$")),             
             tm_shape(uk_tmp3 %>% filter(given==cond_var & method=="one_step")) + tm_dots(col="mu_agg",style="cont",size=0.3,palette="viridis",title=TeX("$\\mu_{AGG}$"))+ tm_layout(main.title="AGG 1 step"),ncol=3)
tmap::tmap_save(tm = pmu, filename = paste0("plots/map_mu_agg3methods",sites[cond_var],".png"),height = 6,width=9)
psig <- tmap_arrange(tm_shape(uk_tmp3 %>% filter(given==cond_var & method=="two_step")) + tm_dots(col="sig_agg",style="cont",size=0.3,palette="viridis",title=TeX("$\\sigma_{AGG}$"))+ tm_layout(main.title="AGG 2 step"),
             tm_shape(uk_tmp3 %>% filter(given==cond_var & method=="sequential")) + tm_dots(col="sig_agg",style="cont",size=0.3,palette="viridis",title=TeX("$\\sigma_{AGG}$")) + tm_layout(main.title=TeX("$\\beta=0 \\rightarrow \\hat{\\alpha} \\rightarrow \\hat{\\beta}$")),             
             tm_shape(uk_tmp3 %>% filter(given==cond_var & method=="one_step")) + tm_dots(col="sig_agg",style="cont",size=0.3,palette="viridis",title=TeX("$\\sigma_{AGG}$"))+ tm_layout(main.title="AGG 1 step"),ncol=3)
#n=10,style="quantile",
tmap::tmap_save(tm = psig, filename = paste0("plots/map_sig_agg3methods",sites[cond_var],".png"),height =6,width=9)
pdeltal <- tmap_arrange(tm_shape(uk_tmp3 %>% filter(given==cond_var & method=="two_step")) + tm_dots(col="deltal",n=10,style="quantile",size=0.3,palette="viridis",title=TeX("$\\delta_l$"))+ tm_layout(main.title="AGG 2 step"),
             tm_shape(uk_tmp3 %>% filter(given==cond_var & method=="sequential")) + tm_dots(col="deltal",size=0.3,n=10,style="quantile",palette="viridis",title=TeX("$\\delta_l$")) + tm_layout(main.title=TeX("$\\beta=0 \\rightarrow \\hat{\\alpha} \\rightarrow \\hat{\\beta}$")),             
             tm_shape(uk_tmp3 %>% filter(given==cond_var & method=="one_step")) + tm_dots(col="deltal",size=0.3,n=10,style="quantile",palette="viridis",title=TeX("$\\delta_l$"))+ tm_layout(main.title="AGG 1 step"),ncol=3)
tmap::tmap_save(tm = pdeltal, filename = paste0("plots/map_deltal_agg3methods",sites[cond_var],".png"),height = 6,width=9)
pdeltau <- tmap_arrange(tm_shape(uk_tmp3 %>% filter(given==cond_var & method=="two_step")) + tm_dots(col="deltau",size=0.3,n=10,style="quantile",palette="viridis",title=TeX("$\\delta_u$"))+ tm_layout(main.title="AGG 2 step"),
             tm_shape(uk_tmp3 %>% filter(given==cond_var & method=="sequential")) + tm_dots(col="deltau",size=0.3,n=10,style="quantile",palette="viridis",title=TeX("$\\delta_u$")) + tm_layout(main.title=TeX("$\\beta=0 \\rightarrow \\hat{\\alpha} \\rightarrow \\hat{\\beta}$")),             
             tm_shape(uk_tmp3 %>% filter(given==cond_var & method=="one_step")) + tm_dots(col="deltau",size=0.3,n=10,style="quantile",palette="viridis",title=TeX("$\\delta_u$"))+ tm_layout(main.title="AGG 1 step"),ncol=3)
tmap::tmap_save(tm = pdeltau, filename = paste0("plots/map_deltau_agg3methods",sites[cond_var],".png"),height = 6,width=9)

uk_tmp3 <- uk_tmp3 %>% mutate(cond_site = factor(cond_site,levels = sites))
# plot parameter estimates with pairwise distance from the conditioning site
p1 <- ggplot(uk_tmp3) + geom_point(aes(x=pair_dist,y=abs(lik),color=cond_site)) + scale_color_manual(values=c("#C11432","#009ADA","#66A64F")) + facet_wrap(~method,nrow=1)
p2 <- ggplot(uk_tmp3) + geom_point(aes(x=pair_dist,y=a,color=cond_site)) + scale_color_manual(values=c("#C11432","#009ADA","#66A64F")) + facet_wrap(~method,nrow=1)
p3 <- ggplot(uk_tmp3) + geom_point(aes(x=pair_dist,y=b,color=cond_site)) + scale_color_manual(values=c("#C11432","#009ADA","#66A64F")) + facet_wrap(~method,nrow=1)
ggsave(grid.arrange(p1,p2,p3,ncol=1), filename = "plots/agglikab.png")

p4 <-  ggplot(uk_tmp3) + geom_point(aes(x=pair_dist,y=mu_agg,color=cond_site)) + scale_color_manual(values=c("#C11432","#009ADA","#66A64F"))+ facet_wrap(~method,nrow=1)
p5 <-  ggplot(uk_tmp3) + geom_point(aes(x=pair_dist,y=sig_agg,color=cond_site)) + scale_color_manual(values=c("#C11432","#009ADA","#66A64F"))+ facet_wrap(~method,nrow=1)
ggsave(grid.arrange(p4,p5,ncol=1), filename = "plots/aggmusig.png")

p6 <-  ggplot(uk_tmp3) + geom_point(aes(x=pair_dist,y=deltal,color=cond_site)) + scale_color_manual(values=c("#C11432","#009ADA","#66A64F")) +ylim(c(0,5))+ facet_wrap(~method,nrow=1)
p7 <-  ggplot(uk_tmp3) + geom_point(aes(x=pair_dist,y=deltau,color=cond_site)) + scale_color_manual(values=c("#C11432","#009ADA","#66A64F"))+ylim(c(0,5)) + facet_wrap(~method,nrow=1)
# p8 <- ggplot(tmp) + geom_point(aes(x=pair_dist,y=(deltal-deltau))) + ylim(c(-2,5))
# p9 <- ggplot(tmp) + geom_point(aes(x=pair_dist,y=(deltal+deltau)/2)) +ylim(c(0,5))
ggsave(grid.arrange(p6,p7,ncol=1), filename = "plots/aggdeltas.png")

# plot parameter estimates with distance from the coast
p1 <- ggplot(uk_tmp3) + geom_point(aes(x=coast_dist,y=abs(lik),color=cond_site)) + scale_color_manual(values=c("#C11432","#009ADA","#66A64F")) + facet_wrap(~method,nrow=1)
p2 <- ggplot(uk_tmp3) + geom_point(aes(x=coast_dist,y=a,color=cond_site)) + scale_color_manual(values=c("#C11432","#009ADA","#66A64F")) + facet_wrap(~method,nrow=1)
p3 <- ggplot(uk_tmp3) + geom_point(aes(x=coast_dist,y=b,color=cond_site)) + scale_color_manual(values=c("#C11432","#009ADA","#66A64F")) + facet_wrap(~method,nrow=1)
ggsave(grid.arrange(p1,p2,p3,ncol=1), filename = "plots/agglikabc.png")

p4 <-  ggplot(uk_tmp3) + geom_point(aes(x=coast_dist,y=mu_agg,color=cond_site)) + scale_color_manual(values=c("#C11432","#009ADA","#66A64F"))+ facet_wrap(~method,nrow=1)
p5 <-  ggplot(uk_tmp3) + geom_point(aes(x=coast_dist,y=sig_agg,color=cond_site)) + scale_color_manual(values=c("#C11432","#009ADA","#66A64F"))+ facet_wrap(~method,nrow=1)
ggsave(grid.arrange(p4,p5,ncol=1), filename = "plots/aggmusigc.png")

p6 <-  ggplot(uk_tmp3) + geom_point(aes(x=coast_dist,y=deltal,color=cond_site)) + scale_color_manual(values=c("#C11432","#009ADA","#66A64F")) +ylim(c(0,5))+ facet_wrap(~method,nrow=1)
p7 <-  ggplot(uk_tmp3) + geom_point(aes(x=coast_dist,y=deltau,color=cond_site)) + scale_color_manual(values=c("#C11432","#009ADA","#66A64F"))+ylim(c(0,5)) + facet_wrap(~method,nrow=1)
# p8 <- ggplot(tmp) + geom_point(aes(x=pair_dist,y=(deltal-deltau))) + ylim(c(-2,5))
# p9 <- ggplot(tmp) + geom_point(aes(x=pair_dist,y=(deltal+deltau)/2)) +ylim(c(0,5))
ggsave(grid.arrange(p6,p7,ncol=1),filename = "plots/aggdeltasc.png")

cond_var <- 50
tmp_est <- par_est(sims,v=0.9,given=c(cond_var),margin = "AGG", method="sequential")
tmp_est$pair_dist <- NA
# tmp_est$pair_dist <- ukcp18 %>% arrange(is_location) %>% filter(is_location != tolower(sites[cond_var])) %>%  dplyr::select(4+cond_var) %>% pull()
#tmp_est$coast_dist <- ukcp18 %>% arrange(is_location) %>% filter(is_location != tolower(sites[cond_var])) %>%  dplyr::select(4) %>% pull()
tmp_est$coast_dist <- NA
tmp <- tmp_est %>% mutate(given=factor(given,levels = cond_var))
tmp1 <- tmp %>% add_row(.before=cond_var)
# match back to spatial locations and plot
uk_tmp <- uk_temp_sf %>% dplyr::select() %>% cbind(ukcp18[,1:8]) %>% 
  arrange(is_location) 
uk_tmp1 <- cbind(uk_tmp,tmp1) %>% mutate(margin="AGG",method=rep("two_step",nrow(uk_tmp)))

cond_var <- 200
tmp_est <- par_est(sims,v=0.9,given=c(cond_var),margin = "AGG", method="sequential")
tmp_est$pair_dist <- NA
# tmp_est$pair_dist <- ukcp18 %>% arrange(is_location) %>% filter(is_location != tolower(sites[cond_var])) %>%  dplyr::select(4+cond_var) %>% pull()
#tmp_est$coast_dist <- ukcp18 %>% arrange(is_location) %>% filter(is_location != tolower(sites[cond_var])) %>%  dplyr::select(4) %>% pull()
tmp_est$coast_dist <- NA
tmp <- tmp_est %>% mutate(given=factor(given,levels = cond_var))
tmp1 <- tmp %>% add_row(.before = cond_var)
# match back to spatial locations and plot
uk_tmp <- uk_temp_sf %>% dplyr::select() %>% cbind(ukcp18[,1:8]) %>% 
  arrange(is_location) 
uk_tmp2 <- rbind(cbind(uk_tmp,tmp1) %>% mutate(margin=rep("AGG",nrow(uk_tmp)),method=rep("one_step",nrow(uk_tmp))),uk_tmp1)

cond_var <- 350
tmp_est <- par_est(sims,v=0.9,given=c(cond_var),margin = "AGG", method="sequential")
tmp_est$pair_dist <- NA
# tmp_est$pair_dist <- ukcp18 %>% arrange(is_location) %>% filter(is_location != tolower(sites[cond_var])) %>%  dplyr::select(4+cond_var) %>% pull()
#tmp_est$coast_dist <- ukcp18 %>% arrange(is_location) %>% filter(is_location != tolower(sites[cond_var])) %>%  dplyr::select(4) %>% pull()
tmp_est$coast_dist <- NA
tmp <- tmp_est %>% mutate(given=factor(given,levels = cond_var))
tmp1 <- tmp %>% add_row(.before=cond_var)
# match back to spatial locations and plot
uk_tmp <- uk_temp_sf %>% dplyr::select() %>% cbind(ukcp18[,1:8]) %>% 
  arrange(is_location) 
uk_tmp3 <- rbind(uk_tmp3,rbind(cbind(uk_tmp,tmp1) %>% mutate(margin=rep("AGG",nrow(uk_tmp)),method=rep("sequential",nrow(uk_tmp))),uk_tmp2) %>% 
                   mutate(cond_site = sites[cond_var]))

pa <- tmap_arrange(tm_shape(uk_tmp3 %>% filter(given==50 & method=="two_step")) + tm_dots(col="a",style="cont",size=0.3,palette="viridis",title=TeX("$\\alpha$")) + tm_layout(main.title="Cond site 50"),
                   tm_shape(uk_tmp3 %>% filter(given==200 & method=="sequential")) + tm_dots(col="a",style="cont",size=0.3,palette="viridis",title=TeX("$\\alpha$")) + tm_layout(main.title="Cond site 200"),             
                   tm_shape(uk_tmp3 %>% filter(given==350 & method=="one_step")) + tm_dots(col="a",style="cont",size=0.3,palette="viridis",title=TeX("$\\alpha$")) + tm_layout(main.title="Cond site 350"),ncol=3)
pa

# use parameteric form for a ----
# calculate distance from the 3 conditioning sites
# transform dataframe to include a vector of x (temperature) and d (distance from the conditioning site)
p <- list()
for (i in 1:3) {
  cond_var <- i
x <- par_est(sims,v=0.9,given=c(cond_var),margin = "AGG", method="sequential")$a
d <- (ukcp18 %>% arrange(is_location))[-cond_var,] %>% select(4+cond_var) %>% pull() %>% units::drop_units()
opt <- optim(par=c(0.01,0.1),fn=NLL_exp_norm_noise,x=x,d=d)
# plot a function of alpha against distance
a <- exp(-opt$par[1]*d)
p[[i]] <- ggplot() + geom_line(data=data.frame(x=d,y=a),aes(x=x,y=y)) +
  geom_point(data=data.frame(x=d,y=x),aes(x=x,y=y)) +
  ylab(TeX("$\\alpha$")) + xlab("Distance")
}
grid.arrange(p[[1]],p[[2]],p[[3]],ncol=3)
# simulate from the parametric form of alpha to compare with marginal fits

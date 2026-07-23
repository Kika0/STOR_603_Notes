library(tidyverse)
library(tmap)
library(tmaptools)
library(sf)
library(parallel)
library(lubridate)
library(gridExtra)
library(LaplacesDemon)
library(latex2exp)
library(evd)
theme_set(theme_bw())
theme_replace(
  panel.spacing = unit(2, "lines"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA) )
folder_name <- "../Documents/final_model_3/"
# load observed data
#source("spatial_parameter_estimation.R") # for spatial_par_est function
file.sources = list.files(pattern="*helpers.R")
sapply(file.sources,source)
load("data_processed/data_mod_Lap.RData",verbose = TRUE)
load("data_processed/spatial_helper.RData", verbose = TRUE)
source("simon/P2q_function_helpers.R")
load("data_processed/P2qselected_helpers.RData", verbose = TRUE)
q <- 0.9 # set quantile threshold
load(paste0("data_processed/N9000_sequential2_AGG_all12sites",q*100,".RData"),verbose = TRUE) # original parameter estimates
load("data_processed/iterative_phi0l_phi0u_estimates_London.RData",verbose=TRUE) # residual margin parameters for delta constants
load("data_processed/residual_dependence_pars.RData", verbose = TRUE) # residual dependence parameters
deltal <- result_new[[12]]$deltal[1]
deltau <- result_new[[12]]$deltau[1]
source("simon/final_model_helpers.R")

model3_wrapper <- function(site_i) {
cond_index <- df_sites[3,site_i]
cond_site_name <- names(df_sites)[site_i]
Nite_2_3 <- 10
Nite_phi <- 10
# Model 3: parameter estimation ------------------------------------------------
par_est_model_3 <- function(cond_index,v=0.9,data_Lap=data_mod_Lap,grid20km=xyUK20_sf,deltal,deltau,Nite_2_3,Nite_phi) {
  # calculate distance from the conditioning site
  dist_tmp <- as.numeric(unlist(st_distance(grid20km[cond_index,],grid20km[-cond_index,])))
  distnorm <- dist_tmp/1000000  # normalise distance using a common constant
  # subset conditioning dataset
  data_Lapv <- data_Lap %>% filter(data_Lap[,cond_index]>quantile(data_Lap[,cond_index],v))
  # initiate objects
  x2 <- x3 <- list() # a list of output at steps 2 and 3
  x2_df <- data.frame("mu_agg"=numeric(),"sigl"=numeric(),"sigu"=numeric(),"phi0u"=numeric(),"phi1u"=numeric(),"phi2u"=numeric(),"phi0l"=numeric(),"phi1l"=numeric(),"phi2l"=numeric(),"iteration"=numeric())
  x3_df <- data.frame("a"=numeric(),"b"=numeric(),"mu"=numeric(),"iteration"=numeric())
  x_time <- data.frame("step2"=numeric(),"step3"=numeric(),"iteration"=numeric())
  # 1. estimate alpha and beta
  x1 <- par_est(df=data_Lap,v = v,given = cond_index,margin = "Normal",method = "sequential2",keef_constraints = c(1,2))
  # iterate between steps 2. and 3.
  for (Nite_i in 1:Nite_2_3) {
  if (Nite_i>1) {
    a <- tmp$a[!is.na(tmp$a)]
    b <- tmp$b[!is.na(tmp$b)]
    pe_res <- pe_phi_ite %>% dplyr::select(sigl,sigu) %>% na.omit() 
  } else {
    a <- x1$a[!is.na(x1$a)]
    b <- x1$b[!is.na(x1$b)]
    pe_res <- data.frame("sigl" = rep(1,length(a)),"sigu" = rep(1,length(a)))
  }
  # update residuals
  Z <- observed_residuals(df=data_Lap,given=cond_index,v = v,a=a,b=b)
  print(v) # check value of quantile threshold
  # 2. estimate phis
  x_time_dummy <- Sys.time()
  x2 <- par_est_ite(z=Z,given=cond_index,cond_site_dist = distnorm, parest_site = pe_res,Nite = Nite_phi,show_ite=TRUE,deltal = deltal,deltau = deltau) 
  x_time2 <- Sys.time()-x_time_dummy
  x2_df <- rbind(x2_df,x2[[12]] %>% mutate("iteration"=Nite_i))
  # 3. reestimate a,b,mu
  pe_phi_ite <- x2[[12]] %>% dplyr::select(mu_agg,sigl,sigu,deltal,deltau)
  x_time_dummy <- Sys.time()
  x3 <- sapply(1:ncol(data_Lap),FUN=NLL_AGG_wrapper,data_Lapv=data_Lapv,cond_index=cond_index,pe_res = pe_phi_ite%>% add_row(.before = cond_index))
  x_time3 <- Sys.time()-x_time_dummy
  tmp <- as.data.frame(do.call(rbind,x3)) %>% mutate("iteration"=Nite_i)
  names(tmp)[1:3] <- c("a","b","mu")
  x3_df <- rbind(x3_df,tmp)
  x_time <- rbind(x_time,data.frame("step2"=x_time2,"step3"=x_time3,"iteration"=Nite_i))
  }
  return(list(x1,x2_df,x3_df,x_time))
}
s <- Sys.time()
y_mod3 <- par_est_model_3(cond_index=cond_index,v=q,data_Lap = data_mod_Lap,deltal = deltal,deltau = deltau,Nite_2_3 = Nite_2_3,Nite_phi = Nite_phi)
Sys.time()-s
return(y_mod3)
}

#par_est_model_3 <- sapply(1:ncol(df_sites),FUN=model3_wrapper)
par_est_model_3 <- mclapply(1:ncol(df_sites),FUN=model3_wrapper,mc.cores=ncol(df_sites))
#par_est_model_3 <- model3_wrapper(site_i=1)

save(par_est_model_3,file="data_processed/final_model_3_parameter_estimates.RData")

# explore east coast definition
legend_text_size <- 0.7
point_size <- 0.3
legend_title_size <- 0.9
coastal_point <- function(grid) {
  sapply(1:nrow(grid),FUN = function(i) {sum(as.vector(st_distance(grid[i,],grid))<20500)<5 & grid$lon[i]>-1})
}

cp <- coastal_point(grid = xyUK20_sf) 
t1 <- tm_shape(cbind(xyUK20_sf,data.frame(cp))) + tm_dots(fill="cp",fill.scale = tm_scale_categorical(values=c("FALSE"="black","TRUE"="#C11432")),size=point_size, fill.legend = tm_legend(title="")) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,legend.show=FALSE,frame=FALSE) + tm_title(text="Longitude > -1")
t1

# try as a function of lon and lat
coastal_point <- function(grid) {
  sapply(1:nrow(grid),FUN = function(i) {sum(as.vector(st_distance(grid[i,],grid))<20500)<5 & grid$lon[i]>-2})
}

cp <- coastal_point(grid = xyUK20_sf) 
t2 <- tm_shape(cbind(xyUK20_sf,data.frame(cp))) + tm_dots(fill="cp",fill.scale = tm_scale_categorical(values=c("FALSE"="black","TRUE"="#C11432")),size=point_size, fill.legend = tm_legend(title="")) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,legend.show=FALSE,frame=FALSE) + tm_title(text="Longitude > -2")
t2

# try as a function of longitude and latitude
coastal_point <- function(grid) {
  sapply(1:nrow(grid),FUN = function(i) {sum(as.vector(st_distance(grid[i,],grid))<20500)<5 & grid$lat[i]+2*grid$lon[i]>48.5})
}

cp <- coastal_point(grid = xyUK20_sf) 
cp[df_sites[3,5]] <- FALSE
t3<- tm_shape(cbind(xyUK20_sf,data.frame(cp))) + tm_dots(fill="cp",fill.scale = tm_scale_categorical(values=c("FALSE"="black","TRUE"="#C11432")),size=point_size, fill.legend = tm_legend(title="")) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,legend.show=FALSE,frame=FALSE) + tm_title(text="East Coast?")
t3

# save as one picture
tmap_save(tmap_arrange(t1,t2,t3,ncol=3),filename="../Documents/east_coast.png",height=6,width=9)

library(tidyverse)
library(tmap)
library(tmaptools)
library(sf)
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
folder_name <- "../Documents/spatial_model_final_steps/"

# load observed data
#source("spatial_parameter_estimation.R") # for spatial_par_est function
load("data_processed/temperature_data.RData",verbose = TRUE)
load("data_processed/spatial_helper.RData", verbose = TRUE)
source("simon/P2q_function_helpers.R")
load("data_processed/P2qselected_helpers.RData", verbose = TRUE)
load(paste0("data_processed/N9000_sequential2_AGG_all12sites",q*100,".RData"),verbose = TRUE) # original parameter estimates
load("data_processed/iterative_phi0l_phi0u_estimates_London.RData",verbose=TRUE) # residual margin parameters
load("data_processed/residual_dependence_pars.RData", verbose = TRUE) # residual dependence parameters
deltal <- result_new[[1]]$deltal[1]
deltau <- result_new[[1]]$deltau[1]

source("simon/final_model_helpers.R")
# Model 3: parameter estimation ------------------------------------------------
par_est_model_3 <- function(cond_index,v=0.9,data_Lap=data_mod_Lap,grid20km=xyUK20_sf,deltal,deltau) {
  # calculate distance from the conditioning site
  dist_tmp <- as.numeric(unlist(st_distance(grid20km[cond_index,],grid20km[-cond_index,])))
  distnorm <- dist_tmp/1000000  # normalise distance using a common constant
  # subset conditioning dataset
  data_Lapv <- data_Lap %>% filter(data_Lap[,cond_index]>quantile(data_Lap[,cond_index],v))
  x1 <- x2 <- x3 <- list()
  # 1. estimate alpha and beta
  x1 <- par_est(df=data_Lapv,v = v,given = cond_index,margin = "Normal",method = "sequential2",keef_constraints = c(1,2))
  a <- discard(x1$a,is.na())
  b <- discard(x1$b,is.na())
  Z <- observed_residuals(df=data_mod_Lap,given=cond_index,v = v,a=a,b=b)
  # iterate between steps 2. and 3.
  # 2. estimate phis
  Nite_phi <- 10
  pe_res <- x1 %>% dplyr::select(sigl,sigu) %>% na.omit()
  try7 <- par_est_ite(z=Z,given=cond_index,cond_site_dist = distnorm, parest_site = pe_res,Nite = Nite_phi,show_ite=TRUE,deltal = deltal,deltau = deltau) 
  # 3. reestimate a,b,mu
  pe <- try7[[12]] %>% dplyr::select(sigl,sigu)
  # x <- sapply(1:ncol(data_Lap),FUN=NLL_AGG_wrapper,data_Lap=data_Lap,cond_index=cond_index,pe_res = pe)
  # tmp <- as.data.frame(do.call(rbind,x))
  # names(tmp) <- c("a","b","mu")
  # a <- tmp$a
  # b <- tmp$b
  # mu <- tmp$mu
#  return(list(x1,x2,x3))
  return(list(x1,try7))
}
q <- 0.9 # set quantile threshold
s <- Sys.time()
y <- par_est_model_3(cond_index=London_index,v=q,data_Lap = data_mod_Lap,deltal = deltal,deltau = deltau)
Sys.time()-s

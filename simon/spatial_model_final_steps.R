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

# 1. recreate simulated fields -----------------------------------------------
# simulate 10 fields overall --------------------------------------------------
# get index for london
x <- july3_obs[London_index]
# transform to Laplace
xL <- unif_laplace_pit(x)
# get parameter values
aest <- discard(est_all_sf %>% filter(cond_site=="London") %>% pull(a),is.na)
best <- discard(est_all_sf %>% filter(cond_site=="London") %>% pull(b),is.na)
# get residual margin parameters
pe <- as.data.frame(result_new[[12]] %>% dplyr::select(mu_agg,sigl,sigu,deltal,deltau))
names(pe) <- c("mu","sigl","sigu","deltal","deltau")
# get fields
# reconstruct the fields
y_sim <- apply(random10N,MARGIN=c(2),FUN=function(xk){xL*aest+xL^best*xk})
# add row for the conditioning site
y_sim <- as.data.frame(y_sim)
#y_sim <- y_sim %>% add_row(.before=London_index)
y_sim[London_index,] <- xL
names(y_sim) <- paste0("random",1:10)
# plot on standard Normal scale
tmpsf <- st_as_sf(cbind(y_sim,xyUK20_sf)) %>% pivot_longer(cols=contains("random"))
t <- tm_shape(tmpsf %>% mutate("name"=factor(name, levels=unique(tmpsf$name)))) + tm_dots(fill="value",size=0.5,fill.scale = tm_scale_continuous(values="-brewer.rd_bu",limits=c(-10,30)),fill.legend = tm_legend(title = "")) + tm_facets(by="name",ncol=5)
tmap_save(tm=t, filename=paste0(folder_name,"random10_y_sim_London.png"),width=10,height=8)

# transform all 10 onto 2 different scales
#unif_orig_P2q(u=july3_obs,P2q=P2q_sites2,gpdpar = gpdpar_sites2)
x_sim1 <- apply(y_sim,MARGIN=c(2),FUN=function(xk){unif_orig_P2q(u=plaplace(xk),P2q=P2q_sites2,gpdpar = gpdpar_sites2)}) %>% as.data.frame()
names(x_sim1) <- paste0("random",1:10)
# plot on standard Normal scale
tmpsf <- st_as_sf(cbind(x_sim1,xyUK20_sf)) %>% pivot_longer(cols=contains("random"))
t <- tm_shape(tmpsf %>% mutate("name"=factor(name, levels=unique(tmpsf$name)))) + tm_dots(fill="value",size=0.5,fill.scale = tm_scale_intervals(values="-brewer.rd_bu",breaks=seq(10,50,5)),fill.legend = tm_legend(title = "")) + tm_facets(by="name",ncol=5)
tmap_save(tm=t, filename=paste0(folder_name,"random10_x_sim2026.png"),width=10,height=8)
t <- tm_shape(tmpsf %>% mutate("name"=factor(name, levels=unique(tmpsf$name)))) + tm_dots(fill="value",size=0.5,fill.scale = tm_scale_continuous(values="-brewer.rd_bu",limits=c(10,50)),fill.legend = tm_legend(title = "")) + tm_facets(by="name",ncol=5)
tmap_save(tm=t, filename=paste0(folder_name,"random10_xcont_sim2026_London.png"),width=10,height=8)

x_sim2 <- apply(y_sim,MARGIN=c(2),FUN=function(xk){unif_orig_P2q(u=plaplace(xk),P2q=P2q_sites3,gpdpar = gpdpar_sites3)}) %>% as.data.frame()
names(x_sim2) <- paste0("random",1:10)
# plot on standard Normal scale
tmpsf <- st_as_sf(cbind(x_sim2,xyUK20_sf)) %>% pivot_longer(cols=contains("random"))
t <- tm_shape(tmpsf %>% mutate("name"=factor(name, levels=unique(tmpsf$name)))) + tm_dots(fill="value",size=0.5,fill.scale = tm_scale_intervals(values="-brewer.rd_bu",breaks=seq(10,50,5)),fill.legend = tm_legend(title = "")) + tm_facets(by="name",ncol=5)
tmap_save(tm=t, filename=paste0(folder_name,"random10_x_sim2076.png"),width=10,height=8)
t <- tm_shape(tmpsf %>% mutate("name"=factor(name, levels=unique(tmpsf$name)))) + tm_dots(fill="value",size=0.5,fill.scale = tm_scale_continuous(values="-brewer.rd_bu",limits=c(10,50)),fill.legend = tm_legend(title = "")) + tm_facets(by="name",ncol=5)
tmap_save(tm=t, filename=paste0(folder_name,"random10_xcont_sim2076_London.png"),width=10,height=8)

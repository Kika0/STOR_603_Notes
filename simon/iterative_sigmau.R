library(tmap) # spatial map plots
library(sf) # for handling spatial sf objects
library(viridis)
library(tidyverse)
library(latex2exp) # latex expressions for plot labels
theme_set(theme_bw())
theme_replace(
  panel.spacing = unit(2, "lines"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA) )
file.sources = list.files(pattern="*helpers.R")
sapply(file.sources,source,.GlobalEnv)
source("spatial_parameter_estimation.R") # for spatial_par_est function
load("data_processed/temperature_data.RData")
load("data_processed/spatial_helper.RData")
Birmingham <- c(-1.9032,52.4806)
Glasgow <- c(-4.258109,55.859112)
London <- c(-0.127676,51.529972)
Inverness <- c(-4.22498,57.48065) # Inverness bus station
Lancaster <- c(-2.78440,54.00871) # PSC building Lancaster University
Newcastle <- c(-1.61682,54.96902) # Newcastle railway station
Cromer <- c(1.28486,53.05349)
Hull <- c(-0.335827,53.767750)
Lowestoft <- c(1.72431,52.48435)
Truro <- c(-5.05125342465549,50.263075821232704)
Dolgellau <- c(-3.8844362867080897,52.74213275545185)
Bournemouth <- c(-1.8650607066137428,50.72173094587856)
#Leeds <- c(-1.5410242288355958,53.80098118214994)
df_sites <- data.frame(Birmingham,Glasgow,London,Inverness,Lancaster,Newcastle,Cromer,Hull,Lowestoft,Truro,Dolgellau,Bournemouth)
#spatial_par_est saves parameter estimates as est_all_sf sf object in ../Documents folder
q <- 0.9 # quantile threshold
# load all three parameter estimates sf objects
load(paste0("data_processed/N9000_sequential2_AGG_all12sites",q*100,".RData"))
est_all <- as.data.frame(est_all_sf)

# load iterative delta estimates
load("data_processed/iterative_delta_estimates.RData")

# repeat for all other sites -------------------------------------------------

Nite <- 5
tmp <- sapply(1:ncol(df_sites),iter_sigmau_site, Nite = 5, simplify = FALSE)

# separate diagnostics to allow for common phi parameters across conditioning sites ------------------------------------------------------------
phi0 <- sapply(1:ncol(df_sites),FUN = function (i) as.numeric(st_drop_geometry( tmp[[i]][1,34])))
phi1 <- sapply(1:ncol(df_sites),FUN = function (i) as.numeric(st_drop_geometry( tmp[[i]][1,35])))
# plot spatially
# get indeces of conditioning sites
get_site_index <- function(j, grid = xyUK20_sf, sites= df_sites) {
  cond_site_name <- names(sites)[j]  
  cond_site_coord <- sites %>% dplyr::select(all_of(cond_site_name)) %>% pull()
 return(find_site_index(cond_site_coord,grid_uk = grid))
}

cond_site_indeces <- sapply(1:ncol(df_sites), get_site_index)
# get the point
phi0tmp <- rep(NA,nrow(xyUK20_sf))
phi1tmp <- rep(NA, nrow(xyUK20_sf))
phi0tmp[cond_site_indeces] <- phi0
phi1tmp[cond_site_indeces] <- phi1
est_phi <- xyUK20_sf[cond_site_indeces,] %>% mutate(phi0=phi0,phi1 = phi1)

toplabel <- c(TeX("$\\phi_0$"),TeX("$\\phi_1$"))
tmphi0 <- tm_shape(xyUK20_sf) + tm_dots() + tm_shape(est_phi) + tm_dots(fill="phi0",size = 2, fill.scale =tm_scale_continuous(values="Blues"),fill.legend = tm_legend(title = TeX("$\\phi_0$"))) + tm_layout(legend.position=c("right","top"),legend.height = 12)
tmphi1 <- tm_shape(xyUK20_sf) + tm_dots() + tm_shape(est_phi) + tm_dots(fill="phi1",size = 2, fill.scale =tm_scale_continuous(values="Blues"),fill.legend = tm_legend(title = TeX("$\\phi_1$"))) + tm_layout(legend.position=c("right","top"),legend.height = 12)
t <- tmap_arrange(tmphi0,tmphi1,ncol=2)
tmap_save(t,filename=paste0("../Documents/iterative_sigmau_res_margin_all/all_phis.png"),width=8,height=6)

# plot sigma_u against distance for all sites
get_sigma_distance <- function(i, grid = xyUK20_sf,sites=df_sites) {
  sigu <- tmp[[i]]$sigu_ite_sigu
  cond_site_index <- get_site_index(i)
  dist_cond_site <- as.numeric(unlist(st_distance(grid[cond_site_index,],grid)))
  sigud <- data.frame(sigu=sigu,dist=dist_cond_site) 
  return(sigud %>% mutate(cond_site = names(sites)[i]))
}

tmp_sigu <- do.call(rbind,lapply(1:ncol(df_sites),FUN=get_sigma_distance)) 
c25 <- c(
  "dodgerblue2", "#E31A1C", # red
               "green4",
               "#6A3D9A", # purple
               "#FF7F00", # orange
               "black", "gold1",
               "skyblue2", "#FB9A99", # lt pink
               "palegreen2",
               "#CAB2D6", # lt purple
               "#FDBF6F", # lt orange
               "gray70", "khaki2",
               "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
               "darkturquoise", "green1", "yellow4", "yellow3",
               "darkorange4", "brown"
)
p <- ggplot(tmp_sigu) + geom_point(aes(y=sigu,x=dist,col=cond_site)) + scale_color_manual(values = sample(c25,ncol(df_sites))) + xlab("Distance [m]") + ylab(TeX("$\\sigma_u$")) + theme(legend.position=c("inside"),legend.position.inside = c(0.8,0.3))
ggsave(p,filename=paste0("../Documents/iterative_sigmau_res_margin/sigu_distance_all.png"),width=10,height=7) 

# save these estimates
iterative_sigmau_estimates <- tmp
save(iterative_sigmau_estimates,file="data_processed/iterative_sigmau_estimates.RData")

# consider the analysis with outliers removed ---------------------------------
x1 <- rep(0,12)
x2 <- c(rep(600000,8),800000,900000,600000,600000)
y1 <- c(0.5,1,0.5,0.8,1,1,1.2,1,1.35,1.1,0.9,1.5)
y2 <- c(3,2.2,2.5,1.5,2.2,1.5,1.1,1.05,1.25,1.8,2.5,2)
iterative_sigmau_no_outliers <- function(i) {
  # extract indeces of outliers
 index_outliers <- sigu_above_below(cond_site_name = names(df_sites)[i], x1 = x1[i], x2 = x2[i], y1 = y1[i], y2 = y2[i])
 iter_sigmau_site(i=i,index_outliers = index_outliers)
}
tmp <- sapply(1:ncol(df_sites),iterative_sigmau_no_outliers, simplify = FALSE)

# explore phi spatially again
phi0 <- sapply(1:ncol(df_sites),FUN = function (i) as.numeric(st_drop_geometry( tmp[[i]][1,34])))
phi1 <- sapply(1:ncol(df_sites),FUN = function (i) as.numeric(st_drop_geometry( tmp[[i]][1,35])))
phi0tmp <- rep(NA,nrow(xyUK20_sf))
phi1tmp <- rep(NA, nrow(xyUK20_sf))
phi0tmp[cond_site_indeces] <- phi0
phi1tmp[cond_site_indeces] <- phi1
est_phi <- xyUK20_sf[cond_site_indeces,] %>% mutate(phi0=phi0,phi1 = phi1)

toplabel <- c(TeX("$\\phi_0$"),TeX("$\\phi_1$"))
phi1_limits <- c(10,60)
tm_shape(xyUK20_sf) + tm_dots() + tm_shape(est_phi %>% pivot_longer(cols = c(phi0,phi1),names_to = "parameter", values_to = "value")) + tm_dots(fill="value",size = 2) + tm_facets(by = "parameter") + tm_layout(legend.position=c("right","top"),legend.height = 12,panels.label=toplabel)
tmphi0 <- tm_shape(xyUK20_sf) + tm_dots() + tm_shape(est_phi) + tm_dots(fill="phi0",size = 2, fill.legend = tm_legend(title = TeX("$\\phi_0$"))) + tm_layout(legend.position=c("right","top"),legend.height = 12)
tmphi1 <- tm_shape(xyUK20_sf) + tm_dots() + tm_shape(est_phi) + tm_dots(fill="phi1",size = 2, fill.legend = tm_legend(title = TeX("$\\phi_1$")),fill.scale = tm_scale_continuous(limits=phi1_limits)) + tm_layout(legend.position=c("right","top"),legend.height = 12)
t <- tmap_arrange(tmphi0,tmphi1,ncol=2)
tmap_save(t,filename=paste0("../Documents/iterative_sigmas_res_margin/all_phis.png"),width=8,height=6)

# look for maxima and minima to set for breaks
summary(est_all_sf[,11:16])
sig_u <- c()
for (i in 1:12) {
  sig_u <- append(sig_u,c(tmp[[i]]$sigu_ite_sigu,tmp[[i]]$sigu_ite,tmp[[i]]$sigu_ite-tmp[[i]]$sigu_ite_sigu))
}
min(sig_u,na.rm = TRUE)
max(sig_u,na.rm = TRUE)



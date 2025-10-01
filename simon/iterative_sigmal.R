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
Leeds <- c(-1.5410242288355958,53.80098118214994)
df_sites <- data.frame(Birmingham,Glasgow,London,Inverness,Lancaster,Newcastle,Cromer,Hull,Lowestoft,Truro,Dolgellau,Bournemouth)
#spatial_par_est saves parameter estimates as est_all_sf sf object in ../Documents folder
q <- 0.9 # quantile threshold
# load all three parameter estimates sf objects
load(paste0("data_processed/N9000_sequential2_AGG_all12sites",q*100,".RData"))
est_all <- as.data.frame(est_all_sf)

# load iterative delta estimates
load("data_processed/iterative_sigmau_estimates.R")

Nite <- 5
tmp <- sapply(1:ncol(df_sites),iter_sigmal_site, Nite = 5, simplify = FALSE)

# separate diagnostics to allow for common phi parameters across conditioning sites ------------------------------------------------------------
phi2 <- sapply(1:ncol(df_sites),FUN = function (i) as.numeric(st_drop_geometry( tmp[[i]][1,41])))
phi3 <- sapply(1:ncol(df_sites),FUN = function (i) as.numeric(st_drop_geometry( tmp[[i]][1,42])))
# plot spatially
cond_site_indeces <- sapply(1:ncol(df_sites), get_site_index)
# get the point
phi2tmp <- rep(NA,nrow(xyUK20_sf))
phi3tmp <- rep(NA, nrow(xyUK20_sf))
phi2tmp[cond_site_indeces] <- phi2
phi3tmp[cond_site_indeces] <- phi3
est_phi <- xyUK20_sf[cond_site_indeces,] %>% mutate(phi2=phi2,phi3 = phi3)

toplabel <- c(TeX("$\\phi_2$"),TeX("$\\phi_3$"))
phi2_limits <- c(0.55,1.8)
tmphi2 <- tm_shape(xyUK20_sf) + tm_dots() + tm_shape(est_phi) + tm_dots(fill="phi2",size = 2, fill.scale =tm_scale_continuous(values="Blues",limits=phi2_limits),fill.legend = tm_legend(title = TeX("$\\phi_2$"))) + tm_layout(legend.position=c("right","top"),legend.height = 12)
tmphi3 <- tm_shape(xyUK20_sf) + tm_dots() + tm_shape(est_phi) + tm_dots(fill="phi3",size = 2, fill.scale =tm_scale_continuous(values="Blues"), fill.legend = tm_legend(title = TeX("$\\phi_3$"))) + tm_layout(legend.position=c("right","top"),legend.height = 12)
t <- tmap_arrange(tmphi2,tmphi3,ncol=2)
tmap_save(t,filename=paste0("../Documents/iterative_sigmal_res_margin/all_phis.png"),width=8,height=6)

# plot sigma_u against distance for all sites
get_sigmal_distance <- function(i, grid = xyUK20_sf,sites=df_sites) {
  sigl <- tmp[[i]]$sigl_ite_sigl
  cond_site_index <- get_site_index(i)
  dist_cond_site <- as.numeric(unlist(st_distance(grid[cond_site_index,],grid)))
  sigld <- data.frame(sigl=sigl,dist=dist_cond_site) 
  return(sigld %>% mutate(cond_site = names(sites)[i]))
}

tmp_sigl <- do.call(rbind,lapply(1:ncol(df_sites),FUN=get_sigmal_distance)) 
c18 <- c(
  "#009ADA", "#C11432", # red
               "#66A64F",
               "#6A3D9A", # purple
               "#FF7F00", # orange
               "black", "#FDD10A",
               "#FB9A99", # lt pink
               "palegreen2",
               "#CAB2D6", # lt purple
               "gray70", 
               "maroon","deeppink1", "blue1",
               "darkturquoise", "green1", "yellow4",
               "darkorange4"
)
p <- ggplot(tmp_sigl) + geom_point(aes(y=sigl,x=dist,col=cond_site)) + scale_color_manual(values = sample(c18,ncol(df_sites))) + xlab("Distance [m]") + ylab(TeX("$\\sigma_l$")) + theme(legend.position=c("inside"),legend.position.inside = c(0.8,0.3))
ggsave(p,filename=paste0("../Documents/iterative_sigmal_res_margin/sigl_distance_all.png"),width=10,height=7) 

# save these estimates
iterative_sigmal_estimates <- tmp
save(iterative_sigmal_estimates,file="data_processed/iterative_sigmal_estimates.RData")

# explore plots of mu against distance
i <- 6 # Newcastle
cond_site_name <- names(df_sites)[i] 
cond_site_coord <- df_sites %>% dplyr::select(all_of(cond_site_name)) %>% pull()
cond_site <- find_site_index(cond_site_coord,grid_uk = xyUK20_sf)
tmp <- tmp[[i]]
mutmp <- data.frame(mu=tmp$mu_agg_ite_sigl,dist=as.numeric(unlist(st_distance(tmp[cond_site,],tmp))))
mud <- mutmp %>% ggplot() + geom_point(aes(y=mu,x=dist))

# plot above and below
sigl_above_below <- function(cond_site_name = "Birmingham",tmp=tmp,sites=df_sites,sigud=mutmp,x1,x2,y1,y2) {
  sigud <- sigud %>% mutate(is.above=is_above(x=dist,y=mu,x1=x1,y1=y1,x2=x2,y2=y2))
  p <- ggplot(sigud) + 
    geom_segment(x=x1,y=y1,xend=x2,yend=y2) +
    geom_point(aes(x=dist,y=mu,col=factor(is.above)),size=0.5) + 
    ylab(TeX("$\\sigma_u$")) + xlab("Distance") + scale_color_manual(values = c("black", "#C11432")) + ggtitle(cond_site_name) + guides(col="none")
  ggsave(p,filename=paste0("../Documents/iterative_sigmal_res_margin/abovebelow_sigl_distance_",cond_site_name,".png"),width=4,height=4) 
  # plot also spatially
  musf <- cbind(tmp,sigud %>% select(dist,is.above))
  t <- tm_shape(musf) + tm_dots(fill="is.above",size=0.6,fill.scale = tm_scale_categorical(values=c("TRUE" = "#C11432", "FALSE" = "black")))  + tm_title(cond_site_name) + tm_layout(legend.position=c("right","top"),legend.height = 12) 
  tmap_save(t,filename=paste0("../Documents/iterative_sigmal_res_margin/abovebelow_sigl_map_",cond_site_name,".png"),width=3,height=6) 
  return(sigud$is.above)
}

sigl_above_below(cond_site_name = "Newcastle", x1=0,y1=0,x2=600000,y2=0.8)

i <- 1 # Birmingham
cond_site_name <- names(df_sites)[i] 
cond_site_coord <- df_sites %>% dplyr::select(all_of(cond_site_name)) %>% pull()
cond_site <- find_site_index(cond_site_coord,grid_uk = xyUK20_sf)
tmp <- iterative_sigmal_estimates[[i]]
mutmp <- data.frame(mu=tmp$mu_agg_ite_sigl,dist=as.numeric(unlist(st_distance(tmp[cond_site,],tmp))))
mud <- mutmp %>% ggplot() + geom_point(aes(y=mu,x=dist))
sigl_above_below(cond_site_name = "Birmingham", x1=0,y1=0,x2=600000,y2=0.8)

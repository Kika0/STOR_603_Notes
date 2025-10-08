library(tmap)
library(sf)
library(tidyverse)
library(latex2exp)

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
load(paste0("data_processed/N9000_sequential2_AGG_all12sites",q*100,".RData"))
est_all <- as.data.frame(est_all_sf)

# identify diagonal sites
# find indeces of start and end sites
Birmingham <- c(-1.9032,52.4806)
Cromer <- c(1.28486,53.05349)
site_start <- find_site_index(site=Birmingham,grid_uk = xyUK20_sf)
site_end <- find_site_index(site=Cromer,grid_uk = xyUK20_sf)

tmap_mode("view")
uk_diag <- xyUK20_sf %>% mutate(siteID=as.numeric(1:nrow(xyUK20_sf)))
tm_shape(uk_diag) + tm_dots(fill="temp",size=1)
sites_index_diagonal <- c(192,193,194,195,196,197,218,219,220,221,242,263) # first is Birmingham and last is Cromer
site_name_diagonal <- c("Birmingham", paste0("diagonal",1:(length(sites_index_diagonal)-2)),"Cromer")
# plot these points on a map
uk_diag <- uk_diag %>% mutate(sites_diagonal=factor(case_match(siteID,c(site_start) ~ "Birmingham",c(site_end)~"Cromer",sites_index_diagonal[2:(length(site_name_diagonal)-1)]~"diagonal_sites")))
tmap_mode("plot")
t <- tm_shape(uk_diag) + tm_dots("sites_diagonal",size=0.5,fill.scale = tm_scale_categorical(values=c("Birmingham"="#C11432","Cromer" = "#009ADA", "diagonal_sites" = "#FDD10A")))
tmap_save(t, filename = "../Documents/Birmingham_Cromer_diagonal/sites_illustrate.png",width=4,height=6)
# estimate parameters along
spatial_par_est(data_Lap = data_mod_Lap,cond_sites = sites_index_diagonal,cond_site_names = site_name_diagonal,dayshift = 0,v=q,Ndays_season = 90,title = paste0("diagonal_sites_Birmingham_Cromer",q*100))

# repeat iterative procedures
# load estimates to calculate observed residuals
load("data_processed/N9000_sequential2_AGG_diagonal_sites_Birmingham_Cromer90.RData")
summary(est_all_sf)

est_all_diag1 <- est_all_sf %>% mutate(cond_site = factor(as.character(cond_site),levels=as.character(site_name_diagonal)))
tm <- map_param(tmp_est=est_all_diag1, method = "AGG", facet_var = c("cond_site"), title_map = "Sites from Birmingham to Cromer")
condmodel_params <- c("a","b","mu","sig","muagg","sigl","sigu","sigdiff","deltal","deltau","deltadiff")

save_map_i <- function(i,tm_list=tm,doc_folder = "Birmingham_Cromer_diagonal",w = 8,h = 8,par_vect = condmodel_params) {
 tmap_save(tm_list[[i]],filename = paste0("../Documents/",doc_folder,"/allsites_",par_vect[i],"_original.png"),width = w,height = h) 
}
sapply(1:length(condmodel_params),FUN = save_map_i)

# move to delta estimates
result <- sapply(1:length(sites_index_diagonal),FUN = iter_delta_site,Nite=100,sites = sites_index_diagonal,cond_site_names = site_name_diagonal,par_est = est_all_diag1,folder_name = "Birmingham_Cromer_diagonal/delta",simplify = FALSE)
#sapply(1,FUN = iter_delta_site,Nite=20,sites = sites_index_diagonal,cond_site_names = site_name_diagonal,par_est = est_all_diag1,folder_name = "Birmingham_Cromer_diagonal/delta",simplify = FALSE)

deltal <- sapply(1:length(sites_index_diagonal),FUN = function (i) as.numeric(st_drop_geometry( result[[i]][1,29])))
deltau <- sapply(1:length(sites_index_diagonal),FUN = function (i) as.numeric(st_drop_geometry( result[[i]][1,30])))
# plot spatially
deltaltmp <- rep(NA,nrow(xyUK20_sf))
deltautmp <- rep(NA, nrow(xyUK20_sf))
deltaltmp[sites_index_diagonal] <- deltal
deltautmp[sites_index_diagonal] <- deltau
est_delta <- xyUK20_sf[sites_index_diagonal,] %>% mutate(deltal=deltal,deltau = deltau)

toplabel <- c(TeX("$\\delta_l$"),TeX("$\\delta_u$"))
size_point <- 0.7
#delta_limits <- c(min(deltal,deltau),max(deltal,deltau))
deltal_limits <- c(min(deltal),max(deltal))
deltau_limits <- c(min(deltau),max(deltau))
tmdeltal <- tm_shape(xyUK20_sf) + tm_dots() + tm_shape(est_delta) + tm_dots(fill="deltal",size = size_point, fill.scale =tm_scale_continuous(values="brewer.blues",limits=deltal_limits),fill.legend = tm_legend(title = toplabel[1])) + tm_layout(legend.position=c("right","top"),legend.height = 12)
tmdeltau <- tm_shape(xyUK20_sf) + tm_dots() + tm_shape(est_delta) + tm_dots(fill="deltau",size = size_point, fill.scale =tm_scale_continuous(values="brewer.blues",limits=deltau_limits), fill.legend = tm_legend(title = toplabel[2])) + tm_layout(legend.position=c("right","top"),legend.height = 12)
t <- tmap_arrange(tmdeltal,tmdeltau,ncol=2)
tmap_save(t,filename=paste0("../Documents/Birmingham_Cromer_diagonal/all_deltas.png"),width=5,height=4)


# save the estimates to pass into iterative sigma_u
#save(result, file="data_processed/iterative_delta_estimates_Birmingham_Cromer_diagonal.RData")

load("data_processed/iterative_delta_estimates_Birmingham_Cromer_diagonal.RData")
result <- sapply(1:length(sites_index_diagonal),FUN = iter_sigmau_site,sites = sites_index_diagonal,cond_site_names = site_name_diagonal,par_est = est_all_diag1,folder_name = "Birmingham_Cromer_diagonal/sigmau",simplify = FALSE)
# save estimates
#save(result, file="data_processed/iterative_sigmau_estimates_Birmingham_Cromer_diagonal.RData")

# plot also all phis 
# separate diagnostics to allow for common phi parameters across conditioning sites ------------------------------------------------------------
phi0 <- sapply(1:length(sites_index_diagonal),FUN = function (i) as.numeric(st_drop_geometry( result[[i]][1,34])))
phi1 <- sapply(1:length(sites_index_diagonal),FUN = function (i) as.numeric(st_drop_geometry( result[[i]][1,35])))
# plot spatially
phi0tmp <- rep(NA,nrow(xyUK20_sf))
phi1tmp <- rep(NA, nrow(xyUK20_sf))
phi0tmp[sites_index_diagonal] <- phi0
phi1tmp[sites_index_diagonal] <- phi1
est_phi <- xyUK20_sf[sites_index_diagonal,] %>% mutate(phi0=phi0,phi1 = phi1)

toplabel <- c(TeX("$\\phi_0$"),TeX("$\\phi_1$"))
size_point <- 0.7
tmphi0 <- tm_shape(xyUK20_sf) + tm_dots() + tm_shape(est_phi) + tm_dots(fill="phi0",size = size_point, fill.scale =tm_scale_continuous(values="brewer.blues"),fill.legend = tm_legend(title = TeX("$\\phi_0$"))) + tm_layout(legend.position=c("right","top"),legend.height = 12)
tmphi1 <- tm_shape(xyUK20_sf) + tm_dots() + tm_shape(est_phi) + tm_dots(fill="phi1",size = size_point, fill.scale =tm_scale_continuous(values="brewer.blues"),fill.legend = tm_legend(title = TeX("$\\phi_1$"))) + tm_layout(legend.position=c("right","top"),legend.height = 12)
t <- tmap_arrange(tmphi0,tmphi1,ncol=2)
tmap_save(t,filename=paste0("../Documents/Birmingham_Cromer_diagonal/phi0_phi1.png"),width=5,height=4)

# plot sigma_u against distance for all sites
get_sigma_distance <- function(i, grid = xyUK20_sf,site_names=site_name_diagonal,tmp = result, cond_site_index = sites_index_diagonal) {
  sigu <- tmp[[i]]$sigu_ite_sigu
  dist_cond_site <- as.numeric(unlist(st_distance(grid[cond_site_index[i],],grid)))
  sigud <- data.frame(sigu=sigu,dist=dist_cond_site) 
  return(sigud %>% mutate(cond_site = site_name_diagonal[i]))
}

tmp_sigu <- do.call(rbind,lapply(1:length(site_name_diagonal),FUN=get_sigma_distance)) %>% mutate(cond_site=factor(cond_site,levels=site_name_diagonal))
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
p <- ggplot(tmp_sigu) + geom_point(aes(y=sigu,x=dist,col=cond_site)) + scale_color_manual(values = sample(c25,length(site_name_diagonal))) + xlab("Distance [m]") + ylab(TeX("$\\sigma_u$")) + theme(legend.position=c("inside"),legend.position.inside = c(0.8,0.3))
ggsave(p,filename=paste0("../Documents/Birmingham_Cromer_diagonal/sigu_distance_all.png"),width=10,height=7) 

# move to sigmal estimates
# load estimates
load("data_processed/iterative_sigmau_estimates_Birmingham_Cromer_diagonal.RData")

result <- sapply(1:length(sites_index_diagonal),FUN = iter_sigmal_site,sites = sites_index_diagonal,cond_site_names = site_name_diagonal,ite_sigu = result,folder_name = "Birmingham_Cromer_diagonal/sigmal",simplify = FALSE)
# save estimates
save(result, file="data_processed/iterative_sigmal_estimates_Birmingham_Cromer_diagonal.RData")

# examine output
phi2 <- sapply(1:length(sites_index_diagonal),FUN = function (i) as.numeric(st_drop_geometry( result[[i]][1,41])))
phi3 <- sapply(1:length(sites_index_diagonal),FUN = function (i) as.numeric(st_drop_geometry( result[[i]][1,42])))
# plot spatially
phi2tmp <- rep(NA,nrow(xyUK20_sf))
phi3tmp <- rep(NA, nrow(xyUK20_sf))
phi2tmp[sites_index_diagonal] <- phi2
phi3tmp[sites_index_diagonal] <- phi3
est_phi <- xyUK20_sf[sites_index_diagonal,] %>% mutate(phi2=phi2,phi3 = phi3)

toplabel <- c(TeX("$\\phi_2$"),TeX("$\\phi_3$"))
size_point <- 0.7
tmphi2 <- tm_shape(xyUK20_sf) + tm_dots() + tm_shape(est_phi) + tm_dots(fill="phi2",size = size_point, fill.scale =tm_scale_continuous(values="brewer.blues"),fill.legend = tm_legend(title = TeX("$\\phi_2$"))) + tm_layout(legend.position=c("right","top"),legend.height = 12)
tmphi3 <- tm_shape(xyUK20_sf) + tm_dots() + tm_shape(est_phi) + tm_dots(fill="phi3",size = size_point, fill.scale =tm_scale_continuous(values="brewer.blues"),fill.legend = tm_legend(title = TeX("$\\phi_3$"))) + tm_layout(legend.position=c("right","top"),legend.height = 12)
t <- tmap_arrange(tmphi2,tmphi3,ncol=2)
tmap_save(t,filename=paste0("../Documents/Birmingham_Cromer_diagonal/phi2_phi3.png"),width=5,height=4)

# plot sigma_u against distance for all sites
get_sigma_distance <- function(i, grid = xyUK20_sf,site_names=site_name_diagonal,tmp = result, cond_site_index = sites_index_diagonal) {
  sigl <- tmp[[i]]$sigl_ite_sigl
  dist_cond_site <- as.numeric(unlist(st_distance(grid[cond_site_index[i],],grid)))
  sigld <- data.frame(sigl=sigl,dist=dist_cond_site) 
  return(sigld %>% mutate(cond_site = site_name_diagonal[i]))
}

tmp_sigl <- do.call(rbind,lapply(1:length(site_name_diagonal),FUN=get_sigma_distance)) %>% mutate(cond_site=factor(cond_site,levels=site_name_diagonal))
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
p <- ggplot(tmp_sigl) + geom_point(aes(y=sigl,x=dist,col=cond_site)) + scale_color_manual(values = sample(c25,length(site_name_diagonal))) + xlab("Distance [m]") + ylab(TeX("$\\sigma_l$")) + theme(legend.position=c("inside"),legend.position.inside = c(0.8,0.3))
ggsave(p,filename=paste0("../Documents/Birmingham_Cromer_diagonal/sigl_distance_all.png"),width=10,height=7) 

# try phi0,phi1,phi2,phi3 joint estimation
load("data_processed/iterative_delta_estimates_Birmingham_Cromer_diagonal.RData")
result <- sapply(1:length(sites_index_diagonal),FUN = iter_sigmas_site,sites = sites_index_diagonal,cond_site_names = site_name_diagonal,par_est = est_all_diag1,folder_name = "Birmingham_Cromer_diagonal/sigmas",simplify = FALSE)
# save these estimates
#save(result, file="data_processed/iterative_sigmas_estimates_Birmingham_Cromer_diagonal.RData")
load("data_processed/iterative_sigmas_estimates_Birmingham_Cromer_diagonal.RData")

# examine output
phi0 <- sapply(1:length(sites_index_diagonal),FUN = function (i) as.numeric(st_drop_geometry( result[[i]][1,34])))
phi1 <- sapply(1:length(sites_index_diagonal),FUN = function (i) as.numeric(st_drop_geometry( result[[i]][1,35])))
phi2 <- sapply(1:length(sites_index_diagonal),FUN = function (i) as.numeric(st_drop_geometry( result[[i]][1,36])))
phi3 <- sapply(1:length(sites_index_diagonal),FUN = function (i) as.numeric(st_drop_geometry( result[[i]][1,37])))

# plot spatially
phi0tmp <- phi1tmp <- rep(NA,nrow(xyUK20_sf))
phi2tmp <- phi3tmp <- rep(NA, nrow(xyUK20_sf))
phi0tmp[sites_index_diagonal] <- phi0
phi1tmp[sites_index_diagonal] <- phi1
phi2tmp[sites_index_diagonal] <- phi2
phi3tmp[sites_index_diagonal] <- phi3
est_phi <- xyUK20_sf[sites_index_diagonal,] %>% mutate(phi0=phi0,phi1=phi1,phi2=phi2,phi3 = phi3)

toplabel <- c(TeX("$\\phi_0$"),TeX("$\\phi_1$"),TeX("$\\phi_2$"),TeX("$\\phi_3$"))
size_point <- 0.7
tmphi0 <- tm_shape(xyUK20_sf) + tm_dots() + tm_shape(est_phi) + tm_dots(fill="phi0",size = size_point, fill.scale =tm_scale_continuous(values="brewer.blues"),fill.legend = tm_legend(title = TeX("$\\phi_0$"))) + tm_layout(legend.position=c("right","top"),legend.height = 12)
tmphi1 <- tm_shape(xyUK20_sf) + tm_dots() + tm_shape(est_phi) + tm_dots(fill="phi1",size = size_point, fill.scale =tm_scale_continuous(values="brewer.blues"),fill.legend = tm_legend(title = TeX("$\\phi_1$"))) + tm_layout(legend.position=c("right","top"),legend.height = 12)
tmphi2 <- tm_shape(xyUK20_sf) + tm_dots() + tm_shape(est_phi) + tm_dots(fill="phi2",size = size_point, fill.scale =tm_scale_continuous(values="brewer.blues"),fill.legend = tm_legend(title = TeX("$\\phi_2$"))) + tm_layout(legend.position=c("right","top"),legend.height = 12)
tmphi3 <- tm_shape(xyUK20_sf) + tm_dots() + tm_shape(est_phi) + tm_dots(fill="phi3",size = size_point, fill.scale =tm_scale_continuous(values="brewer.blues"),fill.legend = tm_legend(title = TeX("$\\phi_3$"))) + tm_layout(legend.position=c("right","top"),legend.height = 12)
t <- tmap_arrange(tmphi0,tmphi1,tmphi2,tmphi3,ncol=4)
tmap_save(t,filename=paste0("../Documents/Birmingham_Cromer_diagonal/phis_all4.png"),width=8,height=4)

# plots against 94 marginal quantile
load("data_processed/94thresholds.RData")
temp94 <- tm_thres %>% filter(data_source=="CPM") %>% pull(temp)
temp94diff <- tm_thres_diff %>% filter(data_source=="CPM_diff") %>% pull(temperature)
data.frame(phi0,temp94=temp94[sites_index_diagonal]) %>%  ggplot() + geom_point(aes(x=temp94,y=phi0))
data.frame(phi0,temp94diff=temp94diff[sites_index_diagonal]) %>%  ggplot() + geom_point(aes(x=temp94diff,y=phi0))

# combine into one plot
tmp <- data.frame(param_value=c(phi0,phi1,phi2,phi3),temp94=rep(temp94[sites_index_diagonal],4),phi = rep(c("phi0","phi1","phi2","phi3"),each=4))
p <- ggplot(tmp) + geom_point(aes(x=temp94,y=param_value,col=phi)) + geom_line(aes(x=temp94,y=param_value,col=phi))
ggsave(p,filename="../Documents/phis_temp94.png",width=5,height=5)
tmp <- data.frame(param_value=c(phi0,phi1,phi2,phi3),temp94diff=rep(temp94diff[sites_index_diagonal],4),phi = rep(c("phi0","phi1","phi2","phi3"),each=4))
p <- ggplot(tmp) + geom_point(aes(x=temp94diff,y=param_value,col=phi)) + geom_line(aes(x=temp94diff,y=param_value,col=phi))
ggsave(p,filename="../Documents/phis_temp94diff.png",width=5,height=5)

# iterative a,b,mu and 
load("data_processed/iterative_sigmal_estimates_Birmingham_Cromer_diagonal.RData")
sites = sites_index_diagonal
cond_site_names = site_name_diagonal
q <- 0.9
if(is.null(cond_site_names)) {
  cond_site_name <- names(sites)[j]
  cond_site_names <- names(sites)
} else {  cond_site_name <- cond_site_names[j] }

if (is.numeric(sites)) {cond_site <- sites[j]} else{
  cond_site_coord <- sites %>% dplyr::select(all_of(cond_site_name)) %>% pull()
  cond_site <- find_site_index(cond_site_coord,grid_uk = grid)    }

# calculate distance from the conditioning site
dist_tmp <- as.numeric(unlist(st_distance(xyUK20_sf[cond_site,],xyUK20_sf)))
# remove zero distance
dist_tmp <- dist_tmp[dist_tmp>0]
# normalise distance using a common constant
distnorm <- dist_tmp/1000000
parest_site <- st_drop_geometry(result[[j]]) %>% dplyr::select(sigl_ite_sigl,sigu_ite_sigu,deltal_ite,deltau_ite) %>% na.omit()

try7 <- par_est_ite(dataLap=data_mod_Lap,given=cond_site,cond_site_dist=distnorm, parest_site = parest_site,Nite=10, show_ite=TRUE)

# separate parameter estimation and analysis


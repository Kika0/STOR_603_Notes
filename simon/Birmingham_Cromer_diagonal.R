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
result <- sapply(1:length(sites_index_diagonal),FUN = iter_delta_site,Nite=50,sites = sites_index_diagonal,cond_site_names = site_name_diagonal,par_est = est_all_diag1,folder_name = "Birmingham_Cromer_diagonal/delta",simplify = FALSE)
# save the estimates to pass into iterative sigma_u
save(result, file="data_processed/iterative_delta_estimates_Birmingham_Cromer_diagonal.RData")

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
c12 <- c(
  "#009ADA", "#C11432", # red
               "green4",
               "#6A3D9A", # purple
               "#FF7F00", # orange
               "black", "gold1",
               "#FB9A99", # lt pink
               "gray70", 
               "darkturquoise", "green1", 
               "darkorange4"
)
p <- ggplot(tmp_sigu) + geom_point(aes(y=sigu,x=dist,col=cond_site)) + scale_color_manual(values = c12) + xlab("Distance [m]") + ylab(TeX("$\\sigma_u$")) + theme(legend.position=c("inside"),legend.position.inside = c(0.8,0.3))
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
  sigl <- result[[i]]$sigl_ite_sigl
  dist_cond_site <- as.numeric(unlist(st_distance(grid[cond_site_index[i],],grid)))
  sigld <- data.frame(sigl=sigl,dist=dist_cond_site) 
  return(sigld %>% mutate(cond_site = site_name_diagonal[i]))
}

tmp_sigl <- do.call(rbind,lapply(1:length(site_name_diagonal),FUN=get_sigma_distance)) %>% mutate(cond_site=factor(cond_site,levels=site_name_diagonal))
p <- ggplot(tmp_sigl) + geom_point(aes(y=sigl,x=dist,col=cond_site)) + scale_color_manual(values = c12) + xlab("Distance [m]") + ylab(TeX("$\\sigma_l$")) + theme(legend.position=c("inside"),legend.position.inside = c(0.8,0.3))
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

# plots against 94 marginal quantile --------
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

# iterative a,b,mu and phis --------------------------------------------------
load("data_processed/iterative_sigmal_estimates_Birmingham_Cromer_diagonal.RData")
abmu_par_est_ite <- function(site,Nite=10,sites = sites_index_diagonal,cond_site_names = site_name_diagonal,q=0.9,grid=xyUK20_sf,result,est_all_sf,deltal=NULL,deltau=NULL) {
if(is.null(cond_site_names)) {
  cond_site_name <- names(sites)[site]
  cond_site_names <- names(sites)
} else {  cond_site_name <- cond_site_names[site] }

if (is.numeric(sites)) {cond_site <- sites[site]} else{
  cond_site_coord <- sites %>% dplyr::select(all_of(cond_site_name)) %>% pull()
  cond_site <- find_site_index(cond_site_coord,grid_uk = grid)    }

# calculate distance from the conditioning site
dist_tmp <- as.numeric(unlist(st_distance(xyUK20_sf[cond_site,],xyUK20_sf)))
# remove zero distance (conditioning site)
dist_tmp <- dist_tmp[dist_tmp>0]
# normalise distance using a common constant
distnorm <- dist_tmp/1000000
parest_site <- st_drop_geometry(result[[site]]) %>% dplyr::select(sigl_ite_sigl,sigu_ite_sigu,deltal_ite,deltau_ite) %>% na.omit()
if (is.null(deltal)) {
  try7 <- par_est_ite(dataLap=data_mod_Lap,given=cond_site,cond_site_dist=distnorm, parest_site = parest_site,Nite=Nite, show_ite=TRUE)
} else {
try7 <- par_est_ite(dataLap=data_mod_Lap,given=cond_site,cond_site_dist=distnorm, parest_site = parest_site,Nite=Nite, show_ite=TRUE,deltal=deltal,deltau=deltau) }
# print summary
sapply(1:12,function(i)print(summary(try7[[i]])),simplify=FALSE)
# separate parameter estimation and analysis
p <- ggplot(try7[[1]] %>% pivot_longer(everything(),names_to = "iteration",values_to = "par") %>% mutate(iteration=factor(iteration,levels=paste0("X",1:Nite)))) + geom_boxplot(aes(x=iteration,y=par)) +ylab(TeX("$\\alpha"))
ggsave(p,file=paste0("../Documents/abmu_iterative/boxplot_alpha_",cond_site_name,".png"),height=5,width=10)
p <- ggplot(try7[[2]] %>% pivot_longer(everything(),names_to = "iteration",values_to = "par")%>% mutate(iteration=factor(iteration,levels=paste0("X",1:Nite)))) + geom_boxplot(aes(x=iteration,y=par)) +ylab(TeX("$\\beta"))
ggsave(p,file=paste0("../Documents/abmu_iterative/boxplot_beta_",cond_site_name,".png"),height=5,width=10)
p <- ggplot(try7[[3]] %>% pivot_longer(everything(),names_to = "iteration",values_to = "par")%>% mutate(iteration=factor(iteration,levels=paste0("X",1:Nite)))) + geom_boxplot(aes(x=iteration,y=par)) +ylab(TeX("$\\mu_{AGG}"))
ggsave(p,file=paste0("../Documents/abmu_iterative/boxplot_mu_",cond_site_name,".png"),height=5,width=10)
p <- ggplot(try7[[4]] %>% pivot_longer(everything(),names_to = "iteration",values_to = "par")%>% mutate(iteration=factor(iteration,levels=paste0("X",1:Nite)))) + geom_boxplot(aes(x=iteration,y=par)) +ylab(TeX("$\\sigma_l"))
ggsave(p,file=paste0("../Documents/abmu_iterative/boxplot_sigmal_",cond_site_name,".png"),height=5,width=10)
p <- ggplot(try7[[5]] %>% pivot_longer(everything(),names_to = "iteration",values_to = "par")%>% mutate(iteration=factor(iteration,levels=paste0("X",1:Nite)))) + geom_boxplot(aes(x=iteration,y=par)) +ylab(TeX("$\\sigma_u"))
ggsave(p,file=paste0("../Documents/abmu_iterative/boxplot_sigmau_",cond_site_name,".png"),height=5,width=10)

p <- ggplot(data.frame("deltal"=try7[[10]][2:(Nite+1)],"iteration"=1:Nite)) + geom_point(aes(x=factor(iteration),y=deltal),size=1.5) +ylab(TeX("$\\delta_l"))
ggsave(p,file=paste0("../Documents/abmu_iterative/plot_deltal_",cond_site_name,".png"),height=5,width=10)

p <- ggplot(data.frame("deltau"=try7[[11]][2:(Nite+1)],"iteration"=1:Nite)) + geom_point(aes(x=factor(iteration),y=deltau),size=1.5) +ylab(TeX("$\\delta_u"))
ggsave(p,file=paste0("../Documents/abmu_iterative/plot_deltau_",cond_site_name,".png"),height=5,width=10)

# look spatially to check
# explore also spatial parameters
est_ite <- try7[[12]] %>% add_row(.before=cond_site)
names(est_ite) <- paste0(names(est_ite),"_ite")
tmpsf <- cbind(est_all_sf %>% filter(cond_site==cond_site_name),est_ite)
# plot parameter estimates against distance
# calculate distance from a conditioning site with st_distance()
folder_name <- "abmu_iterative"
mud <- data.frame(mu=tmpsf$mu_agg_ite,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=mu,x=dist))
ggsave(mud,filename=paste0("../Documents/",folder_name,"/muagg_distance_",cond_site_name,".png")) 

sigld <- data.frame(sigl=tmpsf$sigl_ite,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=sigl,x=dist))
ggsave(sigld,filename=paste0("../Documents/",folder_name,"/sigl_distance_",cond_site_name,".png")) 

sigud <- data.frame(sigu=tmpsf$sigu_ite,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=sigu,x=dist))
ggsave(sigud,filename=paste0("../Documents/",folder_name,"/sigu_distance_",cond_site_name,".png")) 


tmpsf <- tmpsf %>% mutate(adiff=a_ite-a,bdiff=b_ite-b,mudiff = mu_agg_ite - mu_agg, sigldiff = sigl_ite - sigl, sigudiff = sigu_ite - sigu)
toplabel <- c("New iterative approach","Original method","Difference")
t <- tmpsf %>% dplyr::select(a_ite,a,adiff) %>% pivot_longer(cols=c(a_ite,a,adiff),names_to = "parameter", values_to = "value" ) %>% mutate(parameter=factor(parameter,levels=c("a_ite","a","adiff"))) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="-brewer.rd_bu"),fill.legend = tm_legend(title = TeX("$\\alpha$"))) + tm_facets("parameter") +
  tm_layout(legend.position=c("right","top"),legend.height = 12, panel.labels = toplabel,legend.reverse = TRUE) 

tmap_save(t,filename=paste0("../Documents/",folder_name,"/a_",cond_site_name,".png"),width=8,height=6)

t <- tmpsf %>% dplyr::select(b_ite,b,bdiff) %>% pivot_longer(cols=c(b_ite,b,bdiff),names_to = "parameter", values_to = "value" ) %>% mutate(parameter=factor(parameter,levels=c("b_ite","b","bdiff"))) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="-brewer.rd_bu"),fill.legend = tm_legend(title = TeX("$\\beta$"))) + tm_facets("parameter") +
  tm_layout(legend.position=c("right","top"),legend.height = 12, panel.labels = toplabel,legend.reverse = TRUE) 

tmap_save(t,filename=paste0("../Documents/",folder_name,"/b_",cond_site_name,".png"),width=8,height=6)

#mu_limits <- c(-2.61,1)
t <- tmpsf %>% dplyr::select(mu_agg_ite,mu_agg,mudiff) %>% pivot_longer(cols=c(mu_agg_ite,mu_agg,mudiff),names_to = "parameter", values_to = "value" ) %>% mutate(parameter=factor(parameter,levels=c("mu_agg_ite","mu_agg","mudiff"))) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="-brewer.rd_bu"),fill.legend = tm_legend(title = TeX("$\\mu_{AGG}$"))) + tm_facets("parameter") +
  tm_layout(legend.position=c("right","top"),legend.height = 12, panel.labels = toplabel,legend.reverse = TRUE) 

tmap_save(t,filename=paste0("../Documents/",folder_name,"/mu_agg_",cond_site_name,".png"),width=8,height=6)

#sigma_limits <- c(-1.32,2.3)
t <- tmpsf %>% dplyr::select(sigl_ite,sigl,sigldiff) %>% pivot_longer(cols=c(sigl_ite,sigl,sigldiff),names_to = "parameter", values_to = "value" ) %>% mutate(parameter=factor(parameter,levels=c("sigl_ite","sigl","sigldiff"))) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="-brewer.rd_bu"),fill.legend = tm_legend(title = TeX("$\\sigma_l$"))) + tm_facets("parameter") +
  tm_layout(legend.position=c("right","top"),legend.height = 9, panel.labels = toplabel,legend.reverse = TRUE) 
tmap_save(t,filename=paste0("../Documents/",folder_name,"/sigma_lower_",cond_site_name,".png"),width=8,height=6)

t <- tmpsf %>% dplyr::select(sigu_ite,sigu,sigudiff) %>% pivot_longer(cols=c(sigu_ite,sigu,sigudiff),names_to = "parameter", values_to = "value" ) %>% mutate(parameter=factor(parameter,levels=c("sigu_ite","sigu","sigudiff"))) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="-brewer.rd_bu"),fill.legend = tm_legend(title = TeX("$\\sigma_u$"))) + tm_facets("parameter") +
  tm_layout(legend.position=c("right","top"),legend.height = 9, panel.labels = toplabel,legend.reverse = TRUE) 
tmap_save(t,filename=paste0("../Documents/",folder_name,"/sigma_upper_",cond_site_name,".png"),width=8,height=6)
return(try7)
}

tmp <- sapply(1,FUN=abmu_par_est_ite,result=result,est_all_sf=est_all_sf,simplify=FALSE)
#tmp <- sapply(1:length(sites_index_diagonal),FUN=abmu_par_est_ite,result=result,est_all_sf=est_all_sf,simplify=FALSE)
#save(tmp,file="data_processed/Birmingham_Cromer_abmu_iterative.RData")
# plot estimates
load("data_processed/Birmingham_Cromer_abmu_iterative.RData")

# examine output
phi0 <- sapply(1:length(sites_index_diagonal),FUN = function (i) as.numeric(st_drop_geometry( tmp[[i]][[12]][1,6])))
phi1 <- sapply(1:length(sites_index_diagonal),FUN = function (i) as.numeric(st_drop_geometry( tmp[[i]][[12]][1,7])))
phi2 <- sapply(1:length(sites_index_diagonal),FUN = function (i) as.numeric(st_drop_geometry( tmp[[i]][[12]][1,8])))
phi3 <- sapply(1:length(sites_index_diagonal),FUN = function (i) as.numeric(st_drop_geometry( tmp[[i]][[12]][1,9])))

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
tmap_save(t,filename=paste0("../Documents/Birmingham_Cromer_diagonal/phis_all4_abmu.png"),width=8,height=4)

# plots against 94 marginal quantile
load("data_processed/94thresholds.RData")
temp94 <- tm_thres %>% filter(data_source=="CPM") %>% pull(temp)
temp94diff <- tm_thres_diff %>% filter(data_source=="CPM_diff") %>% pull(temperature)
data.frame(phi0,temp94=temp94[sites_index_diagonal]) %>%  ggplot() + geom_point(aes(x=temp94,y=phi0))
data.frame(phi0,temp94diff=temp94diff[sites_index_diagonal]) %>%  ggplot() + geom_point(aes(x=temp94diff,y=phi0))

# combine into one plot
tmp <- data.frame(param_value=c(phi0,phi1,phi2,phi3),temp94=rep(temp94[sites_index_diagonal],4),phi = rep(c("phi0","phi1","phi2","phi3"),each=4))
p <- ggplot(tmp) + geom_point(aes(x=temp94,y=param_value,col=phi)) + geom_line(aes(x=temp94,y=param_value,col=phi))
ggsave(p,filename="../Documents/phis_temp94_abmu.png",width=5,height=5)
tmp <- data.frame(param_value=c(phi0,phi1,phi2,phi3),temp94diff=rep(temp94diff[sites_index_diagonal],4),phi = rep(c("phi0","phi1","phi2","phi3"),each=4))
p <- ggplot(tmp) + geom_point(aes(x=temp94diff,y=param_value,col=phi)) + geom_line(aes(x=temp94diff,y=param_value,col=phi))
ggsave(p,filename="../Documents/phis_temp94diff_abmu.png",width=5,height=5)

deltal <- sapply(1:length(sites_index_diagonal),FUN = function (i) as.numeric(st_drop_geometry( tmp[[i]][[12]][1,10])))
deltau <- sapply(1:length(sites_index_diagonal),FUN = function (i) as.numeric(st_drop_geometry( tmp[[i]][[12]][1,11])))
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
tmap_save(t,filename=paste0("../Documents/Birmingham_Cromer_diagonal/all_deltas_abmu.png"),width=5,height=4)


# plot sigma_u against distance for all sites
get_sigma_distance <- function(i, grid = xyUK20_sf,site_names=site_name_diagonal,tmp = result, cond_site_index = sites_index_diagonal) {
  sigl <- tmp[[i]]$sigl_ite_sigl
  dist_cond_site <- as.numeric(unlist(st_distance(grid[cond_site_index[i],],grid)))
  sigld <- data.frame(sigl=sigl,dist=dist_cond_site) 
  return(sigld %>% mutate(cond_site = site_name_diagonal[i]))
}

tmp_sigl <- do.call(rbind,lapply(1:length(site_name_diagonal),FUN=get_sigma_distance)) %>% mutate(cond_site=factor(cond_site,levels=site_name_diagonal))
p <- ggplot(tmp_sigl) + geom_point(aes(y=sigl,x=dist,col=cond_site)) + scale_color_manual(values = c12) + xlab("Distance [m]") + ylab(TeX("$\\sigma_l$")) + theme(legend.position=c("inside"),legend.position.inside = c(0.8,0.3))
ggsave(p,filename=paste0("../Documents/Birmingham_Cromer_diagonal/sigl_distance_all.png"),width=10,height=7) 

# try deltas fixed as mean of 12 deltas
# load(file="data_processed/iterative_delta_estimates_Birmingham_Cromer_diagonal.RData")
# deltal <- sapply(1:length(sites_index_diagonal),FUN = function (i) as.numeric(st_drop_geometry( result[[i]][1,29])))
# deltau <- sapply(1:length(sites_index_diagonal),FUN = function (i) as.numeric(st_drop_geometry( result[[i]][1,30])))
load("data_processed/iterative_sigmal_estimates_Birmingham_Cromer_diagonal.RData")
deltal <- sapply(1:length(sites_index_diagonal),FUN = function (i) as.numeric(st_drop_geometry( result[[i]][1,29])))
deltau <- sapply(1:length(sites_index_diagonal),FUN = function (i) as.numeric(st_drop_geometry( result[[i]][1,30])))

#tmp <- sapply(1:length(sites_index_diagonal),FUN=function(k) {par_est_ite(dataLap = data_mod_Lap,given = sites_index_diagonal[k],parest_site=result[[k]],deltal=mean(deltal),deltau=mean(deltau))},simplify=FALSE)
#tmp <- sapply(1:length(sites_index_diagonal),FUN=abmu_par_est_ite,result=result,est_all_sf=est_all_sf,simplify=FALSE)
tmp <- sapply(2,FUN=abmu_par_est_ite,result=result,est_all_sf=est_all_sf,deltal=mean(deltal),deltau=mean(deltau),simplify=FALSE)

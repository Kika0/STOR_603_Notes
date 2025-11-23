library(tmap)
library(sf)
library(tidyverse)
library(latex2exp)
library(gridExtra)

theme_set(theme_bw())
theme_replace(
  panel.spacing = unit(2, "lines"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA) )

load("data_processed/spatial_helper.RData") # xyUK20_sf 20km grid
file.sources = list.files(pattern="*helpers.R")
sapply(file.sources,source,.GlobalEnv)

# set folder name to save plots to
folder_name <- "../Documents/GLEN_workshop_plots/"

# 1. plot alphas for a couple of sites ------------------------------
# load estimates
q <- 0.9
load(paste0("data_processed/N9000_sequential2_AGG_all12sites",q*100,".RData"))
#est_all <- as.data.frame(est_all_sf)
# plot alpha values for Lancaster, Birmingham, Cromer
cond_site_names <- c("Lancaster","Birmingham","Cromer")
est_sites <- est_all_sf %>% filter(cond_site %in% cond_site_names) %>% mutate(cond_site=factor(cond_site,levels = cond_site_names))
title_map <- ""
misscol <- "aquamarine"
legend_text_size <- 0.7
point_size <- 0.6
legend_title_size <- 1.2
lims <- c(0,1)
nrow_facet <- 1
p <- tm_shape(est_sites) + tm_dots(fill="a",fill.scale = tm_scale_continuous(limits=lims,values="brewer.blues",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\alpha$"))) + tm_facets(by="cond_site",nrow = nrow_facet) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE) + tm_title(text=title_map) 
tmap_save(p,filename=paste0(folder_name,"alpha_selected_sites.png"))
lims <- c(0,0.5)
p <- tm_shape(est_sites) + tm_dots(fill="b",fill.scale = tm_scale_continuous(limits=lims,values="brewer.blues",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\beta$"))) + tm_facets(by="cond_site",nrow = nrow_facet) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE) + tm_title(text=title_map) 
tmap_save(p,filename=paste0(folder_name,"beta_selected_sites.png"))
# remove objects
rm(est_sites,p)

# 2. spatial and temporal estimates from Birmingham to London
# calculate estimates from Birmingham to London
# identify diagonal sites
# find indeces of start and end sites
source("spatial_parameter_estimation.R") # for spatial_par_est function
load("data_processed/temperature_data.RData") # for data_mod and data_mod_Lap
q <- 0.9
Birmingham <- c(-1.9032,52.4806)
London <- c(-0.127676,51.529972)
site_start <- find_site_index(site=Birmingham,grid_uk = xyUK20_sf)
site_end <- find_site_index(site=London,grid_uk = xyUK20_sf)
sites_index_diagonal <- c(192,193,175,157,137,117,99,100) # first is Birmingham and last is London
site_name_diagonal <- c("Birmingham", paste0("diagonal",1:(length(sites_index_diagonal)-2)),"London")
spatial_par_est(data_Lap = data_mod_Lap,cond_sites = sites_index_diagonal,cond_site_names = site_name_diagonal,dayshift = c(-3:3),v=q,Ndays_season = 90,title = paste0("diagonal_sites_Birmingham_London",q*100))


# 4. estimates of iterative approach compared with original ------------
# show these for Birmingham and Cromer
load("data_processed/iterative_abmu_fixed_delta.RData")
# identify diagonal sites
# find indeces of start and end sites
Birmingham <- c(-1.9032,52.4806)
Cromer <- c(1.28486,53.05349)
site_start <- find_site_index(site=Birmingham,grid_uk = xyUK20_sf)
site_end <- find_site_index(site=Cromer,grid_uk = xyUK20_sf)
sites_index_diagonal <- c(192,193,194,195,196,197,218,219,220,221,242,263) # first is Birmingham and last is Cromer
site_name_diagonal <- c("Birmingham", paste0("diagonal",1:(length(sites_index_diagonal)-2)),"Cromer")
# load also original estimates
load("data_processed/N9000_sequential2_AGG_diagonal_sites_Birmingham_Cromer90.RData")

map_ab_iterative_difference <- function(site,cond_site_names=site_name_diagonal,cond_site_indeces=sites_index_diagonal,estimates=tmp_fixed_deltas) {
legend_text_size <- 0.7
point_size <- 0.6
legend_title_size <- 1.2
misscol <- "aquamarine" 
cond_site <- cond_site_indeces[site]
cond_site_name <- cond_site_names[site]
est_ite <- estimates[[site]][[12]] %>% add_row(.before=cond_site)
names(est_ite) <- paste0(names(est_ite),"_ite")
tmpsf <- cbind(est_all_sf %>% filter(cond_site==cond_site_name),est_ite)

tmpsf <- tmpsf %>% mutate(adiff=a_ite-a,bdiff=b_ite-b,mudiff = mu_agg_ite - mu_agg, sigldiff = sigl_ite - sigl, sigudiff = sigu_ite - sigu)
toplabel <- c("New iterative approach","Original method","Difference")
tmpa <- tmpsf %>% dplyr::select(a_ite,a,adiff) %>% pivot_longer(cols=c(a_ite,a,adiff),names_to = "parameter", values_to = "value" ) %>% mutate(parameter=factor(parameter,levels=c("a_ite","a","adiff"))) 
limits <- c(0,1)
t1 <- tmpa %>% filter(parameter=="a_ite") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="brewer.blues",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\alpha$"))) + tm_facets("parameter") +
  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[1],legend.reverse = TRUE) 
t2 <- tmpa %>% filter(parameter=="a") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="brewer.blues",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\alpha$"))) + tm_facets("parameter") +
  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[2],legend.reverse = TRUE) 
limits <- c(-0.3,0.3)
t3 <- tmpa %>% filter(parameter=="adiff") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="-brewer.rd_bu",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\alpha$"))) + tm_facets("parameter") +
  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[3],legend.reverse = TRUE) 
t <- tmap_arrange(t1,t2,t3,ncol=3)
tmap_save(t,filename=paste0(folder_name,"/a_",cond_site_name,"_iterative_original_difference.png"),width=8,height=6)

tmpb <- tmpsf %>% dplyr::select(b_ite,b,bdiff) %>% pivot_longer(cols=c(b_ite,b,bdiff),names_to = "parameter", values_to = "value" ) %>% mutate(parameter=factor(parameter,levels=c("b_ite","b","bdiff"))) 
limits <- c(0,0.5)
t1 <- tmpb %>% filter(parameter=="b_ite") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="brewer.blues",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\beta$"))) + tm_facets("parameter") +
  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[1],legend.reverse = TRUE) 
t2 <- tmpb %>% filter(parameter=="b") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="brewer.blues",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\beta$"))) + tm_facets("parameter") +
  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[2],legend.reverse = TRUE) 
limits <- c(-0.5,0.5)
t3 <- tmpb %>% filter(parameter=="bdiff") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="-brewer.rd_bu",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\beta$"))) + tm_facets("parameter") +
  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[3],legend.reverse = TRUE) 
t <- tmap_arrange(t1,t2,t3,ncol=3)
tmap_save(t,filename=paste0(folder_name,"/b_",cond_site_name,"_iterative_original_difference.png"),width=8,height=6)
}

map_ab_iterative_difference(site=1)
map_ab_iterative_difference(site=12)

# 5. plot of sigma values comparison -----------------------------------
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
# plot sigmas against distance
get_sigma_distance <- function(i, grid = xyUK20_sf,site_names=site_name_diagonal,tmp = tmp_fixed_deltas,sig_which="sigu", cond_site_index = sites_index_diagonal) {
  sig <- tmp[[i]][[12]] %>% dplyr::select(contains(sig_which)) %>% pull()
  sig <- append(sig,NA,after=(cond_site_index[i]-1))
  dist_cond_site <- as.numeric(unlist(st_distance(grid[cond_site_index[i],],grid) %>%
                                        units::set_units(mi)))
  sigd <- data.frame(sig=sig,dist=dist_cond_site) 
  return(sigd %>% mutate(cond_site = site_name_diagonal[i]))
}

lims <- c(0,2.3)
point_size <- 0.7
tmp_sigl <- do.call(rbind,lapply(1:length(site_name_diagonal),FUN=get_sigma_distance,sig_which="sigl")) %>% mutate(cond_site=factor(cond_site,levels=site_name_diagonal))
p1 <- ggplot(tmp_sigl) + geom_point(aes(y=sig,x=dist,col=cond_site),size=point_size) + scale_color_manual(values = c12) + xlab("Distance [mi]") + ylab(TeX("$\\sigma_l$")) + theme(legend.position="none") + ylim(lims)
tmp_sigu <- do.call(rbind,lapply(1:length(site_name_diagonal),FUN=get_sigma_distance,sig_which="sigu")) %>% mutate(cond_site=factor(cond_site,levels=site_name_diagonal))
p2 <- ggplot(tmp_sigu) + geom_point(aes(y=sig,x=dist,col=cond_site),size=point_size) + scale_color_manual(values = c12) + xlab("Distance [mi]") + ylab(TeX("$\\sigma_u$")) + theme(legend.position="none") + ylim(lims)
tmp_sigmas <- data.frame("sigl"=tmp_sigl$sig, "sigu" = tmp_sigu$sig, "cond_site" = tmp_sigl$cond_site)
p3 <- ggplot(tmp_sigmas) + geom_point(aes(y=sigu,x=sigl,col=cond_site),size=point_size) + scale_color_manual(values = c12) + xlab(TeX("$\\sigma_l$")) + ylab(TeX("$\\sigma_u$")) + xlim(lims) + ylim(lims) + coord_fixed() + labs(col="Conditioning site") + guides(col=guide_legend(override.aes=list(size=4))) + geom_abline(slope=1,linetype="dashed")
p <- grid.arrange(p1,p2,p3,ncol=3)
# save
ggsave(p,filename=paste0(folder_name,"sigmas_distance.png"),width=15,height=3.5)
# remove
rm(tmp_sigl,tmp_sigu,tmp_sigmas,lims,point_size,p1,p2,p3,p)

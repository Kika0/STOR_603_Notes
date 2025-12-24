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
# plot alpha values for Glasgow, Birmingham, London
cond_site_names <- c("Glasgow","Birmingham","London")
est_sites <- est_all_sf %>% filter(cond_site %in% cond_site_names) %>% mutate(cond_site=factor(cond_site,levels = cond_site_names))
title_map <- ""
misscol <- "aquamarine"
legend_text_size <- 0.7
point_size <- 0.6
legend_title_size <- 1.2
lims <- c(0,1)
nrow_facet <- 1
p <- tm_shape(est_sites) + tm_dots(fill="a",fill.scale = tm_scale_continuous(limits=lims,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\alpha$"))) + tm_facets(by="cond_site",nrow = nrow_facet) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE) + tm_title(text=title_map) 
tmap_save(p,filename=paste0(folder_name,"alpha_selected_sites.png"))
tmap_save(p,filename=paste0(folder_name,"alpha_selected_sites.pdf"))
lims <- c(0,0.5)
p <- tm_shape(est_sites) + tm_dots(fill="b",fill.scale = tm_scale_continuous(limits=lims,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\beta$"))) + tm_facets(by="cond_site",nrow = nrow_facet) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE) + tm_title(text=title_map) 
tmap_save(p,filename=paste0(folder_name,"beta_selected_sites.png"))
tmap_save(p,filename=paste0(folder_name,"beta_selected_sites.pdf"))
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
site_name_diagonal <- c("Birmingham", paste0("diag",1:(length(sites_index_diagonal)-2)),"London")
#spatial_par_est(data_Lap = data_mod_Lap,cond_sites = sites_index_diagonal,cond_site_names = site_name_diagonal,dayshift = c(-3:3),v=q,Ndays_season = 90,title = paste0("diagonal_sites_Birmingham_London",q*100))
# load estimates
load("data_processed/N9000_sequential2_AGG_diagonal_sites_Birmingham_London90.RData")
title_map <- ""
misscol <- "aquamarine"
legend_text_size <- 0.4
point_size <- 0.5
legend_title_size <- 0.6
lims <- c(0,1)
nrow_facet <- 1
est_sites <- est_all_sf %>% filter(tau %in% c(0)) %>% mutate(cond_site=factor(cond_site,levels = site_name_diagonal))
p <- tm_shape(est_sites) + tm_dots(fill="a",fill.scale = tm_scale_continuous(limits=lims,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\alpha$"))) + tm_facets(by="cond_site",nrow = nrow_facet) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE) + tm_title(text=title_map) 
tmap_save(p,filename=paste0(folder_name,"alpha_diagonal_sites.png"))

est_sites <- est_all_sf %>% filter(cond_site %in% "Birmingham") %>% mutate(cond_site=factor(cond_site,levels = site_name_diagonal),tau=factor(tau,levels=c(-3:3)))
p <- tm_shape(est_sites) + tm_dots(fill="a",fill.scale = tm_scale_continuous(limits=lims,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\alpha$"))) + tm_facets(by="tau",nrow = nrow_facet) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE) + tm_title(text=title_map) 
tmap_save(p,filename=paste0(folder_name,"alpha_Brimingham_temporal.png"))

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
load("data_processed/N9000_sequential2_AGG_diagonal_sites_Birmingham_Cromer90.RData",verbose = TRUE)

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
t1a <- tmpa %>% filter(parameter=="a_ite") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\alpha$"))) + tm_facets("parameter") +
  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[1],legend.reverse = TRUE) 
t2a <- tmpa %>% filter(parameter=="a") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\alpha$"))) + tm_facets("parameter") +
  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[2],legend.reverse = TRUE) 
limits <- c(-0.3,0.3)
t3a <- tmpa %>% filter(parameter=="adiff") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="-brewer.rd_bu",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\alpha$"))) + tm_facets("parameter") +
  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[3],legend.reverse = TRUE) 
t <- tmap_arrange(t1a,t2a,t3a,ncol=3)
tmap_save(t,filename=paste0(folder_name,"/a_",cond_site_name,"_iterative_original_difference.png"),width=8,height=6)
tmap_save(t,filename=paste0(folder_name,"/a_",cond_site_name,"_iterative_original_difference.pdf"),width=8,height=6)

tmpb <- tmpsf %>% dplyr::select(b_ite,b,bdiff) %>% pivot_longer(cols=c(b_ite,b,bdiff),names_to = "parameter", values_to = "value" ) %>% mutate(parameter=factor(parameter,levels=c("b_ite","b","bdiff"))) 
limits <- c(0,0.5)
t1b <- tmpb %>% filter(parameter=="b_ite") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\beta$"))) + tm_facets("parameter") +
  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[1],legend.reverse = TRUE) 
t2b <- tmpb %>% filter(parameter=="b") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\beta$"))) + tm_facets("parameter") +
  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[2],legend.reverse = TRUE) 
limits <- c(-0.5,0.5)
t3b <- tmpb %>% filter(parameter=="bdiff") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="-brewer.rd_bu",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\beta$"))) + tm_facets("parameter") +
  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[3],legend.reverse = TRUE) 
t <- tmap_arrange(t1b,t2b,t3b,ncol=3)
tmap_save(t,filename=paste0(folder_name,"/b_",cond_site_name,"_iterative_original_difference.png"),width=8,height=6)
tmap_save(t,filename=paste0(folder_name,"/b_",cond_site_name,"_iterative_original_difference.pdf"),width=8,height=6)

tmpmu <- tmpsf %>% dplyr::select(mu_agg_ite,mu_agg,mudiff) %>% pivot_longer(cols=c(mu_agg_ite,mu_agg,mudiff),names_to = "parameter", values_to = "value" ) %>% mutate(parameter=factor(parameter,levels=c("mu_agg_ite","mu_agg","mudiff"))) 
limits <- c(-2,2)
t1mu <- tmpmu %>% filter(parameter=="mu_agg_ite") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\mu_{AGG}$"))) + tm_facets("parameter") +
  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[1],legend.reverse = TRUE) 
t2mu <- tmpmu %>% filter(parameter=="mu_agg") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\mu_{AGG}$"))) + tm_facets("parameter") +
  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[2],legend.reverse = TRUE) 
limits <- c(-2,2)
t3mu <- tmpmu %>% filter(parameter=="mudiff") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="-brewer.rd_bu",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\mu_{AGG}$"))) + tm_facets("parameter") +
  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[3],legend.reverse = TRUE) 
tmu <- tmap_arrange(t1mu,t2mu,t3mu,ncol=3)
tmap_save(tmu,filename=paste0(folder_name,"/mu_agg_",cond_site_name,"_iterative_original_difference.png"),width=8,height=6)
tmap_save(tmu,filename=paste0(folder_name,"/mu_agg_",cond_site_name,"_iterative_original_difference.pdf"),width=8,height=6)

tmpsigl <- tmpsf %>% dplyr::select(sigl_ite,sigl,sigldiff) %>% pivot_longer(cols=c(sigl_ite,sigl,sigldiff),names_to = "parameter", values_to = "value" ) %>% mutate(parameter=factor(parameter,levels=c("sigl_ite","sigl","sigldiff"))) 
limits <- c(0,4)
t1sigl <- tmpsigl %>% filter(parameter=="sigl_ite") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\sigma_l$"))) + tm_facets("parameter") +
  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[1],legend.reverse = TRUE) 
t2sigl <- tmpsigl %>% filter(parameter=="sigl") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\sigma_l$"))) + tm_facets("parameter") +
  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[2],legend.reverse = TRUE) 
limits <- c(-3,3)
t3sigl <- tmpsigl %>% filter(parameter=="sigldiff") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="-brewer.rd_bu",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\sigma_l$"))) + tm_facets("parameter") +
  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[3],legend.reverse = TRUE) 
tsigl <- tmap_arrange(t1sigl,t2sigl,t3sigl,ncol=3)
tmap_save(tsigl,filename=paste0(folder_name,"/sigl_",cond_site_name,"_iterative_original_difference.png"),width=8,height=6)
tmap_save(tsigl,filename=paste0(folder_name,"/sigl_",cond_site_name,"_iterative_original_difference.pdf"),width=8,height=6)

tmpsigu <- tmpsf %>% dplyr::select(sigu_ite,sigu,sigudiff) %>% pivot_longer(cols=c(sigu_ite,sigu,sigudiff),names_to = "parameter", values_to = "value" ) %>% mutate(parameter=factor(parameter,levels=c("sigu_ite","sigu","sigudiff"))) 
limits <- c(0,4)
t1sigu <- tmpsigu %>% filter(parameter=="sigu_ite") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\sigma_u$"))) + tm_facets("parameter") +
  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[1],legend.reverse = TRUE) 
t2sigu <- tmpsigu %>% filter(parameter=="sigu") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\sigma_u$"))) + tm_facets("parameter") +
  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[2],legend.reverse = TRUE) 
limits <- c(-3,3)
t3sigu <- tmpsigu %>% filter(parameter=="sigudiff") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="-brewer.rd_bu",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\sigma_u$"))) + tm_facets("parameter") +
  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[3],legend.reverse = TRUE) 
tsigu <- tmap_arrange(t1sigu,t2sigu,t3sigu,ncol=3)
tmap_save(tsigu,filename=paste0(folder_name,"/sigu_",cond_site_name,"_iterative_original_difference.png"),width=8,height=6)
tmap_save(tsigu,filename=paste0(folder_name,"/sigu_",cond_site_name,"_iterative_original_difference.pdf"),width=8,height=6)

legend_text_size <- 0.3
point_size <- 0.2
legend_title_size <- 0.4
    tmpa <- tmpsf %>% dplyr::select(a_ite,a,adiff) %>% pivot_longer(cols=c(a_ite,a,adiff),names_to = "parameter", values_to = "value" ) %>% mutate(parameter=factor(parameter,levels=c("a_ite","a","adiff"))) 
  limits <- c(0,1)
  t1a <- tmpa %>% filter(parameter=="a_ite") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\alpha$"))) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[1],legend.reverse = TRUE) 
  t2a <- tmpa %>% filter(parameter=="a") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\alpha$"))) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[2],legend.reverse = TRUE) 
  limits <- c(-0.3,0.3)
  t3a <- tmpa %>% filter(parameter=="adiff") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="-brewer.rd_bu",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\alpha$"))) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[3],legend.reverse = TRUE) 

  tmpb <- tmpsf %>% dplyr::select(b_ite,b,bdiff) %>% pivot_longer(cols=c(b_ite,b,bdiff),names_to = "parameter", values_to = "value" ) %>% mutate(parameter=factor(parameter,levels=c("b_ite","b","bdiff"))) 
  limits <- c(0,0.5)
  t1b <- tmpb %>% filter(parameter=="b_ite") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\beta$"))) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[1],legend.reverse = TRUE) 
  t2b <- tmpb %>% filter(parameter=="b") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\beta$"))) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[2],legend.reverse = TRUE) 
  limits <- c(-0.5,0.5)
  t3b <- tmpb %>% filter(parameter=="bdiff") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="-brewer.rd_bu",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\beta$"))) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[3],legend.reverse = TRUE) 

  tmpmu <- tmpsf %>% dplyr::select(mu_agg_ite,mu_agg,mudiff) %>% pivot_longer(cols=c(mu_agg_ite,mu_agg,mudiff),names_to = "parameter", values_to = "value" ) %>% mutate(parameter=factor(parameter,levels=c("mu_agg_ite","mu_agg","mudiff"))) 
  limits <- c(-2,2)
  t1mu <- tmpmu %>% filter(parameter=="mu_agg_ite") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\mu_{AGG}$"))) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[1],legend.reverse = TRUE) 
  t2mu <- tmpmu %>% filter(parameter=="mu_agg") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\mu_{AGG}$"))) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[2],legend.reverse = TRUE) 
  limits <- c(-2,2)
  t3mu <- tmpmu %>% filter(parameter=="mudiff") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="-brewer.rd_bu",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\mu_{AGG}$"))) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[3],legend.reverse = TRUE) 

  tmpsigl <- tmpsf %>% dplyr::select(sigl_ite,sigl,sigldiff) %>% pivot_longer(cols=c(sigl_ite,sigl,sigldiff),names_to = "parameter", values_to = "value" ) %>% mutate(parameter=factor(parameter,levels=c("sigl_ite","sigl","sigldiff"))) 
  limits <- c(0,4)
  t1sigl <- tmpsigl %>% filter(parameter=="sigl_ite") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\sigma_l$"))) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[1],legend.reverse = TRUE) 
  t2sigl <- tmpsigl %>% filter(parameter=="sigl") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\sigma_l$"))) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[2],legend.reverse = TRUE) 
  limits <- c(-3,3)
  t3sigl <- tmpsigl %>% filter(parameter=="sigldiff") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="-brewer.rd_bu",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\sigma_l$"))) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[3],legend.reverse = TRUE) 
  tsigl <- tmap_arrange(t1sigl,t2sigl,t3sigl,ncol=3)

  tmpsigu <- tmpsf %>% dplyr::select(sigu_ite,sigu,sigudiff) %>% pivot_longer(cols=c(sigu_ite,sigu,sigudiff),names_to = "parameter", values_to = "value" ) %>% mutate(parameter=factor(parameter,levels=c("sigu_ite","sigu","sigudiff"))) 
  limits <- c(0,4)
  t1sigu <- tmpsigu %>% filter(parameter=="sigu_ite") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\sigma_u$"))) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[1],legend.reverse = TRUE) 
  t2sigu <- tmpsigu %>% filter(parameter=="sigu") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\sigma_u$"))) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[2],legend.reverse = TRUE) 
  limits <- c(-3,3)
  t3sigu <- tmpsigu %>% filter(parameter=="sigudiff") %>% tm_shape() + tm_dots(fill="value",size=point_size,fill.scale =tm_scale_continuous(limits=limits,values="-brewer.rd_bu",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\sigma_u$"))) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size, panel.labels = toplabel[3],legend.reverse = TRUE) 
  tsigu <- tmap_arrange(t1sigu,t2sigu,t3sigu,ncol=3)

tallpar <- tmap_arrange(t1a,t1b,t1mu,t1sigl,t1sigu,t2a,t2b,t2mu,t2sigl,t2sigu,t3a,t3b,t3mu,t3sigl,t3sigu,ncol=5)
tmap_save(tallpar,filename=paste0(folder_name,"/allpars_",cond_site_name,"_iterative_original_difference.pdf"),width=8,height=8)

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
                                        units::set_units(km)))
  sigd <- data.frame(sig=sig,dist=dist_cond_site) 
  return(sigd %>% mutate(cond_site = site_name_diagonal[i]))
}

lims <- c(0,2.3)
point_size <- 0.7
tmp_sigl <- do.call(rbind,lapply(1:length(site_name_diagonal),FUN=get_sigma_distance,sig_which="sigl")) %>% mutate(cond_site=factor(cond_site,levels=site_name_diagonal))
p1 <- ggplot(tmp_sigl) + geom_point(aes(y=sig,x=dist,col=cond_site),size=point_size) + scale_color_manual(values = c12) + xlab("Distance [km]") + ylab(TeX("$\\sigma_l$")) + theme(legend.position="none") + ylim(lims) + coord_fixed(ratio=300)
tmp_sigu <- do.call(rbind,lapply(1:length(site_name_diagonal),FUN=get_sigma_distance,sig_which="sigu")) %>% mutate(cond_site=factor(cond_site,levels=site_name_diagonal))
p2 <- ggplot(tmp_sigu) + geom_point(aes(y=sig,x=dist,col=cond_site),size=point_size) + scale_color_manual(values = c12) + xlab("Distance [km]") + ylab(TeX("$\\sigma_u$")) + theme(legend.position="none") + ylim(lims) + coord_fixed(ratio=300)
tmp_sigmas <- data.frame("sigl"=tmp_sigl$sig, "sigu" = tmp_sigu$sig, "cond_site" = tmp_sigl$cond_site)
p3 <- ggplot(tmp_sigmas) + geom_point(aes(y=sigu,x=sigl,col=cond_site),size=point_size) + scale_color_manual(values = c12) + xlab(TeX("$\\sigma_l$")) + ylab(TeX("$\\sigma_u$")) + xlim(lims) + ylim(lims) + coord_fixed() + labs(col="Conditioning site") + guides(col=guide_legend(override.aes=list(size=4))) + geom_abline(slope=1,linetype="dashed")
p <- grid.arrange(p1,p2,p3,ncol=3, widths=c(1,1,1.35))
# save
ggsave(p,filename=paste0(folder_name,"sigmas_distance.png"),width=12,height=3.5)
ggsave(p,filename=paste0(folder_name,"sigmas_distance.pdf"),width=12,height=3.5)
# remove
rm(tmp_sigl,tmp_sigu,tmp_sigmas,lims,point_size,p1,p2,p3,p)

# 7. plot of phis against distance
load("data_processed/iterative_abmu_fixed_delta.RData")
tmpdf <- data.frame("value"=numeric(),"cond_site"=character(),"par"=character())
tmp2 <- data.frame("value"=numeric(),"cond_site"=character(),"par"=character())
Nite <- 10
for (subscript in 0:3) {
  
  tmp1 <- sapply(1:length(sites_index_diagonal),FUN=function(i) {
    rbind(tmpdf,data.frame("value"=tmp_fixed_deltas[[i]][[subscript+6]]) %>% mutate("cond_site"=site_name_diagonal[i],"par"=paste0("phi",subscript),"iteration"=1:(Nite+1)))
    
  },simplify=FALSE)
  tmp <- do.call("rbind",tmp1)
  tmp2 <- rbind(tmp2,tmp)
}

tmp2 <- tmp2 %>% mutate(cond_site=factor(cond_site,levels=site_name_diagonal)) %>% filter(iteration==Nite+1) %>% mutate("Parameter"=par)
p1 <- ggplot(tmp2 %>% filter(par %in% c("phi0","phi2")),aes(x=cond_site,y=value)) + geom_line(aes(group=par))+ geom_point(aes(col=par),size=3)+  xlab("") + ylab("") + scale_colour_manual(values = c("#C11432","#009ADA"))
p2 <- ggplot(tmp2 %>% filter(par %in% c("phi1","phi3")),aes(x=cond_site,y=value)) + geom_line(aes(group=par))+ geom_point(aes(col=par),size=3)+  xlab("") + ylab("") + scale_colour_manual(values = c("#FDD10A","#66A64F"))
p <- grid.arrange(p1,p2,ncol=1)
ggsave(p,filename=paste0(folder_name,"allphis_iterativeabmu.png"),width=10,height=5)
ggsave(p,filename=paste0(folder_name,"allphis_iterativeabmu.pdf"),width=10,height=5)

# try do a PP plot
# calculate observed residuals
v <- 0.9
aest <- discard(est_all_sf %>% filter(cond_site==site_name_diagonal[1]) %>% pull(a),is.na)
best <- discard(est_all_sf %>% filter(cond_site==site_name_diagonal[1]) %>% pull(b),is.na)
Zold <- observed_residuals(df=data_mod_Lap,given=sites_index_diagonal[1],v = v,a=aest,b=best)
aest <- tmp_fixed_deltas[[1]][[12]][,c(1)]
best <- tmp_fixed_deltas[[1]][[12]][,c(2)]
Znew <- observed_residuals(df=data_mod_Lap,given=sites_index_diagonal[1],v = v,a=aest,b=best)
U <- c(1:(dim(data_mod_Lap)[1]*(1-v)))/(dim(data_mod_Lap)[1]*(1-v)+1)
# pick a site
sites <- c(1:ncol(data_mod_Lap))[-sites_index_diagonal[1]]
# get estimates new
est_new <- tmp_fixed_deltas[[1]][[12]] %>% add_row(.before=sites_index_diagonal[1])
# set up dataframe
tmp <- data.frame(x=numeric(),y=numeric(),"method"=character(),"res_site"=numeric())
for (i in 1:(ncol(Z)-1)) {
  site <- sites[i]
AGGPars <- as.numeric(unlist(est_new[site,c(3,4,5,10,11)]))
# get estimates old
AGGParsOld <- as.numeric(unlist(st_drop_geometry(est_all_sf)[site,c(13:17)]))
# construct PP plot
xold <- as.numeric(unlist(Zold[,i]))
xnew <- as.numeric(unlist(Znew[,i]))
Um <- pAGG(x=xnew,theta = AGGPars)
Umo <- pAGG(x=xold,theta = AGGParsOld)
tmp_append_old <- data.frame(y=Umo,method="Original_method","res_site"=sites[i]) %>% arrange(y) %>% mutate("x"=U)
tmp_append_new <- data.frame(x=U,y=Um,method="New_iterative_approach","res_site"=sites[i]) %>% arrange(y) %>% mutate("x"=U)
tmp <- rbind(tmp,tmp_append_old,tmp_append_new) 
}

for (i in 1:5) {
  print(tmp_fixed_deltas[[1]][[i]][res_site_over,])
}

p <- ggplot(tmp) + geom_line(aes(x=x,y=y,group=res_site),linewidth=0.01,alpha=0.9)+ geom_abline(slope=1,linetype="dashed") + facet_wrap("method") + xlab("Empirical") + ylab("Model")
ggsave(p,filename=paste0(folder_name,"PP_Birmingham.pdf"),width=10,height=5)

# examine outliers
tmp <- tmp %>% mutate(diag_diff=y-x)
res_site_over <- tmp[tmp$diag_diff==max(tmp$diag_diff),]$res_site
res_site_under <-  tmp[tmp$diag_diff==min(tmp$diag_diff),]$res_site

AGG_underestimate <- append(NA,(tmp %>% group_by(res_site) %>% summarise(n=max(diag_diff,na.rm=TRUE)) %>% pull(n)),after=sites_index_diagonal[1]-1)
tm_shape(xyUK20_sf[1:555,] %>% mutate(diag_diff= ddif)) + tm_dots("AGG_underestimate")
# examine on a map
over_under <- rep(NA,nrow(xyUK20_sf))
over_under[sites_index_diagonal[1]] <- "Conditioning_site"
over_under[res_site_over] <- "Residual_site_underestimate"
over_under[res_site_under] <- "Residual_site_overestimate"
tmp_over_under <- xyUK20_sf %>% mutate("over_under"=over_under)
t <- tm_shape(tmp_over_under) + tm_dots(fill="over_under",fill.scale = tm_scale_categorical(values=c("Residual_site_underestimate"="#C11432","Residual_site_overestimate" = "#009ADA", "Conditioning_site" = "aquamarine"))) 
tmap_save(t,filename=paste0(folder_name,"PP_outliers_map_examine.png"))
# examine over
AGGPars <- as.numeric(unlist(st_drop_geometry(est_all_sf)[res_site_over,c(13:17)]))
Zover <- as.numeric(unlist(Znew[,which(res_site_over==sites)]))
# plot density and kernel smooth
tmp <- data.frame(y=AGG_density(theta = AGGPars,x=seq(-7.5,5,0.01)),x=seq(-7.5,5,0.01))
p1 <- ggplot(data.frame(Zover)) + geom_density(mapping = aes(x=Zover)) + geom_line(data=tmp,mapping = aes(x=x,y=y),col="#C11432") + xlab(TeX(paste0("$Z_{",res_site_over,"}$"))) + ylab("Residual density function") + ggtitle("Underestimation")
# examine under
AGGPars <- as.numeric(unlist(st_drop_geometry(est_all_sf)[res_site_under,c(13:17)]))
Zunder <- as.numeric(unlist(Znew[,which(res_site_under==sites)]))
# plot density and kernel smooth
tmp <- data.frame(y=AGG_density(theta = AGGPars,x=seq(-7.5,5,0.01)),x=seq(-7.5,5,0.01))
p2 <- ggplot(data.frame(Zunder)) + geom_density(mapping = aes(x=Zunder)) + geom_line(data=tmp,mapping = aes(x=x,y=y),col="#009ADA") + xlab(TeX(paste0("$Z_{",res_site_under,"}$"))) + ylab("Residual density function") + ggtitle("Overestimation")
p <- grid.arrange(p1,p2,ncol=1)
ggsave(p,filename=paste0(folder_name,"AGG_Birmingham_outliers.pdf"),width=10,height=10)

# distance from conditioning side vs worst distance


# repeat for new iterative method


# plot conditional quantiles
Lap_max <- max(apply(X=data_mod_Lap,MARGIN = c(2),FUN = max))
conditional_quantiles <- function(x,model_pars,q_res=0.75,Y=data_mod_Lap, cond_site) {
  a <- model_pars$a
  b <- model_pars$b
  mu <- model_pars$mu
  sigl <- model_pars$sigl
  sigu <- model_pars$sigu
  deltal <- model_pars$deltal 
  deltau <- model_pars$deltau
  Zobs <- observed_residuals(df=Y,given = cond_site,a = a[!is.na(a)],b=b[!is.na(b)])
  Zq1 <- apply(X=Zobs,MARGIN = 2,FUN=function(z)quantile(z,p=q_res))
  Zq1 <- append(Zq1,NA,after = cond_site-1)
  tm_shape(xyUK20_sf %>% mutate(Zq1=Zq1)) + tm_dots("Zq1")
  Zq <- sapply(1:length(a),FUN=function(i)qAGG(p=q_res,theta = c(mu[i],sigl[i],sigu[i],deltal[i],deltau[i])))
  #Zq <- rnorm(length(a))
  return(a*x+ x^b*Zq)  
}
cond_quantiles_wrapper <- function(Lap_max,par_est_site,q_res,site_index,cond_site_name = "Birmingham",grid=xyUK20_sf) {
  cond_site <- site_index  
  x <- Lap_max
  model_pars <- list("a"=par_est_site$a, "b" = par_est_site$b, "mu" = par_est_site$mu_agg, "sigl" = par_est_site$sigl, "sigu" = par_est_site$sigu, "deltal" = par_est_site$deltal, "deltau" = par_est_site$deltau)
  return(data.frame(Ycond=conditional_quantiles(x=x,model_pars=model_pars,q_res=q_res,cond_site=cond_site),period=as.character(T)))
}

for (i in 1:length(sites_index_diagonal)) {
par_est_new <- tmp_fixed_deltas[[i]][[12]] %>% add_row(.before=sites_index_diagonal[i])
Q25 <- cond_quantiles_wrapper(Lap_max=Lap_max,par_est_site = par_est_new,q_res=0.25,site_index = sites_index_diagonal[i],cond_site_name = site_name_diagonal[i])
Q75 <- cond_quantiles_wrapper(Lap_max=Lap_max,par_est_site = par_est_new,q_res=0.75,site_index = sites_index_diagonal[i],cond_site_name = site_name_diagonal[i])
tmp <- rbind(xyUK20_sf %>% mutate(Q=Q25$Ycond,"Quantile"="q=0.25"),xyUK20_sf %>% mutate(Q=Q75$Ycond,"Quantile"="q=0.75"))
lims <- c(0,15)
tnew <- tm_shape(tmp)  + tm_dots(fill="Q",fill.scale = tm_scale_continuous(midpoint=Lap_max,limits=lims,values="-brewer.rd_bu",value.na=misscol,label.na = "Conditioning site"),size=point_size, fill.legend = tm_legend(title="Temperature \n (Laplace scale)")) + tm_layout(legend.outside.size=0.3,asp=0.5,legend.text.size = 1,legend.title.size=1.5,legend.reverse = TRUE,legend.position = tm_pos_out("right","center",pos.h="left",pos.v="top")) + tm_title(site_name_diagonal[i]) + tm_facets(by="Quantile")
tmap_save(tnew,filename=paste0(folder_name,"Qnew_q",q*100,"_",site_name_diagonal[i],".png"),width=9,height=7)
tmap_save(tnew,filename=paste0(folder_name,"Qnew_q",q*100,"_",site_name_diagonal[i],".pdf"),width=9,height=7)

  par_est_new <- est_all_sf %>% filter(cond_site==site_name_diagonal[i])
  Q25 <- cond_quantiles_wrapper(Lap_max=Lap_max,par_est_site = par_est_new,q_res=0.25,site_index = sites_index_diagonal[i],cond_site_name = site_name_diagonal[i])
  Q75 <- cond_quantiles_wrapper(Lap_max=Lap_max,par_est_site = par_est_new,q_res=0.75,site_index = sites_index_diagonal[i],cond_site_name = site_name_diagonal[i])
  tmp <- rbind(xyUK20_sf %>% mutate(Q=Q25$Ycond,"Quantile"="q=0.25"),xyUK20_sf %>% mutate(Q=Q75$Ycond,"Quantile"="q=0.75"))
  told <- tm_shape(tmp)  + tm_dots(fill="Q",fill.scale = tm_scale_continuous(midpoint=Lap_max,limits=lims,values="-brewer.rd_bu",value.na=misscol,label.na = "Conditioning site"),size=point_size, fill.legend = tm_legend(title="Temperature \n (Laplace scale)")) + tm_layout(legend.outside.size=0.3,asp=0.5,legend.text.size = 1,legend.title.size=1.5,legend.reverse = TRUE,legend.position = tm_pos_out("right","center",pos.h="left",pos.v="top")) + tm_title(site_name_diagonal[i]) + tm_facets(by="Quantile")
  tmap_save(told,filename=paste0(folder_name,"Qold_q",q*100,"_",site_name_diagonal[i],".png"),width=9,height=7)
  tmap_save(told,filename=paste0(folder_name,"Qold_q",q*100,"_",site_name_diagonal[i],".pdf"),width=9,height=7)
  tnew <- tnew + tm_title("New iterative approach")
  told <- told + tm_title("Original method")
  tboth <- tmap_arrange(told,tnew,ncol=2)
  tmap_save(tboth,filename=paste0(folder_name,"Q_q",q*100,"_",site_name_diagonal[i],".png"),width=15,height=7)
  tmap_save(tboth,filename=paste0(folder_name,"Q_q",q*100,"_",site_name_diagonal[i],".pdf"),width=15,height=7)
  
}

# transform to original margins?
load("../luna/kristina/P2q/ukgd_cpm85_5k_x84y20_MSp2q.RData",verb=TRUE)

load("../luna/kristina/MSGpdParam/ukgd_cpm85_5k_x84y20.MSGpdParam.2025-02-26.RData",verb=TRUE)

str(gpdpar)
qgam.p2q.fn[[22846]](0.99)

gpdpar[22846,]



files <- list.files("../luna/kristina/MSdata01/")
list_of_files <- list() #create empty list
# only subset for x in x20 and y in y20
files_subset <- sapply(1:nrow(xyUK20_sf),function(i){paste0("ukgd_cpm85_5k_x",xyUK20_sf$x[i],"y",xyUK20_sf$y[i],"_MSdata01.RData")})
#loop through the files
files_subset1 <- files_subset[files_subset %in% files]
# could take only x and y divisible by 4 to subset and speed up data loading
for (i in 1:length(files_subset1)) {
  print(files_subset1[i])
  load(paste0("../luna/kristina/MSdata01/", files_subset1[i]))
  list_of_files[[i]] <- data01 #add files to list position
}
xyUK20_sf <- xyUK20_sf[files_subset %in% files_subset1,]

# join the summer data together could try 1960-1999 for non-stationary data -----
# June 1 is 152 doy, August 31 is 243 doy (92 days per year)
# create a dataframe
xyUK20 <- xyUK20_sf %>% dplyr::select(-temp) %>% st_drop_geometry() %>% cbind(as.data.frame(matrix(data=numeric(),ncol=92*40,nrow=nrow(xyUK20_sf))))
names(xyUK20)[5:ncol(xyUK20)] <- paste0(rep(152:243,40),"_",rep(1960:1999,each=92)) 
for (i in 1: length(list_of_files)) {
  xyUK20[i,5:ncol(xyUK20)] <- list_of_files[[i]] %>% mutate(year=floor(time)) %>% filter(class=="obs",doy>=152,doy<=243, year<=1999) %>% pull(x) 
}

# explore issue with deltas Birmingham to Cromer
summary(est_all_sf)
est_all_sf %>% 
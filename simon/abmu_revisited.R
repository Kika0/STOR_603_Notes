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

# load iterative sigmal estimates
load("data_processed/iterative_sigmal_estimates.RData")

# spatial_par_est_abmu(data_Lap = data_mod_Lap,cond_sites = df_sites,res_margin_est = iterative_sigmal_estimates,dayshift = 0,v=q,Ndays_season = 90,title = paste0("abmu_revisited",q*100))

# examine the output
load("data_processed/N9000_abmu_revisited90.RData")
summary(est_all_sf)
est_all_sf <- est_all_sf %>% mutate(cond_site=factor(cond_site,levels=names(df_sites)))
# map parameter estimates for each site
map_param_abmu <- function(tmp_est,facet_var = "cond_site",filename_part="abmu_revisited",title_map="") {
    misscol <- "aquamarine"
    Nsites <- max(tmp_est$res, tmp_est$given,na.rm=TRUE)
    if (identical(facet_var,"cond_site")) {
      Nfacet <- length(unique(tmp_est$cond_site))
      facet_label <- unique(tmp_est$cond_site)
      nrow_facet <- ceiling(length(unique(tmp_est$cond_site))/4)
      legend_outside_size <- 0.3
    } 
    rep_Nsites <- nrow(tmp_est)/Nsites
    uk_tmp1 <- tmp_est 

#      lims <- c(min(as.numeric(uk_tmp1$a),na.rm = TRUE),max(as.numeric(uk_tmp1$a),na.rm = TRUE))
      lims <- c(-1,1)
      pa <- tm_shape(uk_tmp1) + tm_dots(fill="a",fill.scale = tm_scale_continuous(limits=lims,values="-matplotlib.rd_bu",value.na=misscol,label.na = "Conditioning site"),size=0.3, fill.legend = tm_legend(title=TeX("$\\alpha$"))) + tm_facets(by=facet_var,nrow = nrow_facet) +  tm_layout(panel.labels = facet_label,legend.outside.size=legend_outside_size,asp=0.5,legend.text.size = 0.8,legend.title.size=1,legend.reverse=TRUE) + tm_title(text=title_map) 
 #     lims <- c(min(as.numeric(uk_tmp1$b),na.rm = TRUE),max(as.numeric(uk_tmp1$b),na.rm = TRUE))
      lims <- c(0,1)
      pb <- tm_shape(uk_tmp1) + tm_dots(fill="b",fill.scale = tm_scale_continuous(limits=lims,values="Blues",value.na=misscol,label.na = "Conditioning site"),size=0.3, fill.legend = tm_legend(title=TeX("$\\beta$"))) + tm_facets(by=facet_var,nrow = nrow_facet) +  tm_layout(panel.labels = facet_label,legend.outside.size=legend_outside_size,asp=0.5,legend.text.size = 0.8,legend.title.size=1,legend.reverse=TRUE) + tm_title(text=title_map) 
 #    lims <- c(min(as.numeric(uk_tmp1$mu),na.rm = TRUE),max(as.numeric(uk_tmp1$mu),na.rm = TRUE))
      pmu <- tm_shape(uk_tmp1) + tm_dots(fill="mu",fill.scale = tm_scale_continuous(values="-matplotlib.rd_bu",value.na=misscol,label.na = "Conditioning site"),size=0.3, fill.legend = tm_legend(title=TeX("$\\mu$"))) + tm_facets(by=facet_var,nrow = nrow_facet) +  tm_layout(panel.labels = facet_label,legend.outside.size=legend_outside_size,asp=0.5,legend.text.size = 0.8,legend.title.size=1,legend.reverse=TRUE) + tm_title(text=title_map) 
      tmap_save(pa,filename=paste0("../Documents/a_all_",filename_part,".png"),width=8,height=8)
      tmap_save(pb,filename=paste0("../Documents/b_all_",filename_part,".png"),width=8,height=8)
      tmap_save(pmu,filename=paste0("../Documents/mu_all_",filename_part,".png"),width=8,height=8)
      
      return(list(pa,pb,pmu)) 
}

maptit <- TeX("$(\\alpha,\\beta,\\mu)$ estimated at each site for 12 different conditioning sites")
testmap <- map_param_abmu(tmp_est = est_all_sf, title_map = maptit)
testmap[[3]]

# examine outlier
est_all_sfout <- est_all_sf[which.max(est_all_sf$a),]
est_all_sfout

# compare with previous estimates
load("data_processed/N9000_sequential2_AGG_all12sites90.RData")
summary(est_all_sf)
est_all_sf <- est_all_sf %>% mutate(cond_site=factor(cond_site,levels=names(df_sites)))
maptit <- TeX("$(\\alpha,\\beta,\\mu)$ for 12 different conditioning sites (original estimates)")
testmap <- map_param_abmu(tmp_est = est_all_sf,filename_part = "original_estimates",title_map = maptit)

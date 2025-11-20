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

sites_index_diagonal <- c(192,193,194,195,196,197,218,219,220,221,242,263) # first is Birmingham and last is Cromer
site_name_diagonal <- c("Birmingham", paste0("diagonal",1:(length(sites_index_diagonal)-2)),"Cromer")

# load other estimates
load(file="data_processed/iterative_sigmas_estimates_Birmingham_Cromer_diagonal.RData")
deltal <- sapply(1:length(sites_index_diagonal),FUN = function (i) as.numeric(st_drop_geometry( result[[i]][1,29])))[1]
deltau <- sapply(1:length(sites_index_diagonal),FUN = function (i) as.numeric(st_drop_geometry( result[[i]][1,30])))[1]

pe <- as.data.frame(result[[1]] %>% select(a,b,mu_agg_ite_sig,sigl_ite_sig,sigu_ite_sig,deltal_ite,deltau_ite)) %>% drop_na()
names(pe) <- c("a","b","mu","sigl","sigu","deltal","deltau")
# calculate observed residuals
v <- 0.9
Z <- observed_residuals(df=data_mod_Lap,given=sites_index_diagonal[1],v = v,a=pe$a,b=pe$b)

# transform to standard Normal
AGG_Normal_PIT <- function(z,theta) {
  return(qnorm(pAGG(x=z,theta=theta)))
}

ZN <- sapply(1:ncol(Z),FUN=function(k) {AGG_Normal_PIT(z = Z[,k],theta=c(pe$mu[k],pe$sigl[k],pe$sigu[k],pe$deltal[1],pe$deltau[1]))})
# add empty row for Birmingham and transform back to spatial
est_ite <- as.data.frame(t(ZN)) %>% add_row(.before=sites_index_diagonal[1])
# find 10 largest days
hot_temps <- sort(as.numeric(unlist(data_mod_Lap[,sites_index_diagonal[1]][data_mod_Lap[,sites_index_diagonal[1]]>quantile(data_mod_Lap[,sites_index_diagonal[1]],v)])),decreasing=TRUE)[1:10]
# find the corresponding rows
hot_temps_index <- sort(as.numeric(unlist(data_mod_Lap[,sites_index_diagonal[1]][data_mod_Lap[,sites_index_diagonal[1]]>quantile(data_mod_Lap[,sites_index_diagonal[1]],v)])),decreasing=TRUE,index.return=TRUE)$ix[1:10]

tmp <- st_as_sf(cbind(est_ite %>% select(all_of(hot_temps_index)),result[[1]] %>% dplyr::select(geometry)) %>% pivot_longer(cols=starts_with("V")))

t <- tm_shape(tmp %>% mutate("name"=factor(name, levels=unique(tmp$name)))) + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="-brewer.rd_bu",limits=c(-8,8)),fill.legend = tm_legend(title = TeX("$Z^N$"))) + tm_facets(by="name",ncol=5)
tmap_save(tm=t, filename=paste0("../Documents/residual_dependence_Birmingham.png"),width=10,height=8)

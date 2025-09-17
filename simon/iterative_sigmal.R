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

# repeat for all other sites -------------------------------------------------
iter_sigmal_site <- function(i,Nite=5,file_folder="iterative_sigmal_res_margin",sites=df_sites,grid=xyUK20_sf,data=data_mod_Lap,par_est=est_all_sf,ite_sigu = iterative_sigmau_estimates,index_outliers = NULL) {
  q <- 0.9
  cond_site_name <- names(sites)[i] 
  cond_site_coord <- sites %>% dplyr::select(all_of(cond_site_name)) %>% pull()
  
  if (is.null(index_outliers)) {
    est_site <- par_est %>% filter(cond_site==cond_site_name)
    tmpsf <- ite_sigu[[ which(names(df_sites)==cond_site_name) ]]
  } else {
    # change NA to FALSE for subsetting the points
    index_outliers[is.na(index_outliers)] <- FALSE
    est_site <- par_est %>% filter(cond_site==cond_site_name)
    est_site <- est_site[!index_outliers,] 
    tmpsf <- ite_sigu[[ which(names(df_sites)==cond_site_name) ]][!index_outliers,]
    # subset data
    data_mod_Lap <- data_mod_Lap[,!index_outliers]
    # subset grid
    grid <- grid[!index_outliers,]
  }
  cond_site <- find_site_index(cond_site_coord,grid_uk = grid)
  Z <- observed_residuals(df = data_mod_Lap, given = cond_site, v = q,a= discard(as.numeric(est_site$a),is.na),b = discard(as.numeric(est_site$b),is.na))
  # calculate distance from the conditioning site
  dist_tmp <- as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))
  # remove zero distance
  dist_tmp <- dist_tmp[dist_tmp>0]
  # normalise distance using a common constant
  distnorm <- dist_tmp/1000000
  tmp <- sigmal_par_est_ite(data = Z, given = cond_site, cond_site_dist = distnorm, Nite = Nite, show_ite = TRUE, mu_init = discard(as.numeric(tmpsf$mu_agg_ite_sigu),is.na), sigl_init = discard(as.numeric(tmpsf$sigl_ite_sigu),is.na), deltal = as.numeric(tmpsf$deltal_ite[1]), deltau = as.numeric(tmpsf$deltau_ite[1]), sigu = discard(as.numeric(tmpsf$sigu_ite_sigu),is.na))
  
  # explore estimates
  # plot phi estimates
  tmp_phi <- rbind(data.frame("phi"=as.numeric(unlist(tmp[[3]])),iteration=1:(Nite+1),parameter = "phi2"),
                   data.frame("phi"=as.numeric(unlist(tmp[[4]])),iteration=1:(Nite+1),parameter = "phi3"))
  pphi <- ggplot(tmp_phi) + geom_point(aes(x=iteration,y=phi,col=parameter)) + scale_color_manual(values = c("#009ADA","#66A64F"), breaks = c("phi2","phi2"),labels = c(TeX("$\\phi_2$"),TeX("$\\phi_3$"))) + ylab("")
  ggsave(pphi,filename=paste0("../Documents/",file_folder,"/phi_",cond_site_name,".png"),width=5,height=5)
  
  # plot also for a random site as a check for convergence
  # plot across iterations for a selected site
  random_site <- 100
  tmp_all <- rbind(data.frame(delta=as.numeric(unlist(tmp[[1]][random_site,])),iteration=1:(Nite+1),parameter = "mu_agg"),
                   data.frame(delta=as.numeric(unlist(tmp[[2]][random_site,])),iteration=1:(Nite+1),parameter = "sigma_lower"),
                   data.frame(delta=as.numeric(unlist(tmp[[3]])),iteration=1:(Nite+1),parameter = "phi2"),
                   data.frame(delta=as.numeric(unlist(tmp[[4]])),iteration=1:(Nite+1),parameter = "phi3")
                   
  )
  # plot the final iteration
  p1 <- ggplot(tmp_all) + geom_point(aes(x=iteration,y=delta,col=parameter)) + ylab("")
  ggsave(p1,filename=paste0("../Documents/",file_folder,"/random_site_par_",cond_site_name,".png"),width=5,height=5)
  
  # explore also spatial parameters
  est_ite <- tmp[[5]] %>% add_row(.before=cond_site)
  names(est_ite) <- paste0(names(est_ite),"_ite_sigl")
  tmpsf <- cbind(tmpsf,est_ite)
  # plot parameter estimates against distance
  # calculate distance from a conditioning site with st_distance()
  mud <- data.frame(mu=tmpsf$mu_agg_ite_sigl,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=mu,x=dist))
  ggsave(mud,filename=paste0("../Documents/",file_folder,"/muagg_distance_",cond_site_name,".png")) 
  
  sigld <- data.frame(sigl=tmpsf$sigl_ite_sigl,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=sigl,x=dist))
  ggsave(sigld,filename=paste0("../Documents/",file_folder,"/sigl_distance_",cond_site_name,".png")) 
  
 
  tmpsf <- tmpsf %>% mutate(mudiff = mu_agg_ite_sigu - mu_agg_ite_sigl, sigldiff = sigl_ite_sigu - sigl_ite_sigl)
  toplabel <- c(TeX("Iterative $\\sigma_u (d_j)$ method"),TeX("Iterative $\\sigma_l (d_j)$ method"),"Difference")
  mu_limits <- c(-1.75,1.79)
  t <- tmpsf %>% dplyr::select(mu_agg_ite_sigu,mu_agg_ite_sigl,mudiff) %>% pivot_longer(cols=c(mu_agg_ite_sigu,mu_agg_ite_sigl,mudiff),names_to = "parameter", values_to = "value" ) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="-brewer.rd_bu",limits=mu_limits),fill.legend = tm_legend(title = TeX("$\\mu_{AGG}$"))) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 12, panel.labels = toplabel) 
  
  tmap_save(t,filename=paste0("../Documents/",file_folder,"/mu_agg_",cond_site_name,".png"),width=8,height=6)
  
  sigma_limits <- c(-1.32,2.3)
  t <- tmpsf %>% dplyr::select(sigl_ite_sigu,sigl_ite_sigl,sigldiff) %>% pivot_longer(cols=c(sigl_ite_sigu,sigl_ite_sigl,sigldiff),names_to = "parameter", values_to = "value" ) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="-brewer.rd_bu",limits=sigma_limits),fill.legend = tm_legend(title = TeX("$\\sigma_l$"))) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 9, panel.labels = toplabel) 
  tmap_save(t,filename=paste0("../Documents/",file_folder,"/sigma_lower_",cond_site_name,".png"),width=8,height=6)
  
  return(tmpsf)
}

Nite <- 5
tmp <- sapply(1:ncol(df_sites),iter_sigmal_site, Nite = 5, simplify = FALSE)

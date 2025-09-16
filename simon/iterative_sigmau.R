sigmau_par_est_ite <- function(data=Z,given=cond_site,cond_site_dist, Nite=10, show_ite=FALSE,mu_init=NULL,sigl_init=NULL,sigu_init=NULL,deltal,deltau)  {
  d <- ncol(data)
  N <- nrow(data)
  res <- 1:d
  mu_agg <- sigl <- sigu <- data.frame(matrix(ncol=(Nite+1),nrow = d))
  # calculate a with initial values for mu and sigma
  if (is.numeric(mu_init) & is.numeric(sigl_init) & is.numeric(sigl_init)) {
    mu_agg[,1] <- mu_init
    sigl[,1] <- sigl_init
    sigu[,1] <- sigu_init
  } else {
    mu_agg[,1] <- 0
    sigl[,1] <- 1
    sigu[,1] <- 1
  }
phi0. <- phi1. <- c(1) 
  for (i in 1:Nite) {
    # estimate sigu parameters phi0 and phi1
    phi_init <- c(phi0.[i],phi1.[i])
    opt <- optim(fn=NLL_exp_sigmau,x = data,d1j=cond_site_dist,mu1=as.numeric(mu_agg[,i]),sigl1=as.numeric(sigl[,i]),deltal=deltal,deltau=deltau,control=list(maxit=2000),par = phi_init,method = "Nelder-Mead")
    phi0 <- opt$par[1]
    phi1 <- opt$par[2]
    phi0. <- append(phi0.,phi0)
    phi1. <- append(phi1.,phi1)
    sigu[,i+1] <- phi0*(1-exp(-phi1*cond_site_dist))
    # estimate sigl and mu_agg for each site separately
    for (j in 1:d) {
      Z2 <- as.numeric(unlist(data[,j]))
      opt <- optim(fn=NLL_AGG,x=Z2,sigu_hat = as.numeric(unlist(sigu[,i+1])),deltal_hat=deltal,deltau_hat = deltau,par=c(mean(Z2),sd(Z2),sd(Z2)),control=list(maxit=2000),method = "Nelder-Mead")
      mu_agg[j,i+1] <- opt$par[1]
      sigl[j,i+1] <- opt$par[2]
    }

  }
  par_sum <- data.frame("mu_agg" = as.numeric(mu_agg[,Nite+1]),"sigl" = as.numeric(sigl[,Nite+1]),"sigu" = as.numeric(sigu[,Nite+1]),"phi0" = phi0.[Nite+1], "phi1" = phi1.[Nite+1])
  if (show_ite == TRUE) {
    return(list(mu_agg,sigl,sigu,phi0.,phi1.,par_sum))
  } else {return(par_sum)}
}


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
load("data_processed/iterative_delta_estimates.RData")

# repeat for all other sites -------------------------------------------------
iter_sigmau_site <- function(i,Nite=5,sites=df_sites,grid=xyUK20_sf,data=data_mod_Lap,par_est=est_all_sf,ite_delta = result,index_outliers = NULL) {
  q <- 0.9
  cond_site_name <- names(sites)[i] 
  cond_site_coord <- sites %>% dplyr::select(all_of(cond_site_name)) %>% pull()
  
  if (is.null(index_outliers)) {
  est_site <- par_est %>% filter(cond_site==cond_site_name)
  tmpsf <- ite_delta[[ which(names(df_sites)==cond_site_name) ]]
} else {
    # change NA to FALSE for subsetting the points
    index_outliers[is.na(index_outliers)] <- FALSE
    est_site <- par_est %>% filter(cond_site==cond_site_name)
    est_site <- est_site[!index_outliers,] 
    tmpsf <- ite_delta[[ which(names(df_sites)==cond_site_name) ]][!index_outliers,]
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
  tmp <- sigmau_par_est_ite(data = Z, given = cond_site, cond_site_dist = distnorm, Nite = Nite, show_ite = TRUE, mu_init = discard(as.numeric(tmpsf$mu_agg_ite),is.na), sigl_init = discard(as.numeric(tmpsf$sigl_ite),is.na), sigu_init = discard(as.numeric(tmpsf$sigu_ite),is.na), deltal = as.numeric(tmpsf$deltal_ite[1]), deltau = as.numeric(tmpsf$deltau_ite[1]))
  
  # explore estimates
  # plot phi estimates
  tmp_phi <- rbind(data.frame("phi"=as.numeric(unlist(tmp[[4]])),iteration=1:(Nite+1),parameter = "phi0"),
                   data.frame("phi"=as.numeric(unlist(tmp[[5]])),iteration=1:(Nite+1),parameter = "phi1"))
  pphi <- ggplot(tmp_phi) + geom_point(aes(x=iteration,y=phi,col=parameter)) + scale_color_manual(values = c("#009ADA","#66A64F"), breaks = c("phi0","phi1"),labels = c(TeX("$\\phi_0$"),TeX("$\\phi_1$"))) + ylab("")
  ggsave(pphi,filename=paste0("../Documents/iterative_sigmas_res_margin/phi_",cond_site_name,".png"),width=5,height=5)
  
  # plot also for a random site as a check for convergence
  # plot across iterations for a selected site
  random_site <- 100
  tmp_all <- rbind(data.frame(delta=as.numeric(unlist(tmp[[1]][random_site,])),iteration=1:(Nite+1),parameter = "mu_agg"),
                   data.frame(delta=as.numeric(unlist(tmp[[2]][random_site,])),iteration=1:(Nite+1),parameter = "sigma_lower"),
                   data.frame(delta=as.numeric(unlist(tmp[[3]][random_site,])),iteration=1:(Nite+1),parameter = "sigma_upper"),
                   data.frame(delta=as.numeric(unlist(tmp[[4]])),iteration=1:(Nite+1),parameter = "phi0"),
                   data.frame(delta=as.numeric(unlist(tmp[[5]])),iteration=1:(Nite+1),parameter = "phi1")
                   
  )
  # plot the final iteration
  p1 <- ggplot(tmp_all) + geom_point(aes(x=iteration,y=delta,col=parameter)) + ylab("")
  ggsave(p1,filename=paste0("../Documents/iterative_sigmas_res_margin/random_site_par_",cond_site_name,".png"),width=5,height=5)
  
  # explore also spatial parameters
  est_ite <- tmp[[6]] %>% add_row(.before=cond_site)
  names(est_ite) <- paste0(names(est_ite),"_ite_sigu")
  tmpsf <- cbind(tmpsf,est_ite)
  # plot parameter estimates against distance
  # calculate distance from a conditioning site with st_distance()
  mud <- data.frame(mu=tmpsf$mu_agg_ite_sigu,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=mu,x=dist))
  ggsave(mud,filename=paste0("../Documents/iterative_sigmas_res_margin/muagg_distance_",cond_site_name,".png")) 
  
  sigld <- data.frame(sigl=tmpsf$sigl_ite_sigu,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=sigl,x=dist))
  ggsave(sigld,filename=paste0("../Documents/iterative_sigmas_res_margin/sigl_distance_",cond_site_name,".png")) 
  
  sigud <- data.frame(sigu=tmpsf$sigu_ite_sigu,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=sigu,x=dist))
  ggsave(sigud,filename=paste0("../Documents/iterative_sigmas_res_margin/sigu_distance_",cond_site_name,".png")) 
  

  tmpsf <- tmpsf %>% mutate(mudiff = mu_agg_ite - mu_agg_ite_sigu, sigldiff = sigl_ite - sigl_ite_sigu, sigudiff = sigu_ite - sigu_ite_sigu)
  toplabel <- c(TeX("Iterative $\\delta$s method"),TeX("Iterative $\\sigma_u (d_j)$ method"),"Difference")
  mu_limits <- c(-1.75,1.79)
  t <- tmpsf %>% dplyr::select(mu_agg_ite,mu_agg_ite_sigu,mudiff) %>% pivot_longer(cols=c(mu_agg_ite,mu_agg_ite_sigu,mudiff),names_to = "parameter", values_to = "value" ) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="-brewer.rd_bu",limits=mu_limits),fill.legend = tm_legend(title = TeX("$\\mu_{AGG}$"))) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 12, panel.labels = toplabel) 
  
  tmap_save(t,filename=paste0("../Documents/iterative_sigmas_res_margin/mu_agg_",cond_site_name,".png"),width=8,height=6)
  
  sigma_limits <- c(-1.32,2.3)
  t <- tmpsf %>% dplyr::select(sigl_ite,sigl_ite_sigu,sigldiff) %>% pivot_longer(cols=c(sigl_ite,sigl_ite_sigu,sigldiff),names_to = "parameter", values_to = "value" ) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="-brewer.rd_bu",limits=sigma_limits),fill.legend = tm_legend(title = TeX("$\\sigma_l$"))) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 9, panel.labels = toplabel) 
  tmap_save(t,filename=paste0("../Documents/iterative_sigmas_res_margin/sigma_lower_",cond_site_name,".png"),width=8,height=6)
  
  t <- tmpsf %>% dplyr::select(sigu_ite,sigu_ite_sigu,sigudiff) %>% pivot_longer(cols=c(sigu_ite,sigu_ite_sigu,sigudiff),names_to = "parameter", values_to = "value" ) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="-brewer.rd_bu",limits=sigma_limits),fill.legend = tm_legend(title = TeX("$\\sigma_u$"))) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 9, panel.labels = toplabel) 
  tmap_save(t,filename=paste0("../Documents/iterative_sigmas_res_margin/sigma_upper_",cond_site_name,".png"),width=8,height=6)
  
  return(tmpsf)
}
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
tm_shape(xyUK20_sf) + tm_dots() + tm_shape(est_phi %>% pivot_longer(cols = c(phi0,phi1),names_to = "parameter", values_to = "value")) + tm_dots(fill="value",size = 2) + tm_facets(by = "parameter") + tm_layout(legend.position=c("right","top"),legend.height = 12,panels.label=toplabel)
tmphi0 <- tm_shape(xyUK20_sf) + tm_dots() + tm_shape(est_phi) + tm_dots(fill="phi0",size = 2, fill.legend = tm_legend(title = TeX("$\\phi_0$"))) + tm_layout(legend.position=c("right","top"),legend.height = 12)
tmphi1 <- tm_shape(xyUK20_sf) + tm_dots() + tm_shape(est_phi) + tm_dots(fill="phi1",size = 2, fill.legend = tm_legend(title = TeX("$\\phi_1$"))) + tm_layout(legend.position=c("right","top"),legend.height = 12)
t <- tmap_arrange(tmphi0,tmphi1,ncol=2)
tmap_save(t,filename=paste0("../Documents/iterative_sigmas_res_margin_all/all_phis.png"),width=8,height=6)

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
ggsave(p,filename=paste0("../Documents/iterative_sigmas_res_margin/sigu_distance_all.png"),width=10,height=7) 

# save these estimates

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



# save these estimates


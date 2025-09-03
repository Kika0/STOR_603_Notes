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
df_sites <- data.frame(Birmingham,Glasgow,London,Inverness,Lancaster,Newcastle,Cromer,Hull,Lowestoft,Truro,Dolgellau,Bournemouth)
#spatial_par_est saves parameter estimates as est_all_sf sf object in ../Documents folder
q <- 0.9 # quantile threshold
# load all three parameter estimates sf objects
load(paste0("data_processed/N9000_sequential2_AGG_all12sites",q*100,".RData"))
est_all <- as.data.frame(est_all_sf)

# load iterative delta estimates
load("data_processed/iterative_delta_estimates.RData")

# select conditioning site
cond_site_name <- "Birmingham"
cond_site_coord <- df_sites %>% dplyr::select(all_of(cond_site_name)) %>% pull()
cond_site <- find_site_index(cond_site_coord,grid_uk = grid)
tmpsf <- result[[ which(names(df_sites)==cond_site_name) ]]
est_site <- par_est %>% filter(cond_site==cond_site_name)
# calculate the observed residuals (set one of 12 sites as conditioning site)
Z <- observed_residuals(df = data_mod_Lap, given = cond_site, v = q,a= discard(as.numeric(est_site$a),is.na),b = discard(as.numeric(est_site$b),is.na))
# calculate distance from the conditioning site
dist_tmp <- as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))
# remove zero distance
dist_tmp <- dist_tmp[dist_tmp>0]
# normalise distance
distnorm <- dist_tmp/max(dist_tmp)
# check also whether units affect phi estimates
opt <- optim(fn=NLL_exp_sigmau,x = data,d1j=distnorm,mu1=as.numeric(mu_agg[,i]),sigl1=as.numeric(sigl[,i]),deltal=deltal[1],deltau=deltau[1],control=list(maxit=2000),par = phi_init,method = "Nelder-Mead")
# look into the log-likelihood function
NLL_exp_sigmau(phi=c(1,1),x = Z,d1j=dist_tmp,mu1=as.numeric(mu_agg[,i]),sigl1=as.numeric(sigl[,i]),deltal=deltal[1],deltau=deltau[1])

# test the function
Nite <- 50
tmp <- sigmau_par_est_ite(data = Z, given = cond_site, cond_site_dist = distnorm, Nite = Nite, show_ite = TRUE, mu_init = discard(as.numeric(tmpsf$mu_agg_ite),is.na), sigl_init = discard(as.numeric(tmpsf$sigl_ite),is.na), sigu_init = discard(as.numeric(tmpsf$sigu_ite),is.na), deltal = as.numeric(tmpsf$deltal_ite[1]), deltau = as.numeric(tmpsf$deltau_ite[1]))

# explore estimates
# plot delta estimates
tmp_phi <- rbind(data.frame("phi"=as.numeric(unlist(tmp[[4]])),iteration=1:(Nite+1),parameter = "phi0"),
                   data.frame("phi"=as.numeric(unlist(tmp[[5]])),iteration=1:(Nite+1),parameter = "delta_upper"))
pphi <- ggplot(tmp_phi) + geom_point(aes(x=iteration,y=phi,col=parameter))
ggsave(pd,filename=paste0("../Documents/iterative_sigmas_res_margin/deltas_",cond_site_name,".png"))

# plot also for a random site as a check for convergence
# plot across iterations for a selected site
random_site <- 100
tmp_all <- rbind(data.frame(delta=as.numeric(unlist(tmp[[1]][random_site,])),iteration=1:(Nite+1),parameter = "mu_agg"),
                 data.frame(delta=as.numeric(unlist(tmp[[2]][random_site,])),iteration=1:(Nite+1),parameter = "sigma_lower"),
                 data.frame(delta=as.numeric(unlist(tmp[[3]][random_site,])),iteration=1:(Nite+1),parameter = "sigma_upper"),
                 data.frame(delta=as.numeric(unlist(tmp[[4]])),iteration=1:(Nite+1),parameter = "delta_lower"),
                 data.frame(delta=as.numeric(unlist(tmp[[5]])),iteration=1:(Nite+1),parameter = "delta_upper")
                 
)
# plot the final iteration
p1 <- ggplot(tmp_all) + geom_point(aes(x=iteration,y=delta,col=parameter))
ggsave(p1,filename=paste0("../Documents/iterative_sigmas_res_margin/random_site_par_",cond_site_name,".png"))

# explore also spatial parameters
est_ite <- tmp[[6]] %>% add_row(.before=cond_site)
names(est_ite) <- paste0(names(est_ite),"_ite_sigu")
tmpsf <- cbind(tmpsf,est_ite)
# plot parameter estimates against distance
# calculate distance from a conditioning site with st_distance()
mud <- data.frame(mu=tmpsf$mu_agg_ite_sigu,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=mu,x=dist))
ggsave(mud,filename=paste0("../Documents/iterative_sigmas_res_margin/muagg_distance_",cond_site_name,".png")) 

sigld <- data.frame(sigl=tmpsf$sigl_ite,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=sigl,x=dist))
ggsave(sigld,filename=paste0("../Documents/iterative_sigmas_res_margin/sigl_distance_",cond_site_name,".png")) 

sigud <- data.frame(sigu=tmpsf$sigu_ite,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=sigu,x=dist))
ggsave(sigud,filename=paste0("../Documents/iterative_sigmas_res_margin/sigu_distance_",cond_site_name,".png")) 

# compare with iterative delta approach
mud1 <- data.frame(mu=tmpsf$mu_agg_ite,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=mu,x=dist))
mud2 <- data.frame(mu=tmpsf$mu_agg_ite_sigu,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=mu,x=dist))
mud3 <- data.frame(mu=tmpsf$mu_agg_ite-tmpsf$mu_agg_ite_sigu,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=mu,x=dist))
sigld <- data.frame(sigl=tmpsf$sigl_ite,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=sigl,x=dist))
sigld <- data.frame(sigl=tmpsf$sigl_ite,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=sigl,x=dist))
sigld <- data.frame(sigl=tmpsf$sigl_ite,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=sigl,x=dist))
sigud <- data.frame(sigu=tmpsf$sigu_ite,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=sigu,x=dist))
sigud <- data.frame(sigu=tmpsf$sigu_ite,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=sigu,x=dist))
sigud <- data.frame(sigu=tmpsf$sigu_ite,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=sigu,x=dist))
library(gridExtra)
grid.arrange(mud1,mud2,mud3,ncol=3)

tmpsf <- tmpsf %>% mutate(mudiff = mu_agg_ite - mu_agg_ite_sigu, sigldiff = sigl_ite - sigl_ite_sigu, sigudiff = sigu_ite - sigu_ite_sigu)
t <- tmpsf %>% dplyr::select(mu_agg_ite,mu_agg_ite_sigu,mudiff) %>% pivot_longer(cols=c(mu_agg_ite,mu_agg_ite_sigu,mudiff),names_to = "parameter", values_to = "value" ) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="-brewer.rd_bu")) + tm_facets("parameter") +
  tm_layout(legend.position=c("right","top"),legend.height = 12) 

tmap_save(t,filename=paste0("../Documents/iterative_deltas_res_margin/mu_agg_",cond_site_name,".png"),width=8,height=6)

t <- tmpsf %>% dplyr::select(sigl_ite,sigl_ite_sigu,sigldiff) %>% pivot_longer(cols=c(sigl_ite,sigl_ite_sigu,sigldiff),names_to = "parameter", values_to = "value" ) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="brewer.blues")) + tm_facets("parameter") +
  tm_layout(legend.position=c("right","top"),legend.height = 12) 
tmap_save(t,filename=paste0("../Documents/iterative_sigmas_res_margin/sigma_lower_",cond_site_name,".png"),width=8,height=6)

t <- tmpsf %>% dplyr::select(sigu_ite,sigu_ite_sigu,sigudiff) %>% pivot_longer(cols=c(sigu_ite,sigu_ite_sigu,sigudiff),names_to = "parameter", values_to = "value" ) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="brewer.blues")) + tm_facets("parameter") +
  tm_layout(legend.position=c("right","top"),legend.height = 12) 
tmap_save(t,filename=paste0("../Documents/iterative_sigmas_res_margin/sigma_upper_",cond_site_name,".png"),width=8,height=6)

# repeat for all other sites
iter_sigmau_site <- function(j,Nite=10,sites=df_sites,grid=xyUK20_sf,data=data_mod_Lap,par_est=est_all_sf) {
  q <- 0.9
  cond_site_name <- names(sites)[j]  
  est_site <- par_est %>% filter(cond_site==cond_site_name)
  cond_site_coord <- sites %>% dplyr::select(all_of(cond_site_name)) %>% pull()
  cond_site <- find_site_index(cond_site_coord,grid_uk = grid)
  Z <- observed_residuals(df = data_mod_Lap, given = cond_site, v = q,a= discard(as.numeric(est_site$a),is.na),b = discard(as.numeric(est_site$b),is.na))
  # calculate distance from the conditioning site
  dist_tmp <- as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))
  # remove zero distance
  dist_tmp <- dist_tmp[dist_tmp>0]
  # normalise distance
  distnorm <- dist_tmp/max(dist_tmp)
  tmp <- sigmau_par_est_ite(data = Z, given = cond_site, cond_site_dist = distnorm, Nite = Nite, show_ite = TRUE, mu_init = discard(as.numeric(tmpsf$mu_agg_ite),is.na), sigl_init = discard(as.numeric(tmpsf$sigl_ite),is.na), sigu_init = discard(as.numeric(tmpsf$sigu_ite),is.na), deltal = as.numeric(tmpsf$deltal_ite[1]), deltau = as.numeric(tmpsf$deltau_ite[1]))
  
  # explore estimates
  # plot delta estimates
  tmp_phi <- rbind(data.frame("phi"=as.numeric(unlist(tmp[[4]])),iteration=1:(Nite+1),parameter = "phi0"),
                   data.frame("phi"=as.numeric(unlist(tmp[[5]])),iteration=1:(Nite+1),parameter = "delta_upper"))
  pphi <- ggplot(tmp_phi) + geom_point(aes(x=iteration,y=phi,col=parameter))
  ggsave(pd,filename=paste0("../Documents/iterative_sigmas_res_margin/deltas_",cond_site_name,".png"))
  
  # plot also for a random site as a check for convergence
  # plot across iterations for a selected site
  random_site <- 100
  tmp_all <- rbind(data.frame(delta=as.numeric(unlist(tmp[[1]][random_site,])),iteration=1:(Nite+1),parameter = "mu_agg"),
                   data.frame(delta=as.numeric(unlist(tmp[[2]][random_site,])),iteration=1:(Nite+1),parameter = "sigma_lower"),
                   data.frame(delta=as.numeric(unlist(tmp[[3]][random_site,])),iteration=1:(Nite+1),parameter = "sigma_upper"),
                   data.frame(delta=as.numeric(unlist(tmp[[4]])),iteration=1:(Nite+1),parameter = "delta_lower"),
                   data.frame(delta=as.numeric(unlist(tmp[[5]])),iteration=1:(Nite+1),parameter = "delta_upper")
                   
  )
  # plot the final iteration
  p1 <- ggplot(tmp_all) + geom_point(aes(x=iteration,y=delta,col=parameter))
  ggsave(p1,filename=paste0("../Documents/iterative_sigmas_res_margin/random_site_par_",cond_site_name,".png"))
  
  # explore also spatial parameters
  est_ite <- tmp[[6]] %>% add_row(.before=cond_site)
  names(est_ite) <- paste0(names(est_ite),"_ite_sigu")
  tmpsf <- cbind(tmpsf,est_ite)
  # plot parameter estimates against distance
  # calculate distance from a conditioning site with st_distance()
  mud <- data.frame(mu=tmpsf$mu_agg_ite_sigu,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=mu,x=dist))
  ggsave(mud,filename=paste0("../Documents/iterative_sigmas_res_margin/muagg_distance_",cond_site_name,".png")) 
  
  sigld <- data.frame(sigl=tmpsf$sigl_ite,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=sigl,x=dist))
  ggsave(sigld,filename=paste0("../Documents/iterative_sigmas_res_margin/sigl_distance_",cond_site_name,".png")) 
  
  sigud <- data.frame(sigu=tmpsf$sigu_ite,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=sigu,x=dist))
  ggsave(sigud,filename=paste0("../Documents/iterative_sigmas_res_margin/sigu_distance_",cond_site_name,".png")) 
  

  tmpsf <- tmpsf %>% mutate(mudiff = mu_agg_ite - mu_agg_ite_sigu, sigldiff = sigl_ite - sigl_ite_sigu, sigudiff = sigu_ite - sigu_ite_sigu)
  t <- tmpsf %>% dplyr::select(mu_agg_ite,mu_agg_ite_sigu,mudiff) %>% pivot_longer(cols=c(mu_agg_ite,mu_agg_ite_sigu,mudiff),names_to = "parameter", values_to = "value" ) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="-brewer.rd_bu")) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 12) 
  
  tmap_save(t,filename=paste0("../Documents/iterative_sigmas_res_margin/mu_agg_",cond_site_name,".png"),width=8,height=6)
  
  t <- tmpsf %>% dplyr::select(sigl_ite,sigl_ite_sigu,sigldiff) %>% pivot_longer(cols=c(sigl_ite,sigl_ite_sigu,sigldiff),names_to = "parameter", values_to = "value" ) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="brewer.blues")) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 12) 
  tmap_save(t,filename=paste0("../Documents/iterative_sigmas_res_margin/sigma_lower_",cond_site_name,".png"),width=8,height=6)
  
  t <- tmpsf %>% dplyr::select(sigu_ite,sigu_ite_sigu,sigudiff) %>% pivot_longer(cols=c(sigu_ite,sigu_ite_sigu,sigudiff),names_to = "parameter", values_to = "value" ) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="brewer.blues")) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 12) 
  tmap_save(t,filename=paste0("../Documents/iterative_sigmas_res_margin/sigma_upper_",cond_site_name,".png"),width=8,height=6)
  
  return(tmpsf)
}

iter_sigmau_site_wrapper <- function(i) {
  iter_sigmau_site(j=i)
}

result <- sapply(1:ncol(df_sites),iter_delta_site,simplify = FALSE)

deltal <- deltau <-  c()
for (i in 1:12) {
  deltal[i] <- result[[i]]$deltal_ite[1]
  deltau[i] <- result[[i]]$deltau_ite[1]
}


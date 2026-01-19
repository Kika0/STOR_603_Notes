# change functions to allow extra parameters phi0u (and optional also phi0l)

# start with the likelihood function
NLL_exp_phis <- function(phi,x = Z, d1j, mu1=as.numeric(unlist(mu_agg[,1])),deltal=2,deltau=2,phi0l=NULL,phi0u=NULL) {
  N <- nrow(x)
  mu1 <- rep(mu1, each = N)
  z <- as.numeric(unlist(x))
  dij. <- rep(d1j, each = N)
  log_lik <- rep(NA,length(z))
  if (is.null(deltal)==FALSE) {
    phi <- append(phi,deltal,after=4)
  } 
  if (is.null(deltau)==FALSE) {
    phi <- append(phi,deltau,after=5)
  } 
  if (is.null(phi0l)==FALSE) {
    phi <- append(phi,phi0l,after=6)
  }
  if (is.null(phi0u)==FALSE) {
    phi <- append(phi,phi0u,after=7)
  }
  deltal <- phi[5]
  deltau <- phi[6]
  phi0l <- phi[7]
  phi0u <- phi[8]
  if(deltal<1 |deltau<1 | phi[1]<0 | phi[2] < 0 | phi[3]<0 | phi[4]<0 | phi0l<0 | phi0u<0 ){return(10e10)}
  
  sigu <- phi0u + phi[1]*(1-exp(-(phi[2]*dij.)))
  sigl <- phi0l + phi[3]*(1-exp(-(phi[4]*dij.)))
  C_AGG <-  (sigl/deltal*gamma(1/deltal) + sigu/deltau*gamma(1/deltau)  )^(-1)
  log_lik[x<mu1] <- log(C_AGG[x<mu1])-((mu1[x<mu1]-x[x<mu1])/sigl[x<mu1])^deltal 
  log_lik[x>=mu1] <- log(C_AGG[x>=mu1])-((x[x>=mu1]-mu1[x>=mu1])/sigu[x>=mu1])^deltau 
  return(-sum(log_lik))
}

par_est_mu <- function(z,v=0.99,given=c(1),res_margin_est = res_margin_est) {
  lik <- mu_hat <- res_var <- c()
 # names(z) <- paste0("Z",res)
  d <- ncol(z)+1
  for (j in given) {
    res <- c(1:d)[-j]
    init_par <- c()
    init_lik <- c()
    sigl <- res_margin_est$sigl
    sigu <- res_margin_est$sigu
    deltal <- res_margin_est$deltal
    deltau <- res_margin_est$deltau
    for (i in 2:d) {
      # optimise using the initial parameters
      init_par <- c(0)
      opt <- optim(fn=NLL_AGG,x=as.numeric(unlist(z[,i-1])),sigl_hat = sigl[i-1], sigu_hat = sigu[i-1], deltal_hat = deltal, deltau_hat = deltau,par=init_par,control=list(maxit=2000))
      mu_hat <- append(mu_hat,opt$par[1])
      lik <- append(lik,opt$value)
      res_var <- append(res_var,res[i-1])
    }
  }
  par_sum <- data.frame("lik"=lik, 
                        "a" = a,"b" = b,
                        "mu" = mu_hat,
                        "sigl" = sigl,"sigu" = sigu,
                        "deltal" = deltal,"deltau" = deltau,
                        "given" = rep(given,each=(d-1)),"res" = res_var)  
  
  return(par_sum)
}

par_est_ite <- function(z=Z,v=q,given=cond_site,cond_site_dist, parest_site = result[[1]],Nite=10, show_ite=FALSE,deltal=NULL,deltau=NULL,phi0u=NULL,phi0l=NULL)  {
  d <- ncol(z)
  N <- nrow(z)
  res <- 1:d[-given]
  mu_agg <- sigl <- sigu <- data.frame(matrix(ncol=(Nite),nrow = d))
  
 phi1l. <- phi2l. <- phi1u. <- phi2u. <- c(1) 
 phi0l. <- phi0u. <- c(0)
  
  if (is.null(deltal)) {
    residual_pars <- list(sigl = parest_site$sigl_ite_sigl,
                          sigu = parest_site$sigu_ite_sigu,
                          deltal = parest_site$deltal_ite[1],
                          deltau = parest_site$deltau_ite[1])
    deltal. <- 2
    deltau. <- 2
    
  } else {
    residual_pars <- list("sigl" = parest_site$sigl_ite_sigl,
                          "sigu" = parest_site$sigu_ite_sigu,
                          "deltal" = deltal,
                          "deltau" = deltau)
    deltal. <- deltal
    deltau. <- deltau   
  }
  for (k in 1:Nite) {
    if (k >1) {
      residual_pars <- list(sigl=phi0l + phi1l *(1-exp(-phi2l*cond_site_dist)),sigu=phi0u + phi1u*(1-exp(-phi2u*cond_site_dist)),deltal=deltal,deltau=deltau)
    }
    # estimate mu
    pe <- par_est_mu(z=z,v=v,given=given,res_margin_est=residual_pars)
    # update mu parameters
    mu_agg[,k] <- pe$mu
    # estimate phi parameters
    phi_init <- c(phi0u.[k],phi1u.[k],phi2u.[k],phi1l.[k],phi2l.[k],deltal.[k],deltau.[k])
    opt <- optim(fn=NLL_exp_phis,x = Z,d1j=cond_site_dist,mu1=as.numeric(mu_agg[,k]),control=list(maxit=2000),par = phi_init,deltal=deltal,deltau=deltau,phi0l=0,method = "Nelder-Mead")
    phi1u <- opt$par[1]
    phi2u <- opt$par[2]
    phi1l <- opt$par[3]
    phi2l <- opt$par[4] 
    if (!is.numeric(deltal) & !is.numeric(deltau)) {
      deltal <- opt$par[5]
      deltau <- opt$par[6] } 
    if (!is.numeric(phi0u) & is.numeric(phi0l)) {
      phi0u <- opt$par[5]
    }
    if (!is.numeric(phi0u) & !is.numeric(phi0l)) {
      phi0l <- opt$par[5]
      phi0u <- opt$par[5]
    }
    phi0u. <- append(phi0u.,phi0u)
    phi1u. <- append(phi1u.,phi1u)
    phi2u. <- append(phi2u.,phi2u)
    phi0l. <- append(phi0l.,phi0l)
    phi1l. <- append(phi1l.,phi1l)
    phi2l. <- append(phi2l.,phi2l)
    deltal. <- append(deltal.,deltal)
    deltau. <- append(deltau.,deltau)
    sigu[,k] <- phi0u + phi1u*(1-exp(-phi2u*cond_site_dist))
    sigl[,k] <- phi0l + phi1l*(1-exp(-phi2l*cond_site_dist))
  }
  par_sum <- data.frame("mu_agg" = as.numeric(mu_agg[,Nite]),"sigl" = as.numeric(sigl[,Nite]),"sigu" = as.numeric(sigu[,Nite]),"phi1u" = phi1u.[Nite], "phi2u" = phi2u.[Nite],"phi1l" = phi1l.[Nite], "phi2l" = phi2l.[Nite], "deltal" = deltal.[Nite], "deltau" = deltau.[Nite])
  if (show_ite == TRUE) {
    return(list(mu_agg=mu_agg,sigl=sigl,sigu=sigu, phi0u=phi0u.,phi1u=phi1u., phi2u=phi2u.,phi0l=phi0l., phi1l=phi1l., phi2l=phi2l.,deltal=deltal., deltau=deltau.,par_sum))
  } else {return(par_sum)}
}


AGG_par_est_ite <- function(data_mod_Lap,site,Nite=10,sites = sites_index_diagonal,cond_site_names = site_name_diagonal,q=0.9,grid=xyUK20_sf,result,est_all_sf,deltal=NULL,deltau=NULL,folder_name=NULL) {
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
  # calculate observed residuals
  aest <- discard(est_all_sf %>% filter(cond_site==cond_site) %>% pull(a),is.na)
  best <- discard(est_all_sf %>% filter(cond_site==cond_site) %>% pull(b),is.na)
  Z <- observed_residuals(df=data_mod_Lap,given=cond_site,v = v,a=aest,b=best)
  if (is.null(deltal)) {
    try7 <- par_est_ite(z=Z,given=cond_site,cond_site_dist=distnorm, parest_site = parest_site,Nite=Nite,phi0l=0, show_ite=TRUE)
  } else {
    try7 <- par_est_ite(z=Z,given=cond_site,cond_site_dist=distnorm, parest_site = parest_site,Nite=Nite,phi0l=0, show_ite=TRUE,deltal=parest_site$deltal_ite[1],deltau= parest_site$deltau_ite[1]) 
    }

  if (is.null(folder_name)) {
    folder_name <- "mu_iterative" 
  }
  
  # separate parameter estimation and analysis
  p <- ggplot(try7[[1]] %>% pivot_longer(everything(),names_to = "iteration",values_to = "par") %>% mutate(iteration=factor(iteration,levels=paste0("X",1:Nite)))) + geom_boxplot(aes(x=iteration,y=par)) +ylab(TeX("$\\alpha"))
  ggsave(p,file=paste0("../Documents/",folder_name,"/boxplot_alpha_",cond_site_name,".png"),height=5,width=10)
  p <- ggplot(try7[[2]] %>% pivot_longer(everything(),names_to = "iteration",values_to = "par")%>% mutate(iteration=factor(iteration,levels=paste0("X",1:Nite)))) + geom_boxplot(aes(x=iteration,y=par)) +ylab(TeX("$\\beta"))
  ggsave(p,file=paste0("../Documents/",folder_name,"/boxplot_beta_",cond_site_name,".png"),height=5,width=10)
  p <- ggplot(try7[[3]] %>% pivot_longer(everything(),names_to = "iteration",values_to = "par")%>% mutate(iteration=factor(iteration,levels=paste0("X",1:Nite)))) + geom_boxplot(aes(x=iteration,y=par)) +ylab(TeX("$\\mu_{AGG}"))
  ggsave(p,file=paste0("../Documents/",folder_name,"/boxplot_mu_",cond_site_name,".png"),height=5,width=10)
  p <- ggplot(try7[[4]] %>% pivot_longer(everything(),names_to = "iteration",values_to = "par")%>% mutate(iteration=factor(iteration,levels=paste0("X",1:Nite)))) + geom_boxplot(aes(x=iteration,y=par)) +ylab(TeX("$\\sigma_l"))
  ggsave(p,file=paste0("../Documents/",folder_name,"/boxplot_sigmal_",cond_site_name,".png"),height=5,width=10)
  p <- ggplot(try7[[5]] %>% pivot_longer(everything(),names_to = "iteration",values_to = "par")%>% mutate(iteration=factor(iteration,levels=paste0("X",1:Nite)))) + geom_boxplot(aes(x=iteration,y=par)) +ylab(TeX("$\\sigma_u"))
  ggsave(p,file=paste0("../Documents/",folder_name,"/boxplot_sigmau_",cond_site_name,".png"),height=5,width=10)
  
  
  # look spatially to check
  # explore also spatial parameters
  est_ite <- try7[[12]] %>% add_row(.before=cond_site)
  names(est_ite) <- paste0(names(est_ite),"_ite")
  tmpsf <- cbind(est_all_sf %>% filter(cond_site==cond_site_name),est_ite)
  # plot parameter estimates against distance
  # calculate distance from a conditioning site with st_distance()
  
  mud <- data.frame(mu=tmpsf$mu_agg_ite,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=mu,x=dist))
  ggsave(mud,filename=paste0("../Documents/",folder_name,"/muagg_distance_",cond_site_name,".png")) 
  
  sigld <- data.frame(sigl=tmpsf$sigl_ite,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=sigl,x=dist))
  ggsave(sigld,filename=paste0("../Documents/",folder_name,"/sigl_distance_",cond_site_name,".png")) 
  
  sigud <- data.frame(sigu=tmpsf$sigu_ite,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=sigu,x=dist))
  ggsave(sigud,filename=paste0("../Documents/",folder_name,"/sigu_distance_",cond_site_name,".png")) 
  
  tmpsf <- tmpsf %>% mutate(mudiff = mu_agg_ite - mu_agg, sigldiff = sigl_ite - sigl, sigudiff = sigu_ite - sigu)
  toplabel <- c("New iterative approach","Original method","Difference")
  
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

file.sources = list.files(pattern="*helpers.R")
sapply(file.sources,source,.GlobalEnv)
#source("spatial_parameter_estimation.R") # for spatial_par_est function
load("data_processed/temperature_data.RData",verbose=TRUE) # for data_mod_Lap
load("data_processed/spatial_helper.RData",verbose = TRUE) # for xyUK20_sf
q <- 0.9
load("data_processed/N9000_sequential2_AGG_diagonal_sites_Birmingham_Cromer90.RData",verbose=TRUE)
est_all <- as.data.frame(est_all_sf)
# load previous iterative estimates for initial values
load("data_processed/iterative_sigmal_estimates_Birmingham_Cromer_diagonal.RData",verbose=TRUE)

# identify diagonal sites
# find indeces of start and end sites
Birmingham <- c(-1.9032,52.4806)
Cromer <- c(1.28486,53.05349)
site_start <- find_site_index(site=Birmingham,grid_uk = xyUK20_sf)
site_end <- find_site_index(site=Cromer,grid_uk = xyUK20_sf)
result <- sapply(1:1,FUN = AGG_par_est_ite,data_mod_Lap = data_mod_Lap,sites = sites_index_diagonal,cond_site_names = site_name_diagonal,est_all_sf = est_all_sf,result=result,folder_name = "Birmingham_Cromer_diagonal/new_iterative_sigmas_mu",simplify = FALSE)


summary(result[[1]])

#' Estimate conditional model parameters for spatio-temporal data
#'
#' `spatial_par_est()` estimates $\alpha$ and $\beta$ dependence parameters,
#'  residual margin parameters and links back to spatial locations
#'
#' @param data_Lap A dataframe of observed data on Laplace scale Y
#' @param cond_sites A dataframe of conditioning sites (a subset of df_sites)
#' @param dayshift A numeric vector of temporal lags (days)
#' @param v A quantile
#' @param ab_method A method to estimate dependence parameters $\alpha$ and $\beta$
#' @param res_margin A parametric method to estimate residual margins
#' @param grid_uk An sf point object of site locations
#'
#' @return `spatial_par_est()` returns an sf object of spatial points with parameter estimates
#' @export
#'
#' @examples
spatial_par_est <- function(data_Lap,cond_sites,cond_site_names=NULL,dayshift=c(0),Ndays_season=90,v=0.9,ab_method="sequential2",res_margin="AGG",grid_uk=xyUK20_sf,title="") {
  est_all <- data.frame("lik" = numeric(), "lika"= numeric(),"likb"= numeric(),"lik2"= numeric(),
                            "a" = numeric(), "b" = numeric(),
                            "mu" = numeric(),"mu_agg" = numeric(),
                            "sig" = numeric(),"sig_agg" = numeric(),"sigl" = numeric(),"sigu" = numeric(),
                            "delta" = numeric(),"deltal" = numeric(), "deltau" = numeric(),
                            "given" = numeric(), "res" = numeric(), "cond_site" = character(), "tau" = numeric())
  if(is.null(cond_site_names)) {
    cond_site_names <- names(cond_sites)
  }
  for (i in 1:length(cond_site_names)) {
    if (is.numeric(cond_sites)) {cond_site <- cond_sites[i]} else{
      cond_site <- find_site_index(as.numeric(cond_sites[,i]),grid_uk = grid_uk)
    }
    for (j in 1:length(dayshift)) {
      sims_tau <- shift_time(sims=data_Lap,cond_site=cond_site,tau=dayshift[j],Ndays_season = Ndays_season)
      pe <- par_est(df=sims_tau,v=v,given=cond_site,margin = "Normal", method=ab_method,keef_constraints = c(1,2))
      # calculate observed residuals
      obsr <- observed_residuals(df=sims_tau,given=cond_site,v = v,a=pe$a,b=pe$b)
      # estimate residual margin parameters
      pe_res <- res_margin_par_est(obs_res = obsr,method="AGG")
      est_all <- rbind(est_all,cbind(pe[,-c(8,10:12,14:15)],pe_res) %>% add_row(.before=cond_site) %>%  mutate("cond_site"=cond_site_names[i],tau=as.character(dayshift[j])))
       }
  }
  est_all <- est_all %>% mutate(tau=factor(as.character(tau),levels = as.character(dayshift))) %>% mutate(cond_site=factor(cond_site))
  est_all_sf <- est_join_spatial(tmp_est=est_all,grid_uk=grid_uk)
  save(est_all_sf,file=paste0("data_processed/N",nrow(data_Lap),"_",ab_method,"_",res_margin,"_",title,".RData"))
}

#' Iterative estimation for common deltal and deltau parameters
#'
#' @param data A dataframe of observed residuals
#' @param given An index of a conditioning variable
#' @param N A number of iterations
#' @param show_ite A logical to return also values from each iteration
#' @param mu_init A vector of initial value for mu
#' @param sigl_init A vector of initial value for sigl
#' @param sigu_init A vector of initial value for sigu
#' @param deltal_init A vector of initial value for deltal
#' @param deltau_init A vector of initial value for deltau
#'
#' @return A dataframe of final iteration estimates or a list for each parameter when show_ite=TRUE
#' @export
#'
#' @examples
deltas_par_est_ite <- function(data=Z,given=cond_site,N=10, show_ite=FALSE,mu_init=NULL,sigl_init=NULL,sigu_init=NULL,deltal_init=NULL,deltau_init=NULL)  {
  d <- ncol(data)
  nv <- nrow(data)
  res <- 1:d
  mu_agg <- sigl <- sigu <- deltal <- deltau <- data.frame(matrix(ncol=(N+1),nrow = d))
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
  if (is.numeric(deltal_init) & is.numeric(deltau_init)) {
    deltal[,1] <- deltal_init
    deltau[,1] <- deltau_init
  } else {
    deltal[,1] <- 2
    deltau[,1] <- 2
  }
  
  for (i in 1:N) {
    for (j in 1:d) {
      Z2 <- as.numeric(unlist(data[,j]))
      opt <- optim(fn=NLL_AGG,x=Z2,deltal_hat=as.numeric(deltal[1,i]),deltau_hat=as.numeric(deltau[1,i]),par=c(mean(Z2),sd(Z2),sd(Z2)),control=list(maxit=2000),method = "Nelder-Mead")
      mu_agg[j,i+1] <- opt$par[1]
      sigl[j,i+1] <- opt$par[2]
      sigu[j,i+1] <- opt$par[3]
    }
    # calculate deltal and deltau
    opt <- optim(fn=NLL_AGG_deltas,x = data,mu1=as.numeric(mu_agg[,i+1]),sigl1=as.numeric(sigl[,i+1]),sigu1=as.numeric(sigu[,i+1]),control=list(maxit=2000),par = c(as.numeric(deltal[1,i]),as.numeric(deltau[1,i])),method = "Nelder-Mead")
    deltal[,i+1] <- opt$par[1]
    deltau[,i+1] <- opt$par[2]
  }
  par_sum <- data.frame("mu_agg" = as.numeric(mu_agg[,N+1]),"sigl" = as.numeric(sigl[,N+1]),"sigu" = as.numeric(sigu[,N+1]),"deltal" = as.numeric(deltal[,N+1]), "deltau" = as.numeric(deltau[,N+1]))
  if (show_ite == TRUE) {
    return(list(mu_agg,sigl,sigu,deltal,deltau,par_sum))
  } else {return(par_sum)}
}

# write into a function ----------------------------------------------------
iter_delta_site <- function(j,Nite=50,sites=df_sites,cond_site_names=NULL,grid=xyUK20_sf,data=data_mod_Lap,par_est=est_all_sf,folder_name = "iterative_deltas_res_margin") {
  q <- 0.9
  if(is.null(cond_site_names)) {
    cond_site_name <- names(sites)[j]
  } else {  cond_site_name <- cond_site_names[j] }

    if (is.numeric(sites)) {cond_site <- sites[j]} else{
      cond_site_coord <- sites %>% dplyr::select(all_of(cond_site_name)) %>% pull()
      cond_site <- find_site_index(cond_site_coord,grid_uk = grid)    }
 
  est_site <- par_est %>% filter(cond_site==cond_site_name)
  Z <- observed_residuals(df = data_mod_Lap, given = cond_site, v = q,a= discard(as.numeric(est_site$a),is.na),b = discard(as.numeric(est_site$b),is.na))
  
  # estimate parameters iteratively --------------------------------------------
  tmp <- deltas_par_est_ite(data=Z,show_ite = TRUE,N=Nite, mu_init = discard(as.numeric(est_site$mu_agg),is.na), sigl_init = discard(as.numeric(est_site$sigl),is.na), sigu_init = discard(as.numeric(est_site$sigu),is.na), deltal = discard(as.numeric(est_site$deltal),is.na), deltau = discard(as.numeric(est_site$deltau),is.na))
  # plot delta estimates
  tmp_delta <- rbind(data.frame(delta=as.numeric(unlist(tmp[[4]][1,])),iteration=1:(Nite+1),parameter = "delta_lower"),
                     data.frame(delta=as.numeric(unlist(tmp[[5]][1,])),iteration=1:(Nite+1),parameter = "delta_upper"))
  pd <- ggplot(tmp_delta) + geom_point(aes(x=iteration,y=delta,col=parameter))
  ggsave(pd,filename=paste0("../Documents/",folder_name,"/deltas_",cond_site_name,".png"))
  
  # plot across iterations for a selected site
  random_site <- 100
  tmp_all <- rbind(data.frame(delta=as.numeric(unlist(tmp[[1]][random_site,])),iteration=1:(Nite+1),parameter = "mu_agg"),
                   data.frame(delta=as.numeric(unlist(tmp[[2]][random_site,])),iteration=1:(Nite+1),parameter = "sigma_lower"),
                   data.frame(delta=as.numeric(unlist(tmp[[3]][random_site,])),iteration=1:(Nite+1),parameter = "sigma_upper"),
                   data.frame(delta=as.numeric(unlist(tmp[[4]][random_site,])),iteration=1:(Nite+1),parameter = "delta_lower"),
                   data.frame(delta=as.numeric(unlist(tmp[[5]][random_site,])),iteration=1:(Nite+1),parameter = "delta_upper")
                   
  )
  # plot the final iteration
  p1 <- ggplot(tmp_all) + geom_point(aes(x=iteration,y=delta,col=parameter))
  ggsave(p1,filename=paste0("../Documents/",folder_name,"/random_site_par_",cond_site_name,".png"))
  
  est_ite <- tmp[[6]] %>% add_row(.before=cond_site)
  names(est_ite) <- paste0(names(est_ite),"_ite")
  tmpsf <- cbind(est_site,est_ite)
  # plot parameter estimates against distance
  # calculate distance from a conditioning site with st_distance()
  mud <- data.frame(mu=tmpsf$mu_agg_ite,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=mu,x=dist))
  ggsave(mud,filename=paste0("../Documents/",folder_name,"/muagg_distance_",cond_site_name,".png")) 
  
  sigld <- data.frame(sigl=tmpsf$sigl_ite,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=sigl,x=dist))
  ggsave(sigld,filename=paste0("../Documents/",folder_name,"/sigl_distance_",cond_site_name,".png")) 
  
  sigud <- data.frame(sigu=tmpsf$sigu_ite,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=sigu,x=dist))
  ggsave(sigud,filename=paste0("../Documents/",folder_name,"/sigu_distance_",cond_site_name,".png")) 
  
  # map final iteration 
  t <- tmpsf %>% dplyr::select(mu_agg,mu_agg_ite) %>% pivot_longer(cols=c(mu_agg,mu_agg_ite),names_to = "parameter", values_to = "value" ) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="-brewer.rd_bu")) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 12) 
  
  tmap_save(t,filename=paste0("../Documents/",folder_name,"/mu_agg_",cond_site_name,".png"),width=8,height=6)
  
  t <- tmpsf %>% dplyr::select(sigl,sigl_ite) %>% pivot_longer(cols=c(sigl,sigl_ite),names_to = "parameter", values_to = "value" ) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="brewer.blues")) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 12) 
  tmap_save(t,filename=paste0("../Documents/",folder_name,"/sigma_lower_",cond_site_name,".png"),width=8,height=6)
  
  t <- tmpsf %>% dplyr::select(sigu,sigu_ite) %>% pivot_longer(cols=c(sigu,sigu_ite),names_to = "parameter", values_to = "value" ) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="brewer.blues")) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 12) 
  tmap_save(t,filename=paste0("../Documents/",folder_name,"/sigma_upper_",cond_site_name,".png"),width=8,height=6)
  
  return(tmpsf)
}

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

iter_sigmau_site <- function(j,Nite=5,sites=df_sites,cond_site_names=NULL,grid=xyUK20_sf,data=data_mod_Lap,par_est=est_all_sf,ite_delta = result,index_outliers = NULL, folder_name = "iterative_sigmau_res_margin") {
  q <- 0.9
  if(is.null(cond_site_names)) {
    cond_site_name <- names(sites)[j]
    cond_site_names <- names(sites)
  } else {  cond_site_name <- cond_site_names[j] }
  
  if (is.numeric(sites)) {cond_site <- sites[j]} else{
    cond_site_coord <- sites %>% dplyr::select(all_of(cond_site_name)) %>% pull()
    cond_site <- find_site_index(cond_site_coord,grid_uk = grid)    }
  
  if (is.null(index_outliers)) {
    est_site <- par_est %>% filter(cond_site==cond_site_name)
    tmpsf <- ite_delta[[ which(cond_site_names==cond_site_name) ]]
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
  ggsave(pphi,filename=paste0("../Documents/",folder_name,"/phi_",cond_site_name,".png"),width=5,height=5)
  
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
  ggsave(p1,filename=paste0("../Documents/",folder_name,"/random_site_par_",cond_site_name,".png"),width=5,height=5)
  
  # explore also spatial parameters
  est_ite <- tmp[[6]] %>% add_row(.before=cond_site)
  names(est_ite) <- paste0(names(est_ite),"_ite_sigu")
  tmpsf <- cbind(tmpsf,est_ite)
  # plot parameter estimates against distance
  # calculate distance from a conditioning site with st_distance()
  mud <- data.frame(mu=tmpsf$mu_agg_ite_sigu,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=mu,x=dist))
  ggsave(mud,filename=paste0("../Documents/",folder_name,"/muagg_distance_",cond_site_name,".png")) 
  
  sigld <- data.frame(sigl=tmpsf$sigl_ite_sigu,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=sigl,x=dist))
  ggsave(sigld,filename=paste0("../Documents/",folder_name,"/sigl_distance_",cond_site_name,".png")) 
  
  sigud <- data.frame(sigu=tmpsf$sigu_ite_sigu,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=sigu,x=dist))
  ggsave(sigud,filename=paste0("../Documents/",folder_name,"/sigu_distance_",cond_site_name,".png")) 
  
  
  tmpsf <- tmpsf %>% mutate(mudiff = mu_agg_ite - mu_agg_ite_sigu, sigldiff = sigl_ite - sigl_ite_sigu, sigudiff = sigu_ite - sigu_ite_sigu)
  toplabel <- c(TeX("Iterative $\\delta$s method"),TeX("Iterative $\\sigma_u (d_j)$ method"),"Difference")
  mu_limits <- c(-1.75,1.79)
  t <- tmpsf %>% dplyr::select(mu_agg_ite,mu_agg_ite_sigu,mudiff) %>% pivot_longer(cols=c(mu_agg_ite,mu_agg_ite_sigu,mudiff),names_to = "parameter", values_to = "value" ) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="-brewer.rd_bu",limits=mu_limits),fill.legend = tm_legend(title = TeX("$\\mu_{AGG}$"))) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 12, panel.labels = toplabel) 
  
  tmap_save(t,filename=paste0("../Documents/",folder_name,"/mu_agg_",cond_site_name,".png"),width=8,height=6)
  
  sigma_limits <- c(-1.32,2.3)
  t <- tmpsf %>% dplyr::select(sigl_ite,sigl_ite_sigu,sigldiff) %>% pivot_longer(cols=c(sigl_ite,sigl_ite_sigu,sigldiff),names_to = "parameter", values_to = "value" ) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="-brewer.rd_bu",limits=sigma_limits),fill.legend = tm_legend(title = TeX("$\\sigma_l$"))) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 9, panel.labels = toplabel) 
  tmap_save(t,filename=paste0("../Documents/",folder_name,"/sigma_lower_",cond_site_name,".png"),width=8,height=6)
  
  t <- tmpsf %>% dplyr::select(sigu_ite,sigu_ite_sigu,sigudiff) %>% pivot_longer(cols=c(sigu_ite,sigu_ite_sigu,sigudiff),names_to = "parameter", values_to = "value" ) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="-brewer.rd_bu",limits=sigma_limits),fill.legend = tm_legend(title = TeX("$\\sigma_u$"))) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 9, panel.labels = toplabel) 
  tmap_save(t,filename=paste0("../Documents/",folder_name,"/sigma_upper_",cond_site_name,".png"),width=8,height=6)
  
  return(tmpsf)
}


sigmal_par_est_ite <- function(data=Z,given=cond_site,cond_site_dist, Nite=10, show_ite=FALSE,mu_init=NULL,sigl_init=NULL,deltal,deltau,sigu)  {
  d <- ncol(data)
  N <- nrow(data)
  res <- 1:d
  mu_agg <- sigl <- data.frame(matrix(ncol=(Nite+1),nrow = d))
  # calculate a with initial values for mu and sigma
  if (is.numeric(mu_init) & is.numeric(sigl_init) ) {
    mu_agg[,1] <- mu_init
    sigl[,1] <- sigl_init
  } else {
    mu_agg[,1] <- 0
    sigl[,1] <- 1
  }
  phi2. <- phi3. <- c(1) 
  for (i in 1:Nite) {
    # estimate sigu parameters phi0 and phi1
    phi_init <- c(phi2.[i],phi3.[i])
    opt <- optim(fn=NLL_exp_sigmal,x = data,d1j=cond_site_dist,mu1=as.numeric(mu_agg[,i]),sigu1=sigu,deltal=deltal,deltau=deltau,control=list(maxit=2000),par = phi_init,method = "Nelder-Mead")
    phi2 <- opt$par[1]
    phi3 <- opt$par[2]
    phi2. <- append(phi2.,phi2)
    phi3. <- append(phi3.,phi3)
    sigl[,i+1] <- phi2*(1-exp(-phi3*cond_site_dist))
    # estimate sigl and mu_agg for each site separately
    for (j in 1:d) {
      Z2 <- as.numeric(unlist(data[,j]))
      opt <- optim(fn=NLL_AGG,x=Z2,sigl_hat = as.numeric(unlist(sigl[,i+1])),sigu_hat = sigu,deltal_hat=deltal,deltau_hat = deltau,par=c(mean(Z2),sd(Z2)),control=list(maxit=2000),method = "Nelder-Mead")
      mu_agg[j,i+1] <- opt$par[1]
    }
    
  }
  par_sum <- data.frame("mu_agg" = as.numeric(mu_agg[,Nite+1]),"sigl" = as.numeric(sigl[,Nite+1]),"phi2" = phi2.[Nite+1], "phi3" = phi3.[Nite+1])
  if (show_ite == TRUE) {
    return(list(mu_agg,sigl,phi2.,phi3.,par_sum))
  } else {return(par_sum)}
}

# repeat for all other sites -------------------------------------------------
iter_sigmal_site <- function(j,Nite=5,sites=df_sites,cond_site_names=NULL,grid=xyUK20_sf,data=data_mod_Lap,par_est=est_all_sf,ite_sigu = iterative_sigmau_estimates,index_outliers = NULL,folder_name="iterative_sigmal_res_margin") {
  q <- 0.9
  if(is.null(cond_site_names)) {
    cond_site_name <- names(sites)[j]
    cond_site_names <- names(sites)
  } else {  cond_site_name <- cond_site_names[j] }
  
  if (is.numeric(sites)) {cond_site <- sites[j]} else{
    cond_site_coord <- sites %>% dplyr::select(all_of(cond_site_name)) %>% pull()
    cond_site <- find_site_index(cond_site_coord,grid_uk = grid)    }
  
  
  if (is.null(index_outliers)) {
    est_site <- par_est %>% filter(cond_site==cond_site_name)
    tmpsf <- ite_sigu[[ which(cond_site_names==cond_site_name) ]]
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
  ggsave(pphi,filename=paste0("../Documents/",folder_name,"/phi_",cond_site_name,".png"),width=5,height=5)
  
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
  ggsave(p1,filename=paste0("../Documents/",folder_name,"/random_site_par_",cond_site_name,".png"),width=5,height=5)
  
  # explore also spatial parameters
  est_ite <- tmp[[5]] %>% add_row(.before=cond_site)
  names(est_ite) <- paste0(names(est_ite),"_ite_sigl")
  tmpsf <- cbind(tmpsf,est_ite)
  # plot parameter estimates against distance
  # calculate distance from a conditioning site with st_distance()
  mud <- data.frame(mu=tmpsf$mu_agg_ite_sigl,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=mu,x=dist))
  ggsave(mud,filename=paste0("../Documents/",folder_name,"/muagg_distance_",cond_site_name,".png")) 
  
  sigld <- data.frame(sigl=tmpsf$sigl_ite_sigl,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=sigl,x=dist))
  ggsave(sigld,filename=paste0("../Documents/",folder_name,"/sigl_distance_",cond_site_name,".png")) 
  
  
  tmpsf <- tmpsf %>% mutate(mudiff = mu_agg_ite_sigu - mu_agg_ite_sigl, sigldiff = sigl_ite_sigl - sigl_ite_sigu)
  toplabel <- c(TeX("Iterative $\\sigma_u (d_j)$ method"),TeX("Iterative $\\sigma_l (d_j)$ method"),"Difference")
  mu_limits <- c(-1.75,1.79)
  t <- tmpsf %>% dplyr::select(mu_agg_ite_sigu,mu_agg_ite_sigl,mudiff) %>% pivot_longer(cols=c(mu_agg_ite_sigu,mu_agg_ite_sigl,mudiff),names_to = "parameter", values_to = "value") %>% mutate(parameter=factor(parameter,levels=c("mu_agg_ite_sigu","mu_agg_ite_sigl","mudiff"))) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="-brewer.rd_bu",limits=mu_limits),fill.legend = tm_legend(title = TeX("$\\mu_{AGG}$"))) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 12, panel.labels = toplabel) 
  
  tmap_save(t,filename=paste0("../Documents/",folder_name,"/mu_agg_",cond_site_name,".png"),width=8,height=6)
  
  sigma_limits <- c(-1.32,2.3)
  t <- tmpsf %>% dplyr::select(sigl_ite_sigu,sigl_ite_sigl,sigldiff) %>% pivot_longer(cols=c(sigl_ite_sigu,sigl_ite_sigl,sigldiff),names_to = "parameter", values_to = "value" ) %>% mutate(parameter=factor(parameter,levels=c("sigl_ite_sigu","sigl_ite_sigl","sigldiff"))) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="-brewer.rd_bu",limits=sigma_limits),fill.legend = tm_legend(title = TeX("$\\sigma_l$"))) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 9, panel.labels = toplabel) 
  tmap_save(t,filename=paste0("../Documents/",folder_name,"/sigma_lower_",cond_site_name,".png"),width=8,height=6)
  
  return(tmpsf)
}


#' Get index of a conditioning site
#'
#' @param j Column index of sites
#' @param grid An sf object of site points grid
#' @param sites A dataframe of conditioning sites df_sites
#'
#' @return An index on the grid of a conditioning site j
#' @export
#'
#' @examples
get_site_index <- function(j, grid = xyUK20_sf, sites= df_sites) {
  cond_site_name <- names(sites)[j]  
  cond_site_coord <- sites %>% dplyr::select(all_of(cond_site_name)) %>% pull()
  return(find_site_index(cond_site_coord,grid_uk = grid))
}

#' Estimate conditional model parameters for spatio-temporal data
#'
#' `spatial_par_est_abmu()` estimates $\alpha$, $\beta$ and $\mu$ parameters,
#'  and links back to spatial locations
#'
#' @param data_Lap A dataframe of data on Laplace scale Y
#' @param cond_sites A dataframe of conditioning sites ((a subset of) df_sites)
#' @param dayshift A numeric vector of temporal lags (days)
#' @param Ndays_season A number of days in a season (90 for CPM data)
#' @param v A quantile threshold
#' @param res_margin_est A list of residual parameter estimates sigl,sigu,deltal,deltau
#' @param grid_uk An sf point object of site locations
#' @param title A character string added to end of filename of .RData
#'
#' @return `spatial_par_est()` returns an sf object of spatial points with parameter estimates
#' @export
#'
#' @examples
spatial_par_est_abmu <- function(data_Lap,cond_sites,dayshift=c(0),Ndays_season=90,v=0.9,res_margin_est,grid_uk=xyUK20_sf,title="") {
  est_all <- data.frame("lik" = numeric(), 
                        "a" = numeric(), "b" = numeric(),
                        "mu" = numeric(),
                        "given" = numeric(), "res" = numeric(), "cond_site" = character(), "tau" = numeric())
  for (i in 1:ncol(cond_sites)) {
    cond_site <- find_site_index(as.numeric(cond_sites[,i]),grid_uk = grid_uk)
    parest_site <- st_drop_geometry(res_margin_est[[i]]) %>% dplyr::select(sigl_ite_sigl,sigu_ite_sigu,deltal_ite,deltau_ite) %>% na.omit()
    residual_pars <- list(sigl = parest_site$sigl_ite_sigl,
                          sigu = parest_site$sigu_ite_sigu,
                          deltal = parest_site$deltal_ite[1],
                          deltau = parest_site$deltau_ite[1])
    
    for (j in 1:length(dayshift)) {
      sims_tau <- shift_time(sims=data_Lap,cond_site=cond_site,tau=dayshift[j],Ndays_season = Ndays_season)
      pe <- par_est_abmu(df=sims_tau,v=q,given=cond_site,res_margin_est = residual_pars)
      # estimate residual margin parameters
      est_all <- rbind(est_all,pe %>% add_row(.before=cond_site) %>%  mutate(cond_site=names(cond_sites[i]),tau=as.character(dayshift[j])))
    }
  }
  est_all <- est_all %>% mutate(tau=factor(as.character(tau),levels = as.character(dayshift))) %>% mutate(cond_site=factor(cond_site))
  est_all_sf <- est_join_spatial(tmp_est=est_all,grid_uk=grid_uk)
  save(est_all_sf,file=paste0("data_processed/N",nrow(data_Lap),"_",title,".RData"))
}


sigmas_par_est_ite <- function(data=Z,given=cond_site,cond_site_dist, Nite=10, show_ite=FALSE,mu_init=NULL,sigl_init=NULL,sigu_init=NULL,deltal,deltau)  {
  d <- ncol(data)
  N <- nrow(data)
  res <- 1:d
  mu_agg <- sigl <- sigu <- data.frame(matrix(ncol=(Nite+1),nrow = d))
  # calculate a with initial values for mu and sigma
  if (is.numeric(mu_init) & is.numeric(sigl_init) & is.numeric(sigu_init)) {
    mu_agg[,1] <- mu_init
    sigl[,1] <- sigl_init
    sigu[,1] <- sigu_init
  } else {
    mu_agg[,1] <- 0
    sigl[,1] <- 1
    sigu[,1] <- 1
  }
  phi0. <- phi1. <- phi2. <- phi3. <- c(1) 
  for (i in 1:Nite) {
    # estimate sigu parameters phi0 and phi1
    phi_init <- c(phi0.[i],phi1.[i],phi2.[i],phi3.[i])
    opt <- optim(fn=NLL_exp_phis,x = data,d1j=cond_site_dist,mu1=as.numeric(mu_agg[,i]),deltal=deltal,deltau=deltau,control=list(maxit=2000),par = phi_init,method = "Nelder-Mead")
    phi0 <- opt$par[1]
    phi1 <- opt$par[2]
    phi2 <- opt$par[3]
    phi3 <- opt$par[4]    
    phi0. <- append(phi0.,phi0)
    phi1. <- append(phi1.,phi1)
    phi2. <- append(phi2.,phi2)
    phi3. <- append(phi3.,phi3)
    sigu[,i+1] <- phi0*(1-exp(-phi1*cond_site_dist))
    sigl[,i+1] <- phi2*(1-exp(-phi3*cond_site_dist))
    # estimate sigl and mu_agg for each site separately
    for (j in 1:d) {
      Z2 <- as.numeric(unlist(data[,j]))
      opt <- optim(fn=NLL_AGG,x=Z2,sigl_hat = as.numeric(unlist(sigl[,i+1])),sigu_hat = as.numeric(unlist(sigu[,i+1])),deltal_hat=deltal,deltau_hat = deltau,par=c(mean(Z2),sd(Z2),sd(Z2)),control=list(maxit=2000),method = "Nelder-Mead")
      mu_agg[j,i+1] <- opt$par[1]
    }
    
  }
  par_sum <- data.frame("mu_agg" = as.numeric(mu_agg[,Nite+1]),"sigl" = as.numeric(sigl[,Nite+1]),"sigu" = as.numeric(sigu[,Nite+1]),"phi0" = phi0.[Nite+1], "phi1" = phi1.[Nite+1],"phi2" = phi2.[Nite+1], "phi3" = phi3.[Nite+1])
  if (show_ite == TRUE) {
    return(list(mu_agg,sigl,sigu,phi0.,phi1.,phi2.,phi3.,par_sum))
  } else {return(par_sum)}
}

iter_sigmas_site <- function(j,Nite=5,sites=df_sites,cond_site_names=NULL,grid=xyUK20_sf,data=data_mod_Lap,par_est=est_all_sf,ite_delta = result,index_outliers = NULL, folder_name = "iterative_sigmas_res_margin") {
  q <- 0.9
  if(is.null(cond_site_names)) {
    cond_site_name <- names(sites)[j]
    cond_site_names <- names(sites)
  } else {  cond_site_name <- cond_site_names[j] }
  
  if (is.numeric(sites)) {cond_site <- sites[j]} else{
    cond_site_coord <- sites %>% dplyr::select(all_of(cond_site_name)) %>% pull()
    cond_site <- find_site_index(cond_site_coord,grid_uk = grid)    }
  
  if (is.null(index_outliers)) {
    est_site <- par_est %>% filter(cond_site==cond_site_name)
    tmpsf <- ite_delta[[ which(cond_site_names==cond_site_name) ]]
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
  
  Z <- observed_residuals(df = data_mod_Lap, given = cond_site, v = q,a= discard(as.numeric(est_site$a),is.na),b = discard(as.numeric(est_site$b),is.na))
  # calculate distance from the conditioning site
  dist_tmp <- as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))
  # remove zero distance
  dist_tmp <- dist_tmp[dist_tmp>0]
  # normalise distance using a common constant
  distnorm <- dist_tmp/1000000
  tmp <- sigmas_par_est_ite(data = Z, given = cond_site, cond_site_dist = distnorm, Nite = Nite, show_ite = TRUE, mu_init = discard(as.numeric(tmpsf$mu_agg_ite),is.na), sigl_init = discard(as.numeric(tmpsf$sigl_ite),is.na), sigu_init = discard(as.numeric(tmpsf$sigu_ite),is.na), deltal = as.numeric(tmpsf$deltal_ite[1]), deltau = as.numeric(tmpsf$deltau_ite[1]))
  
  # explore estimates
  # plot phi estimates
  tmp_phi <- rbind(data.frame("phi"=as.numeric(unlist(tmp[[4]])),iteration=1:(Nite+1),parameter = "phi0"),
                   data.frame("phi"=as.numeric(unlist(tmp[[5]])),iteration=1:(Nite+1),parameter = "phi1"),
                   data.frame("phi"=as.numeric(unlist(tmp[[6]])),iteration=1:(Nite+1),parameter = "phi2"),
                   data.frame("phi"=as.numeric(unlist(tmp[[7]])),iteration=1:(Nite+1),parameter = "phi3"))
  pphi <- ggplot(tmp_phi) + geom_point(aes(x=iteration,y=phi,col=parameter)) + scale_color_manual(values = c("#009ADA","#C11432","#66A64F","#FDD10A"), breaks = c("phi0","phi1","phi2","phi3"),labels = c(TeX("$\\phi_0$"),TeX("$\\phi_1$"),TeX("$\\phi_2$"),TeX("$\\phi_3$"))) + ylab("")
  ggsave(pphi,filename=paste0("../Documents/",folder_name,"/phi_",cond_site_name,".png"),width=5,height=5)
  
  # plot also for a random site as a check for convergence
  # plot across iterations for a selected site
  random_site <- 100
  tmp_all <- rbind(data.frame(delta=as.numeric(unlist(tmp[[1]][random_site,])),iteration=1:(Nite+1),parameter = "mu_agg"),
                   data.frame(delta=as.numeric(unlist(tmp[[2]][random_site,])),iteration=1:(Nite+1),parameter = "sigma_lower"),
                   data.frame(delta=as.numeric(unlist(tmp[[3]][random_site,])),iteration=1:(Nite+1),parameter = "sigma_upper"),
                   data.frame(delta=as.numeric(unlist(tmp[[4]])),iteration=1:(Nite+1),parameter = "phi0"),
                   data.frame(delta=as.numeric(unlist(tmp[[5]])),iteration=1:(Nite+1),parameter = "phi1"),
                   data.frame(delta=as.numeric(unlist(tmp[[6]])),iteration=1:(Nite+1),parameter = "phi2"),
                   data.frame(delta=as.numeric(unlist(tmp[[7]])),iteration=1:(Nite+1),parameter = "phi3")
                   
  )
  # plot the final iteration
  p1 <- ggplot(tmp_all) + geom_point(aes(x=iteration,y=delta,col=parameter)) + ylab("")
  ggsave(p1,filename=paste0("../Documents/",folder_name,"/random_site_par_",cond_site_name,".png"),width=5,height=5)
  
  # explore also spatial parameters
  est_ite <- tmp[[8]] %>% add_row(.before=cond_site)
  names(est_ite) <- paste0(names(est_ite),"_ite_sig")
  tmpsf <- cbind(tmpsf,est_ite)
  # plot parameter estimates against distance
  # calculate distance from a conditioning site with st_distance()
  mud <- data.frame(mu=tmpsf$mu_agg_ite_sig,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=mu,x=dist))
  ggsave(mud,filename=paste0("../Documents/",folder_name,"/muagg_distance_",cond_site_name,".png")) 
  
  sigld <- data.frame(sigl=tmpsf$sigl_ite_sig,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=sigl,x=dist))
  ggsave(sigld,filename=paste0("../Documents/",folder_name,"/sigl_distance_",cond_site_name,".png")) 
  
  sigud <- data.frame(sigu=tmpsf$sigu_ite_sig,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=sigu,x=dist))
  ggsave(sigud,filename=paste0("../Documents/",folder_name,"/sigu_distance_",cond_site_name,".png")) 
  
  
  tmpsf <- tmpsf %>% mutate(mudiff = mu_agg_ite - mu_agg_ite_sig, sigldiff = sigl_ite - sigl_ite_sig, sigudiff = sigu_ite - sigu_ite_sig)
  toplabel <- c(TeX("Iterative $\\delta$s method"),TeX("Iterative $\\phi$s method"),"Difference")
  mu_limits <- c(-1.75,1.79)
  t <- tmpsf %>% dplyr::select(mu_agg_ite,mu_agg_ite_sig,mudiff) %>% pivot_longer(cols=c(mu_agg_ite,mu_agg_ite_sig,mudiff),names_to = "parameter", values_to = "value" ) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="-brewer.rd_bu",limits=mu_limits),fill.legend = tm_legend(title = TeX("$\\mu_{AGG}$"))) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 12, panel.labels = toplabel) 
  
  tmap_save(t,filename=paste0("../Documents/",folder_name,"/mu_agg_",cond_site_name,".png"),width=8,height=6)
  
  sigma_limits <- c(-1.32,2.3)
  t <- tmpsf %>% dplyr::select(sigl_ite,sigl_ite_sig,sigldiff) %>% pivot_longer(cols=c(sigl_ite,sigl_ite_sig,sigldiff),names_to = "parameter", values_to = "value" ) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="-brewer.rd_bu",limits=sigma_limits),fill.legend = tm_legend(title = TeX("$\\sigma_l$"))) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 9, panel.labels = toplabel) 
  tmap_save(t,filename=paste0("../Documents/",folder_name,"/sigma_lower_",cond_site_name,".png"),width=8,height=6)
  
  t <- tmpsf %>% dplyr::select(sigu_ite,sigu_ite_sig,sigudiff) %>% pivot_longer(cols=c(sigu_ite,sigu_ite_sig,sigudiff),names_to = "parameter", values_to = "value" ) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="-brewer.rd_bu",limits=sigma_limits),fill.legend = tm_legend(title = TeX("$\\sigma_u$"))) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 9, panel.labels = toplabel) 
  tmap_save(t,filename=paste0("../Documents/",folder_name,"/sigma_upper_",cond_site_name,".png"),width=8,height=6)
  
  return(tmpsf)
}


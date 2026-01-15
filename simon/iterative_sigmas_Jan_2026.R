# change function to allow extra parameter

# start with the likelihood function
NLL_exp_phis <- function(phi,x = Z, d1j, mu1=as.numeric(unlist(mu_agg[,1])),deltal=2,deltau=2,phiBl=0,phiBu=0) {
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
  if (is.null(phiBl)==FALSE) {
    phi <- append(phi,phiBl,after=6)
  }
  if (is.null(phiBu)==FALSE) {
    phi <- append(phi,phiBu,after=7)
  }
  deltal <- phi[5]
  deltau <- phi[6]
  phiBl <- phi[7]
  phiBu <- phi[7]
  if(deltal<1 |deltau<1 | phi[1]<0 | phi[2] < 0 | phi[3]<0 | phi[4]<0 | phiBl<0 | phiBu<0 ){return(10e10)}
  
  sigu <- phiBu + phi[1]*(1-exp(-(phi[2]*dij.)))
  sigl <- phiBl + phi[3]*(1-exp(-(phi[4]*dij.)))
  C_AGG <-  (sigl/deltal*gamma(1/deltal) + sigu/deltau*gamma(1/deltau)  )^(-1)
  log_lik[x<mu1] <- log(C_AGG[x<mu1])-((mu1[x<mu1]-x[x<mu1])/sigl[x<mu1])^deltal 
  log_lik[x>=mu1] <- log(C_AGG[x>=mu1])-((x[x>=mu1]-mu1[x>=mu1])/sigu[x>=mu1])^deltau 
  return(-sum(log_lik))
}

par_est_ite <- function(dataLap=data_Laplace,v=q,given=cond_site,a,b,cond_site_dist, parest_site = result[[1]],Nite=10, show_ite=FALSE,deltal=NULL,deltau=NULL)  {
  d <- (ncol(dataLap)-1)
  N <- nrow(dataLap)
  res <- 1:d
  mu_agg <- sigl <- sigu <- data.frame(matrix(ncol=(Nite),nrow = d))
  
 phi1l. <- phi2l. <- phi1u. <- phi2u. <- c(1) 
 phi0l. <- phi0u. <- c(0)
  # calculate observed residuals
  Z <- observed_residuals(df=dataLap,given=given,v = v,a=a,b=b)
  
  if (is.null(deltal)) {
    residual_pars <- list(sigl = parest_site$sigl_ite_sigl,
                          sigu = parest_site$sigu_ite_sigu,
                          deltal = parest_site$deltal_ite[1],
                          deltau = parest_site$deltau_ite[1])
    deltal. <- 2
    deltau. <- 2
    
  } else {
    residual_pars <- list(sigl = parest_site$sigl_ite_sigl,
                          sigu = parest_site$sigu_ite_sigu,
                          "deltal" = deltal,
                          "deltau" = deltau)
    deltal. <- deltal
    deltau. <- deltau   
  }
  for (k in 1:Nite) {
    if (k >1) {
      residual_pars <- list(sigl=phi2*(1-exp(-phi3*cond_site_dist)),sigu=phi0*(1-exp(-phi1*cond_site_dist)),deltal=deltal,deltau=deltau)
    }
    # estimate mu
    pe <- par_est_mu(df=dataLap,v=v,given=given,res_margin_est=residual_pars,a=a,b=b)
     # update parameters
    mu_agg[,k] <- pe$mu
    
    
    # estimate sigu parameters phi0 and phi1
    phi_init <- c(phi0.[k],phi1.[k],phi2.[k],phi3.[k],deltal.[k],deltau.[k])
    opt <- optim(fn=NLL_exp_phis,x = Z,d1j=cond_site_dist,mu1=as.numeric(mu_agg[,k]),control=list(maxit=2000),par = phi_init,deltal=deltal,deltau=deltau,method = "Nelder-Mead")
    phi0 <- opt$par[1]
    phi1 <- opt$par[2]
    phi2 <- opt$par[3]
    phi3 <- opt$par[4] 
    if (!is.numeric(deltal) & !is.numeric(deltau)) {
      deltal <- opt$par[5]
      deltau <- opt$par[6] } 
    phi0. <- append(phi0.,phi0)
    phi1. <- append(phi1.,phi1)
    phi2. <- append(phi2.,phi2)
    phi3. <- append(phi3.,phi3)
    deltal. <- append(deltal.,deltal)
    deltau. <- append(deltau.,deltau)
    sigu[,k] <- phi0*(1-exp(-phi1*cond_site_dist))
    sigl[,k] <- phi2*(1-exp(-phi3*cond_site_dist))
  }
  par_sum <- data.frame("a" = as.numeric(a[,Nite]),"b" = as.numeric(b[,Nite]),"mu_agg" = as.numeric(mu_agg[,Nite]),"sigl" = as.numeric(sigl[,Nite]),"sigu" = as.numeric(sigu[,Nite]),"phi0" = phi0.[Nite], "phi1" = phi1.[Nite],"phi2" = phi2.[Nite], "phi3" = phi3.[Nite], "deltal" = deltal.[Nite], "deltau" = deltau.[Nite])
  if (show_ite == TRUE) {
    return(list(a,b,mu_agg,sigl,sigu,phi0.,phi1.,phi2.,phi3.,deltal.,deltau.,par_sum))
  } else {return(par_sum)}
}

par_est_mu <- function(df=sims,v=0.99,given=c(1),res_margin_est = res_margin_est,a,b) {
  lik <- mu_hat <- res_var <- c()
  names(df) <- paste0("Y",1:ncol(df))
  d <- ncol(df)
  for (j in given) {
    Y_given1extreme <- df %>% filter(df[,j]>quantile(df[,j],v))
    res <- c(1:d)[-j]
    init_par <- c()
    init_lik <- c()
    sigl <- res_margin_est$sigl
    sigu <- res_margin_est$sigu
    deltal <- res_margin_est$deltal
    deltau <- res_margin_est$deltau
    for (i in 2:d) {
      # optimise using the initial parameters
      Y1 <- Y_given1extreme[,j]
      Y2 <- Y_given1extreme[,res[i-1]]
      init_par <- c(0.8,0.2,0)
      opt <- optim(fn=NLL_AGG_onestep,x=data.frame(Y1,Y2),sigl_hat = sigl[i-1], sigu_hat = sigu[i-1], deltal_hat = deltal, deltau_hat = deltau,par=init_par,control=list(maxit=2000))
      a_hat <- append(a_hat,opt$par[1])
      b_hat <- append(b_hat,opt$par[2])
      mu_hat <- append(mu_hat,opt$par[3])
      lik <- append(lik,opt$value)
      res_var <- append(res_var,res[i-1])
    }
  }
  nas <- rep(NA,length(a_hat)) # NA values for parameters not used by a given method
  par_sum <- data.frame("lik"=lik, 
                        "a" = a_hat,"b" = b_hat,
                        "mu" = mu_hat,
                        "sigl" = sigl,"sigu" = sigu,
                        "deltal" = deltal,"deltau" = deltau,
                        "given" = rep(given,each=(d-1)),"res" = res_var)  
  
  return(par_sum)
}


abmu_par_est_ite <- function(site,Nite=10,sites = sites_index_diagonal,cond_site_names = site_name_diagonal,q=0.9,grid=xyUK20_sf,result,est_all_sf,deltal=NULL,deltau=NULL,folder_name=NULL) {
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

  if (is.null(folder_name)) {
    folder_name <- "abmu_iterative" 
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

result <- sapply(1:1,FUN = iter_sigmas_site,sites = sites_index_diagonal,cond_site_names = site_name_diagonal,par_est = est_all_diag1,folder_name = "Birmingham_Cromer_diagonal/sigmas",simplify = FALSE)



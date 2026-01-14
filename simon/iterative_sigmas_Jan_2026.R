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
  if(deltal<1 |deltau<1 ){return(10e10)}
  
  sigu <- phiBu + phi[1]*(1-exp(-(phi[2]*dij.)))
  sigl <- phiBl + phi[3]*(1-exp(-(phi[4]*dij.)))
  C_AGG <-  (sigl/deltal*gamma(1/deltal) + sigu/deltau*gamma(1/deltau)  )^(-1)
  log_lik[x<mu1] <- log(C_AGG[x<mu1])-((mu1[x<mu1]-x[x<mu1])/sigl[x<mu1])^deltal 
  log_lik[x>=mu1] <- log(C_AGG[x>=mu1])-((x[x>=mu1]-mu1[x>=mu1])/sigu[x>=mu1])^deltau 
  return(-sum(log_lik))
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

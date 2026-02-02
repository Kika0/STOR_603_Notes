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
result_previous <- result 
# identify diagonal sites
# find indeces of start and end sites
Birmingham <- c(-1.9032,52.4806)
Cromer <- c(1.28486,53.05349)
site_start <- find_site_index(site=Birmingham,grid_uk = xyUK20_sf)
site_end <- find_site_index(site=Cromer,grid_uk = xyUK20_sf)

# change functions to allow extra parameters phi0u (and optional also phi0l)

# start with the likelihood function
#' Title
#'
#' @param phi A vector of parameters
#' @param x A dataframe of observed residuals
#' @param d1j A vector of distance from the conditioning site
#' @param mu1 A vector of mu parameters
#' @param deltal Optional fixed value of deltal
#' @param deltau Optional fixed value of deltau
#' @param phi0l Optional fixed value of phi0l
#' @param phi0u Optional fixed value of phi0u
#'
#' @return
#' @export
#'
#' @examples
NLL_exp_phis <- function(phi,x, d1j, mu1=as.numeric(unlist(mu_agg[,1])),deltal=NULL,deltau=NULL,phi0l=NULL,phi0u=NULL) {
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
  # parameter value constraints
  if(deltal<1 |deltau<1 | phi[1]<0 | phi[2] < 0 | phi[3]<0 | phi[4]<0 | phi0l<0 | phi0u<0 ){return(10e10)}
 # if(deltal<1 |deltau<1 | phi[1]<0 | phi[2] < 0 | phi[3]<0 | phi[4]<0 | phi0l<0.1 | phi0u<0 ){return(10e10)}
  
  sigu <- phi0u + phi[1]*(1-exp(-(phi[2]*dij.)))
  sigl <- phi0l + phi[3]*(1-exp(-(phi[4]*dij.)))
  C_AGG <-  (sigl/deltal*gamma(1/deltal) + sigu/deltau*gamma(1/deltau)  )^(-1)
  log_lik[x<mu1] <- log(C_AGG[x<mu1])-((mu1[x<mu1]-x[x<mu1])/sigl[x<mu1])^deltal 
  log_lik[x>=mu1] <- log(C_AGG[x>=mu1])-((x[x>=mu1]-mu1[x>=mu1])/sigu[x>=mu1])^deltau 
  return(-sum(log_lik))
}

#' Title
#'
#' @param z A dataframe of observed residuals
#' @param v A quantile threshold
#' @param given A numeric index of conditioning variable
#' @param res_margin_est A list of parameters: sigl,sigu,optional:deltal,deltau
#'
#' @return A list of residual margin estimates plus likelihood, conditioning index and residual site
#' @export
#'
#' @examples
par_est_mu <- function(z,v,given,res_margin_est) {
  lik <- mu_hat <- res <- c()
  d <- ncol(z)+1
  res <- c(1:d)[-given]
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
  }
  par_sum <- data.frame("lik"=lik, 
                        "mu" = mu_hat,
                        "sigl" = sigl,"sigu" = sigu,
                        "deltal" = deltal,"deltau" = deltau,
                        "given" = rep(given,each=(d-1)),"res" = res)  
  
  return(par_sum)
}

#' Title
#'
#' @param z A dataframe of observed residuals
#' @param v A quantile threshold
#' @param given A numeric index of conditioning site
#' @param cond_site_dist A numeric vector of distance from conditioning site
#' @param parest_site A list of parameters: sigl,sigu,optional:deltal,deltau
#' @param Nite A numeric number of iterations
#' @param show_ite A logical, if TRUE, estimates are saved for each iteration
#' @param deltal A numeric fixed deltal value, if NULL deltal is also estimated
#' @param deltau A numeric fixed deltau value, if NULL deltau is also estimated
#' @param phi0u A numeric fixed phi0u value, if NULL phi0u is estimated
#' @param phi0l A numeric fixed phi0l value, if NULL phi0l is estimated
#'
#' @return A list of estimated, depends on show_ite, last item par_sum are final iteration estimates
#' @export
#'
#' @examples
par_est_ite <- function(z,v=q,given=cond_site,cond_site_dist, parest_site, Nite=10, show_ite=FALSE,deltal=NULL,deltau=NULL,phi0u=NULL,phi0l=0)  {
  d <- ncol(z)
  N <- nrow(z)
  res <- 1:d[-given]
  mu_agg <- sigl <- sigu <- data.frame(matrix(ncol=(Nite),nrow = d))
  
  phi1l. <- phi2l. <- phi1u. <- phi2u. <- c(1) 
  phi0l. <- phi0u. <- c(0.2)
  
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
  phi0l_init <- phi0l
  for (k in 1:Nite) {
    if (k >1) {
      residual_pars <- list(sigl=phi0l + phi1l *(1-exp(-phi2l*cond_site_dist)),sigu=phi0u + phi1u*(1-exp(-phi2u*cond_site_dist)),deltal=deltal,deltau=deltau)
    }
    # estimate mu
    pe <- par_est_mu(z=z,v=v,given=given,res_margin_est=residual_pars)
    # update mu parameters
    mu_agg[,k] <- pe$mu
    # estimate phi parameters
    if (is.null(phi0l_init)) {
    phi_init <- c(phi1u.[k],phi2u.[k],phi1l.[k],phi2l.[k],phi0l.[k],phi0u.[k])
    } else {
      phi_init <- c(phi1u.[k],phi2u.[k],phi1l.[k],phi2l.[k],phi0u.[k])
    }
      
    opt <- optim(fn=NLL_exp_phis,x = Z,d1j=cond_site_dist,mu1=as.numeric(mu_agg[,k]),control=list(maxit=2000),par = phi_init,deltal=deltal,deltau=deltau,phi0l=phi0l_init,method = "Nelder-Mead")
    print(opt)
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
      phi0u <- opt$par[6]
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
    return(list( "mu_agg" = mu_agg, "sigl" = sigl, "sigu" = sigu, "phi0u" = phi0u., "phi1u" = phi1u., "phi2u" = phi2u., "phi0l" = phi0l., "phi1l" = phi1l., "phi2l" = phi2l., "deltal" = deltal., "deltau" = deltau.,"par_sum" = par_sum))
  } else {return(par_sum)}
}


AGG_par_est_ite <- function(data_mod_Lap,site,v=0.9,Nite=10,sites = sites_index_diagonal,cond_site_names = site_name_diagonal,q=0.9,grid20km=xyUK20_sf,result,est_all_sf,deltal=NULL,deltau=NULL,phi0l=NULL,folder_name=NULL) {
  # filter conditioning site name from a vector of names
  if(is.null(cond_site_names)) {
    cond_site_name <- names(sites)[site]
    cond_site_names <- names(sites)
  } else {  cond_site_name <- cond_site_names[site] }
  # filter conditioning site index from a vector of indeces
  if (is.numeric(sites)) {cond_site <- sites[site]} else{
    cond_site_coord <- sites %>% dplyr::select(all_of(cond_site_name)) %>% pull()
    cond_site <- find_site_index(cond_site_coord,grid_uk = grid20km)    }
  print(cond_site)
  # calculate distance from the conditioning site
  dist_tmp <- as.numeric(unlist(st_distance(grid20km[cond_site,],grid20km)))
  # remove zero distance (conditioning site)
  dist_tmp <- dist_tmp[dist_tmp>0]
  # normalise distance using a common constant
  distnorm <- dist_tmp/1000000
  print(summary(distnorm))
  parest_site <- st_drop_geometry(result[[site]]) %>% dplyr::select(sigl_ite_sigl,sigu_ite_sigu,deltal_ite,deltau_ite) %>% na.omit()
  # calculate observed residuals
  aest <- discard(est_all_sf %>% filter(cond_site==cond_site) %>% pull(a),is.na)
  best <- discard(est_all_sf %>% filter(cond_site==cond_site) %>% pull(b),is.na)
  Z <- observed_residuals(df=data_mod_Lap,given=cond_site,v = v,a=aest,b=best)
print(names(Z))
    try7 <- par_est_ite(z=Z,given=cond_site,cond_site_dist=distnorm, parest_site = parest_site,Nite=Nite,show_ite=TRUE,deltal=deltal,deltau= deltau,phi0l=phi0l) 

  
 
  # separate parameter estimation and analysis
  p <- ggplot(try7$mu_agg %>% pivot_longer(everything(),names_to = "iteration",values_to = "par") %>% mutate(iteration=factor(iteration,levels=paste0("X",1:Nite)))) + geom_boxplot(aes(x=iteration,y=par)) +ylab(TeX("$\\mu_{AGG}$"))
  ggsave(p,file=paste0("../Documents/",folder_name,"/boxplot_mu_agg_",cond_site_name,".png"),height=5,width=10)
  p <- ggplot(try7$sigl %>% pivot_longer(everything(),names_to = "iteration",values_to = "par")%>% mutate(iteration=factor(iteration,levels=paste0("X",1:Nite)))) + geom_boxplot(aes(x=iteration,y=par)) +ylab(TeX("$\\sigma_l$"))
  ggsave(p,file=paste0("../Documents/",folder_name,"/boxplot_sigmal_",cond_site_name,".png"),height=5,width=10)
  p <- ggplot(try7$sigu %>% pivot_longer(everything(),names_to = "iteration",values_to = "par")%>% mutate(iteration=factor(iteration,levels=paste0("X",1:Nite)))) + geom_boxplot(aes(x=iteration,y=par)) +ylab(TeX("$\\sigma_u$"))
  ggsave(p,file=paste0("../Documents/",folder_name,"/boxplot_sigmau_",cond_site_name,".png"),height=5,width=10)
  
  tmp_phis <- data.frame("param_value"=c(try7$phi0l,try7$phi1l,try7$phi2l,try7$phi0u,try7$phi1u,try7$phi2u),"iteration"=paste0("X",0:Nite), "phi" = rep(c("phi0l","phi1l","phi2l","phi0u","phi1u","phi2u"),each=Nite+1))
  p <- ggplot(tmp_phis %>% mutate(iteration=factor(iteration,levels=paste0("X",0:Nite)))) + geom_point(aes(x=iteration,y=param_value,col=phi)) +ylab(TeX("$\\phi$"))
  ggsave(p,file=paste0("../Documents/",folder_name,"/iteration_phis_",cond_site_name,".png"),height=5,width=10)
  
  
  # look spatially to check
  # explore also spatial parameters
  est_ite <- try7$par_sum %>% add_row(.before=cond_site)
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
  
  misscol <- "aquamarine"
  #mu_limits <- c(-2.61,1)
  t <- tmpsf %>% dplyr::select(mu_agg_ite,mu_agg,mudiff) %>% pivot_longer(cols=c(mu_agg_ite,mu_agg,mudiff),names_to = "parameter", values_to = "value" ) %>% mutate(parameter=factor(parameter,levels=c("mu_agg_ite","mu_agg","mudiff"))) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="-brewer.rd_bu",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\mu_{AGG}$"))) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 12, panel.labels = toplabel,legend.reverse = TRUE) 
  
  tmap_save(t,filename=paste0("../Documents/",folder_name,"/mu_agg_",cond_site_name,".png"),width=8,height=6)
  
  #sigma_limits <- c(-1.32,2.3)
  t <- tmpsf %>% dplyr::select(sigl_ite,sigl,sigldiff) %>% pivot_longer(cols=c(sigl_ite,sigl,sigldiff),names_to = "parameter", values_to = "value" ) %>% mutate(parameter=factor(parameter,levels=c("sigl_ite","sigl","sigldiff"))) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="-brewer.rd_bu",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\sigma_l$"))) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 9, panel.labels = toplabel,legend.reverse = TRUE) 
  tmap_save(t,filename=paste0("../Documents/",folder_name,"/sigma_lower_",cond_site_name,".png"),width=8,height=6)
  
  t <- tmpsf %>% dplyr::select(sigu_ite,sigu,sigudiff) %>% pivot_longer(cols=c(sigu_ite,sigu,sigudiff),names_to = "parameter", values_to = "value" ) %>% mutate(parameter=factor(parameter,levels=c("sigu_ite","sigu","sigudiff"))) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="-brewer.rd_bu",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\sigma_u$"))) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 9, panel.labels = toplabel,legend.reverse = TRUE) 
  tmap_save(t,filename=paste0("../Documents/",folder_name,"/sigma_upper_",cond_site_name,".png"),width=8,height=6)
  return(try7)
}
folder_name <- "Birmingham_Cromer_diagonal/new_iterative_sigmas_mu"
result_new <- sapply(1:12,FUN = function(site_order){AGG_par_est_ite(site=site_order,data_mod_Lap = data_mod_Lap,sites = sites_index_diagonal,cond_site_names = site_name_diagonal,est_all_sf = est_all_sf,Nite=10,result=result_previous,deltal=deltal,deltau=deltau,folder_name = folder_name,phi0l=0)},simplify = FALSE)

summary(result_new[[1]])

# examine outliers
plot_AGG_diagnostics_method <- function(site=1,method_name = "Original_method",result=result_new) {
  # try do a PP plot
  # calculate observed residuals
  aest <- discard(est_all_sf %>% filter(cond_site==site_name_diagonal[site]) %>% pull(a),is.na)
  best <- discard(est_all_sf %>% filter(cond_site==site_name_diagonal[site]) %>% pull(b),is.na)
  Z <- observed_residuals(df=data_mod_Lap,given=sites_index_diagonal[site],v = q,a=aest,b=best)
  U <- c(1:(dim(data_mod_Lap)[1]*(1-q)))/(dim(data_mod_Lap)[1]*(1-q)+1)
  # pick a site
  sites_residual <- c(1:ncol(data_mod_Lap))[-sites_index_diagonal[site]]
  # get estimates new
  est_new <- result[[site]][[12]] %>% add_row(.before=sites_index_diagonal[site])
  # set up dataframe
  tmp <- data.frame(x=numeric(),y=numeric(),"method"=character(),"res_site"=numeric())
  for (i in 1:(ncol(Z))) {
    site_res <- sites_residual[i]
    AGGPars <- as.numeric(unlist(est_new[site_res,c(1:3,8,9)]))
    # get estimates old
    AGGParsOld <- as.numeric(unlist(st_drop_geometry(est_all_sf%>% filter(cond_site==site_name_diagonal[site]))[site_res,c(13:17)]))
    # construct PP plot
    x <- as.numeric(unlist(Z[,i]))
    Um <- pAGG(x=x,theta = AGGPars)
    Umo <- pAGG(x=x,theta = AGGParsOld)
    tmp_append_old <- data.frame(y=Umo,method="Original_method","res_site"=site_res) %>% arrange(y) %>% mutate("x"=U)
    tmp_append_new <- data.frame(x=U,y=Um,method="New_iterative_approach","res_site"=site_res) %>% arrange(y) %>% mutate("x"=U)
    tmp <- rbind(tmp,tmp_append_old,tmp_append_new) 
  }
  # PP plot of the new iterative approach and the original method
  p <- ggplot(tmp) + geom_line(aes(x=x,y=y,group=res_site),linewidth=0.01,alpha=0.9)+ geom_abline(slope=1,linetype="dashed") + facet_wrap("method") + xlab("Empirical") + ylab("Model")
  ggsave(p,filename=paste0("../Documents/",folder_name,"/PP_",site_name_diagonal[site],".pdf"),width=10,height=5)
  
  tmp <- tmp %>% dplyr::filter(method==method_name) %>% mutate(diag_diff=y-x)
  res_site_over <- tmp[tmp$diag_diff==max(tmp$diag_diff),]$res_site
  res_site_under <-  tmp[tmp$diag_diff==min(tmp$diag_diff),]$res_site
  
  AGG_underestimate <- append((tmp %>% group_by(res_site) %>% summarise(n=max(diag_diff,na.rm=TRUE)) %>% pull(n)),NA,after=sites_index_diagonal[site]-1)
  AGG_overestimate <- append((tmp %>% group_by(res_site) %>% summarise(n=min(diag_diff,na.rm=TRUE)) %>% pull(n)),NA,after=sites_index_diagonal[site]-1)
  # map worst distances from the diagonal
  tmpsf <- xyUK20_sf %>% mutate(AGG_underestimate,AGG_overestimate) %>% pivot_longer(cols=c(AGG_underestimate,AGG_overestimate),names_to="under_over",values_to = "diag_diff")
  title_map <- ""
  misscol <- "aquamarine"
    legend_text_size <- 0.7
    point_size <- 0.6
    legend_title_size <- 1.2
    lims_under <- c(-0.2,0)
    lims_over <- c(0,0.2)
    nrow_facet <- 1
    t1 <- tm_shape(tmpsf %>% filter(under_over=="AGG_overestimate")) + tm_dots(fill="diag_diff",fill.scale = tm_scale_continuous(limits=lims_under,values="-Blues",value.na=misscol,label.na = "Conditioning\n site"),size=point_size)  + tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE) + tm_title(text="AGG overestimate") 
    t2 <- tm_shape(tmpsf%>% filter(under_over=="AGG_underestimate")) + tm_dots(fill="diag_diff",fill.scale = tm_scale_continuous(limits=lims_over,values="Blues",value.na=misscol,label.na = "Conditioning\n site"),size=point_size)  + tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE) + tm_title(text="AGG underestimate") 
    t <- tmap_arrange(t1,t2,ncol=2)
    tmap_save(t,filename=paste0("../Documents/",folder_name,"/PP_all_map_examine_",method_name,"_",site_name_diagonal[site],".png"))
    
    # examine on a map only the worst case
    over_under <- rep(NA,nrow(xyUK20_sf))
    over_under[sites_index_diagonal[site]] <- "Conditioning_site"
    over_under[res_site_over] <- "Residual_site_underestimate"
    over_under[res_site_under] <- "Residual_site_overestimate"
    tmp_over_under <- xyUK20_sf %>% mutate("over_under"=over_under)
    t <- tm_shape(tmp_over_under) + tm_dots(fill="over_under",fill.scale = tm_scale_categorical(values=c("Residual_site_underestimate"="#C11432","Residual_site_overestimate" = "#009ADA", "Conditioning_site" = "aquamarine"))) 
    tmap_save(t,filename=paste0("../Documents/",folder_name,"/PP_outliers_map_examine_",method_name,"_",site_name_diagonal[site],".png"))
    # examine over
    if (method_name=="Original_method") {
      AGGPars <- as.numeric(unlist(st_drop_geometry(est_all_sf)[res_site_over,c(13:17)]))
      Zover <- as.numeric(unlist(Z[,which(res_site_over==sites)]))
    }
    if (method_name=="New_iterative_approach") {
      AGGPars <- as.numeric(unlist(est_new[res_site_over,c(1:3,8,9)]))
      Zover <- as.numeric(unlist(Z[,which(res_site_over==sites)]))
    }
    # plot density and kernel smooth
    tmpd <- data.frame(y=AGG_density(theta = AGGPars,x=seq(-7.5,5,0.01)),x=seq(-7.5,5,0.01))
    p1 <- ggplot(data.frame(Zover)) + geom_density(mapping = aes(x=Zover)) + geom_line(data=tmpd,mapping = aes(x=x,y=y),col="#C11432") + xlab(TeX(paste0("$Z_{",res_site_over,"}$"))) + ylab("Residual density function") + ggtitle("Underestimation")
    # examine under
    if (method_name=="Original_method") {
      AGGPars <- as.numeric(unlist(st_drop_geometry(est_all_sf)[res_site_under,c(13:17)]))
      Zunder <- as.numeric(unlist(Z[,which(res_site_under==sites)]))
    }
    if (method_name=="New_iterative_approach") {
      AGGPars <- as.numeric(unlist(est_new[res_site_under,c(1:3,8,9)]))
      Zunder <- as.numeric(unlist(Z[,which(res_site_under==sites)]))
    }
    
    # plot density and kernel smooth
    tmpd <- data.frame(y=AGG_density(theta = AGGPars,x=seq(-7.5,5,0.01)),x=seq(-7.5,5,0.01))
    p2 <- ggplot(data.frame(Zunder)) + geom_density(mapping = aes(x=Zunder)) + geom_line(data=tmpd,mapping = aes(x=x,y=y),col="#009ADA") + xlab(TeX(paste0("$Z_{",res_site_under,"}$"))) + ylab("Residual density function") + ggtitle("Overestimation")
    p <- grid.arrange(p1,p2,ncol=1)
    ggsave(p,filename=paste0("../Documents/",folder_name,"/AGG_outliers_",method_name,".pdf"),width=10,height=10)
    # compare worst fit with original approach
    # plot density and kernel smooth
    AGGParsorg <- as.numeric(unlist(st_drop_geometry(est_all_sf)[res_site_under,c(13:17)]))
    tmpdorg <- data.frame(y=AGG_density(theta = AGGParsorg,x=seq(-7.5,5,0.01)),x=seq(-7.5,5,0.01))
    tmpd <- data.frame(y=AGG_density(theta = AGGPars,x=seq(-7.5,5,0.01)),x=seq(-7.5,5,0.01))
    p3 <- ggplot(data.frame(Zunder)) + geom_density(mapping = aes(x=Zunder)) + geom_line(data=tmpd,mapping = aes(x=x,y=y),col="#FDD10A") + geom_line(data=tmpdorg,mapping = aes(x=x,y=y),col="#66A64F")+ xlab(TeX(paste0("$Z_{",res_site_under,"}$"))) + ylab("Residual density function") + ggtitle("Overestimation")
    ggsave(p3,filename=paste0("../Documents/",folder_name,"/AGG_outliers_compare",method_name,"_",site_name_diagonal[site],".pdf"),width=10,height=10)
    
    
      
    # distance from conditioning side vs worst distance of a given residual site
    dist_cond_site <- as.numeric(unlist(st_distance(xyUK20_sf[sites_index_diagonal[site],],xyUK20_sf) %>%
                                          units::set_units(km)))
    p1 <- ggplot(data.frame(AGG_underestimate,dist_cond_site)) + geom_point(aes(x=dist_cond_site,y=AGG_underestimate))
    p2 <- ggplot(data.frame(AGG_overestimate,dist_cond_site)) + geom_point(aes(x=dist_cond_site,y=AGG_overestimate))
    p <- grid.arrange(p1,p2,ncol=1)
    ggsave(p,filename=paste0("../Documents/",folder_name,"/AGG_distance_",method_name,"_",site_name_diagonal[site],".pdf"),width=5,height=5)
}

# run diagnostics for both methods
sapply(1:12,FUN=plot_AGG_diagnostics_method,method_name="Original_method",result=result_new)
sapply(1:12,FUN=plot_AGG_diagnostics_method,method_name="New_iterative_approach",result=result_new)

# repeat for phi0l also being estimated
folder_name <- "Birmingham_Cromer_diagonal/new_iterative_sigmas_mu_phi0l_phiul"
result_new <- sapply(1:1,FUN = AGG_par_est_ite,data_mod_Lap = data_mod_Lap,sites = sites_index_diagonal,cond_site_names = site_name_diagonal,phi0l=NULL,est_all_sf = est_all_sf,Nite=10,result=result_previous,folder_name = folder_name,simplify = FALSE)
sapply(1:12,FUN=plot_AGG_diagnostics_method,method_name="Original_method",result=result_new)
sapply(1:12,FUN=plot_AGG_diagnostics_method,method_name="New_iterative_approach",result=result_new)

# run diagnostics for both methods (Birmingham only)
#plot_AGG_diagnostics_method(site=1,method_name="Original_method",result=result_new_phi0l)
#plot_AGG_diagnostics_method(site=1,method_name="New_iterative_approach",)

# plot outliers 260 and 321 against Birmingham
p1 <- ggplot(data_mod_Lap %>% filter(Y192>quantile(Y192,q))) + geom_point(aes(x=Y192,y=Y260),size=0.5)
p2 <- ggplot(data_mod_Lap %>% filter(Y192>quantile(Y192,q))) + geom_point(aes(x=Y192,y=Y321),size=0.5)
p <- grid.arrange(p1,p2,ncol=2)
ggsave(p,filename=paste0("../Documents/",folder_name,"/AGG_Birmingham_outliers_Y.pdf"),width=10,height=5)

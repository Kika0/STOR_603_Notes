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
spatial_par_est <- function(data_Lap,cond_sites,dayshift=c(0),Ndays_season=90,v=0.9,ab_method="sequential2",res_margin="AGG",grid_uk=xyUK20_sf,title="") {
  est_all <- data.frame("lik" = numeric(), "lika"= numeric(),"likb"= numeric(),"lik2"= numeric(),
                            "a" = numeric(), "b" = numeric(),
                            "mu" = numeric(),"mu_agg" = numeric(),
                            "sig" = numeric(),"sig_agg" = numeric(),"sigl" = numeric(),"sigu" = numeric(),
                            "delta" = numeric(),"deltal" = numeric(), "deltau" = numeric(),
                            "given" = numeric(), "res" = numeric(), "cond_site" = character(), "tau" = numeric())
  for (i in 1:ncol(cond_sites)) {
    cond_site <- find_site_index(as.numeric(cond_sites[,i]),grid_uk = grid_uk)
    for (j in 1:length(dayshift)) {
      sims_tau <- shift_time(sims=data_Lap,cond_site=cond_site,tau=dayshift[j],Ndays_season = Ndays_season)
      pe <- par_est(df=sims_tau,v=v,given=cond_site,margin = "Normal", method=ab_method,keef_constraints = c(1,2))
      # calculate observed residuals
      obsr <- observed_residuals(df=sims_tau,given=cond_site,v = v,a=pe$a,b=pe$b)
      # estimate residual margin parameters
      pe_res <- res_margin_par_est(obs_res = obsr,method="AGG")
      est_all <- rbind(est_all,cbind(pe[,-c(8,10:12,14:15)],pe_res) %>% add_row(.before=cond_site) %>%  mutate(cond_site=names(cond_sites[i]),tau=as.character(dayshift[j])))
       }
  }
  est_all <- est_all %>% mutate(tau=factor(as.character(tau),levels = as.character(dayshift))) %>% mutate(cond_site=factor(cond_site))
  est_all_sf <- est_join_spatial(tmp_est=est_all,grid_uk=grid_uk)
  save(est_all_sf,file=paste0("data_processed/N",nrow(data_Lap),"_",ab_method,"_",res_margin,"_",title,".RData"))
}

res_margin_par_est_ite <- function(data=Z,given=cond_site,N=10, show_ite=FALSE,mu_init=NULL,sigl_init=NULL,sigu_init=NULL,deltal_init=NULL,deltau_init=NULL,method="onephi",SN=NULL, b_inc=FALSE)  {
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
    opt <- optim(fn=NLL_AGG_deltas,df = data,mu1=as.numeric(mu_agg[,i+1]),sigl1=as.numeric(sigl[,i+1]),sigu1=as.numeric(sigu[,i+1]),control=list(maxit=2000),par = c(as.numeric(deltal[1,i]),as.numeric(deltau[1,i])),method = "Nelder-Mead")
    deltal[,i+1] <- opt$par[1]
    deltau[,i+1] <- opt$par[2]
  }
  par_sum <- data.frame("mu_agg" = as.numeric(mu_agg[,N+1]),"sigl" = as.numeric(sigl[,N+1]),"sigu" = as.numeric(sigu[,N+1]),"deltal" = as.numeric(deltal[,N+1]), "deltau" = as.numeric(deltau[,N+1]))
  if (show_ite == TRUE) {
    return(list(mu_agg,sigl,sigu,deltal,deltau,par_sum))
  } else {return(par_sum)}
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


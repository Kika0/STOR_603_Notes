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
  log_lik[z<mu1] <- log(C_AGG[z<mu1])-((mu1[z<mu1]-z[z<mu1])/sigl[z<mu1])^deltal 
  log_lik[z>=mu1] <- log(C_AGG[z>=mu1])-((z[z>=mu1]-mu1[z>=mu1])/sigu[z>=mu1])^deltau 
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
par_est_ite <- function(z,v,given,cond_site_dist, parest_site, Nite=10, show_ite=FALSE,deltal=NULL,deltau=NULL,phi0u=NULL,phi0l=0)  {
  d <- ncol(z)
  N <- nrow(z)
  res <- 1:d[-given]
  mu_agg <- sigl <- sigu <- data.frame(matrix(ncol=(Nite),nrow = d))
  
  phi1l. <- phi2l. <- phi1u. <- phi2u. <- c(1) 
  phi0l. <- phi0u. <- c(0.2)
  
  if (is.null(deltal)) {
    residual_pars <- list(sigl = parest_site$sigl,
                          sigu = parest_site$sigu,
                          deltal = parest_site$deltal[1],
                          deltau = parest_site$deltau[1])
    deltal. <- 2
    deltau. <- 2
    
  } else {
    residual_pars <- list("sigl" = parest_site$sigl,
                          "sigu" = parest_site$sigu,
                          "deltal" = deltal,
                          "deltau" = deltau)
    deltal. <- deltal
    deltau. <- deltau   
  }
  phi0l_init <- phi0l
  phi0u_init <- phi0u
  for (k in 1:Nite) {
    if (k >1) {
      residual_pars <- list(sigl=phi0l + phi1l *(1-exp(-phi2l*cond_site_dist)),sigu=phi0u + phi1u*(1-exp(-phi2u*cond_site_dist)),deltal=deltal,deltau=deltau)
    }
    # estimate mu
    pe <- par_est_mu(z=z,v=v,given=given,res_margin_est=residual_pars)
    print(summary(pe))
    print(glimpse(pe))
    # update mu parameters
    mu_agg[,k] <- pe$mu
    # estimate phi parameters
    if (is.null(phi0l_init)) {
      phi_init <- c(phi1u.[k],phi2u.[k],phi1l.[k],phi2l.[k],phi0l.[k],phi0u.[k])
    } else {
      phi_init <- c(phi1u.[k],phi2u.[k],phi1l.[k],phi2l.[k],phi0u.[k])
    }
    
    opt <- optim(fn=NLL_exp_phis,x = z,d1j=cond_site_dist,mu1=as.numeric(mu_agg[,k]),control=list(maxit=2000),par = phi_init,deltal=deltal,deltau=deltau,phi0l=phi0l_init,method = "Nelder-Mead")
    print(opt)
    phi1u <- opt$par[1]
    phi2u <- opt$par[2]
    phi1l <- opt$par[3]
    phi2l <- opt$par[4] 
    if (!is.numeric(deltal) & !is.numeric(deltau)) {
      deltal <- opt$par[5]
      deltau <- opt$par[6] } 
    if (!is.numeric(phi0u_init) & is.numeric(phi0l_init)) {
      phi0u <- opt$par[5]
    }
    if (!is.numeric(phi0u_init) & !is.numeric(phi0l_init)) {
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

# Step 2. helper functions
to_opt <- function(data_Lapv,theta,i,cond_index=London_index,pe_i) {
  # a <- theta[1]
  # b <- theta[2]
  # if (a<0 | a>1 | b<0 | b>1) {return(10^6)}
  y <-  NLL_AGG_onestep(x=data_Lapv,theta=theta,sigl_hat = pe_i[2], sigu_hat = pe_i[3], deltal_hat = pe_i[4], deltau_hat = pe_i[5])
  return(y)
}
NLL_AGG_wrapper <- function(data_Lapv,i,pe_res,cond_index) {
  pe_i <- as.numeric(unlist(pe_res[i,]))
  if (is.na(pe_i[1])) {return(NA)}
  y <- optim(par=c(0.8,0.3,pe_i[1]),fn= to_opt,data_Lapv=data_Lapv,i=i,cond_index=cond_index,pe_i=pe_i)
  return(y$par) 
}

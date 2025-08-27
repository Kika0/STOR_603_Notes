sigmau_par_est_ite <- function(data=Z,given=cond_site,cond_site_dist, N=10, show_ite=FALSE,mu_init=NULL,sigl_init=NULL,sigu_init=NULL,deltal,deltau)  {
  d <- ncol(data)
  nv <- nrow(data)
  res <- 1:d
  mu_agg <- sigl <- sigu <- data.frame(matrix(ncol=(N+1),nrow = d))
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
phi0. <- phi1. <- c(1,1) 
  for (i in 1:N) {
    # estimate sigu parameters phi0 and phi1
    phi_init <- c(phi0.[i],phi1.[i])
    opt <- optim(fn=NLL_exp_sigmau,x = data,d1j=cond_site_distance,mu1=as.numeric(mu_agg[,i]),sigl1=as.numeric(sigl[,i]),deltal=deltal,deltau=deltau,control=list(maxit=2000),par = phi_init,method = "Nelder-Mead")
    phi0 <- opt$par[1]
    phi1 <- opt$par[2]
    phi0. <- append(phi0.,phi0)
    phi1. <- append(phi1.,phi1)
    a[,1] <- exp(-(phi1*as.numeric(SN)+phi0*as.numeric(!SN))*d1j)
    
    sigmau[,i+1] <- phi0*(1-exp(-phi1*cond_site_distance))
    # estimate sigl and mu_agg for each site separately
    for (j in 1:d) {
      Z2 <- as.numeric(unlist(data[,j]))
      opt <- optim(fn=NLL_AGG,x=Z2,deltal_hat=as.numeric(deltal[1,i]),deltau_hat=as.numeric(deltau[1,i]),par=c(mean(Z2),sd(Z2),sd(Z2)),control=list(maxit=2000),method = "Nelder-Mead")
      mu_agg[j,i+1] <- opt$par[1]
      sigl[j,i+1] <- opt$par[2]
      sigu[j,i+1] <- opt$par[3]
    }

  }
  par_sum <- data.frame("mu_agg" = as.numeric(mu_agg[,N+1]),"sigl" = as.numeric(sigl[,N+1]),"sigu" = as.numeric(sigu[,N+1]),"deltal" = as.numeric(deltal[,N+1]), "deltau" = as.numeric(deltau[,N+1]))
  if (show_ite == TRUE) {
    return(list(mu_agg,sigl,sigu,deltal,deltau,par_sum))
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
# calculate distance from the conditioning site

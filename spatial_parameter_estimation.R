#' Estimate conditional model parameters for spatio-temporal data
#'
#' `spatial_par_est()` estimates $\alpha$ and $\beta$ dependence parameters,
#'  residual margin parameters and links back to spatial locations
#'
#' @param data_Lap A dataframe of observed data on Laplace scale Y
#' @param cond_sites A vector of conditioning sites (a subset of 12 possible sites)
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
spatial_par_est <- function(data_Lap,cond_sites,dayshift=c(0),v=0.9,ab_method="sequential2",res_margin="AGG",grid_uk=xyUK20_sf) {
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
  
  est_all <- data.frame("lik" = numeric(), "lika"= numeric(),"likb"= numeric(),"lik2"= numeric(),
                            "a" = numeric(), "b" = numeric(),
                            "mu" = numeric(),"mu_agg" = numeric(),
                            "sig" = numeric(),"sig_agg" = numeric(),"sigl" = numeric(),"sigu" = numeric(),
                            "delta" = numeric(),"deltal" = numeric(), "deltau" = numeric(),
                            "given" = numeric(), "res" = numeric(), "cond_site" = character(), "tau" = numeric())
  for (i in 1:ncol(cond_sites)) {
    cond_site <- find_site_index(as.numeric(cond_sites[,i]),grid_uk = grid_uk)
    for (j in 1:length(dayshift)) {
      sims_tau <- shift_time(sims=data_Lap,cond_site=cond_site,tau=dayshift[j],Ndays_season = 90)
      pe <- par_est(df=sims_tau,v=v,given=cond_site,margin = "Normal", method=ab_method,keef_constraints = c(1,2))
      # calculate observed residuals
      obsr <- observed_residuals(df=data_Lap,given=cond_site,v = v,a=pe$a,b=pe$b)
      # estimate residual margin parameters
      pe_res <- res_margin_par_est(obs_res = obsr,method="AGG")
      est_all <- rbind(est_all,cbind(pe_res,pe) %>% add_row(.before=cond_site) %>%  mutate(cond_site=names(df_sites[i]),tau=as.character(dayshift[j])))
       }
  }
  est_all <- est_all %>% mutate(tau=factor(as.character(tau),levels = as.character(dayshift)))
  est_all_sf <- est_join_spatial(tmp_est=est_all,grid_uk=grid_uk)
  save(est_all_sf,file=paste0("data_processed/",data_Lap,"_",ab_method,"_",res_margin,".RData"))
}

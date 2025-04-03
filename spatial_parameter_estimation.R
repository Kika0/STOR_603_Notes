#' Estimate conditional model parameters for spatio-temporal data
#'
#' `spatial_par_est()` estimates $\alpha$ and $\beta$ dependence parameters,
#'  residual margin parameters and links back to spatial locations
#'
#' @param data_Lap A dataframe of observed data on Laplace scale Y
#' @param cond_sites A vector of conditioning sites (a subset of 12 possible sites)
#' @param taus A numeric vector of temporal lags (days)
#' @param v A quantile
#' @param ab_method A method to estimate dependence parameters $\alpha$ and $\beta$
#' @param res_margin A parametric method to estimate residual margins
#'
#' @return `spatial_par_est()` returns an sf object of spatial points with parameter estimates
#' @export
#'
#' @examples
spatial_par_est <- function(data_Lap,cond_sites, taus=c(0),v=0.9,ab_method="sequential2",res_margin="AGG") {
  est_all <- data.frame("lik" = numeric(), "lika"= numeric(),"likb"= numeric(),"lik2"= numeric(),
                            "a" = numeric(), "b" = numeric(),
                            "mu" = numeric(),"mu_agg" = numeric(),
                            "sig" = numeric(),"sig_agg" = numeric(),"sigl" = numeric(),"sigu" = numeric(),
                            "delta" = numeric(),"deltal" = numeric(), "deltau" = numeric(),
                            "given" = numeric(), "res" = numeric(), "cond_site" = character(), "tau" = numeric())
  for (i in 1:ncol(df_sites)) {
    cond_site <- find_site_index(as.numeric(df_sites[,i]),grid_uk = xyUK20_sf)
    for (j in 1:length(dayshift)) {
      sims_tau <- shift_time(sims=sims,cond_site=cond_site,tau=dayshift[j],Ndays_season = 92)
      est_all <- rbind(est_all,par_est(df=sims_tau,v=v,given=cond_site,margin = "AGG", method="sequentialGG",keef_constraints = c(1,2)) %>% add_row(.before=cond_site) %>%  mutate(cond_site=names(df_sites[i]),tau=as.character(dayshift[j])))
    }
  }
  
  save(est_all,file=paste0("data_processed/",data_Lap,"_",ab_method,"_",res_margin,".RData"))
}
library(tmap)
library(sf)
library(tidyverse)
library(latex2exp)

# identify diagonal sites
# find indeces of start and end sites
Birmingham <- c(-1.9032,52.4806)
Cromer <- c(1.28486,53.05349)
site_start <- find_site_index(site=Birmingham,grid_uk = xyUK20_sf)
site_end <- find_site_index(site=Cromer,grid_uk = xyUK20_sf)

tmap_mode("view")
uk_diag <- xyUK20_sf %>% mutate(siteID=as.numeric(1:nrow(xyUK20_sf)))
tm_shape(uk_diag) + tm_dots(fill="temp",size=1)
sites_index_diagonal <- c(192,193,194,195,196,197,218,219,220,221,242,263) # first is Birmingham and last is Cromer
site_name_diagonal <- c("Birmingham", paste0("diagonal",1:(length(sites_index_diagonal)-2)),"Cromer")
# plot these points on a map
uk_diag <- uk_diag %>% mutate(sites_diagonal=factor(case_match(siteID,c(site_start) ~ "Birmingham",c(site_end)~"Cromer",sites_index_diagonal[2:(length(site_name_diagonal)-1)]~"diagonal_sites")))
tmap_mode("plot")
tm_shape(uk_diag) + tm_dots("sites_diagonal",size=0.1,fill.scale = tm_scale_categorical(values=c("Birmingham"="#C11432","Cromer" = "#009ADA", "diagonal_sites" = "#FDD10A")))

# estimate parameters along
spatial_par_est(data_Lap = data_mod_Lap,cond_sites = sites_index_diagonal,cond_site_names = site_name_diagonal,dayshift = 0,v=q,Ndays_season = 90,title = paste0("diagonal_sites_Birmingham_Cromer",q*100))

# repeat iterative procedures
# load estimates to calculate observed residuals
load("data_processed/N9000_sequential2_AGG_diagonal_sites_Birmingham_Cromer90.RData")
summary(est_all_sf)

est_all_diag1 <- est_all_sf %>% mutate(cond_site = factor(as.character(cond_site),levels=as.character(site_name_diagonal)))
tm <- map_param(tmp_est=est_all_diag1, method = "AGG", facet_var = c("cond_site"), title_map = "Sites from Birmingham to Cromer")
tm[[4]]
condmodel_params <- c("a","b","mu","sig","muagg","sigl","sigu","sigdiff","deltal","deltau","deltadiff")

save_map_i <- function(i,tm_list=tm,doc_folder = "Birmingham_Cromer_diagonal",w = 8,h = 8,par_vect = condmodel_params) {
 tmap_save(tm_list[[i]],filename = paste0("../Documents/",doc_folder,"/allsites_",par_vect[i],"_original.png"),width = w,height = h) 
}
sapply(1:length(condmodel_params),FUN = save_map_i)

# move to delta estimates
result <- sapply(1:length(sites_index_diagonal),FUN = iter_delta_site,sites = sites_index_diagonal,cond_site_names = site_name_diagonal,par_est = est_all_diag1,folder_name = "Birmingham_Cromer_diagonal/delta",simplify = FALSE)

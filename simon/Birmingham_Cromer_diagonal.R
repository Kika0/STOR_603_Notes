library(tmap)
library(sf)
library(dplyr)


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
#spatial_par_est(data_Lap = data_mod_Lap,cond_sites = sites_index_diagonal,cond_site_names = site_name_diagonal,dayshift = 0,v=q,Ndays_season = 90,title = paste0("diagonal_sites_Bormingham_Cromer",q*100))

# repeat iterative procedures
# load estimates to calculate observed residuals




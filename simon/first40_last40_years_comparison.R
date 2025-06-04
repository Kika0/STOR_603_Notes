library(tmap)
library(tidyverse)
library(latex2exp)
theme_set(theme_bw())
theme_replace(
  panel.spacing = unit(2, "lines"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA) )
file.sources = list.files(pattern="*helpers.R")
sapply(file.sources,source,.GlobalEnv)
# source("spatial_parameter_estimation.R") # for spatial_par_est function
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
#spatial_par_est(data_Lap = data_mod_Lap,cond_sites = df_sites,dayshift = 0,Ndays_season = 90,title = "all12sites")
#spatial_par_est(data_Lap = data_mod_Lap[1:3600,],cond_sites = df_sites,dayshift = 0,Ndays_season = 90,title = "first40_12sites")
#spatial_par_est(data_Lap = data_mod_Lap[5401:9000,],cond_sites = df_sites,dayshift = 0,Ndays_season = 90,title = "last40_12sites")

# load all three parameter estimates sf objects ----
load("data_processed/N9000_sequential2_AGG_all12sites.RData")
est_all <- est_all_sf

load("data_processed/N3600_sequential2_AGG_first40_12sites.RData")
est_allfirst40 <- est_all_sf

load("data_processed/N3600_sequential2_AGG_last40_12sites.RData")
est_alllast40 <- est_all_sf

# combine the first and the last 40 years
est_comb <- bind_rows(est_allfirst40 %>% mutate("Year_range"=c("1981_2020")),est_alllast40 %>% mutate("Year_range"=c("2041_2080")))

p1 <- ggplot(est_comb) + geom_density(aes(x=a,fill=Year_range),alpha=0.5)  +  scale_fill_manual(name="Year_range",labels=c(TeX("$1980-2020$"),
                                                                                                             TeX("$2041-2080$")),
                                                                                       values=c("black","#C11432")) + facet_wrap(~cond_site)


p2 <- ggplot(est_comb) + geom_density(aes(x=b,fill=Year_range),alpha=0.5)  +  scale_fill_manual(name="Years",labels=c(TeX("$1980-2020$"),
                                                                                                             TeX("$2041-2080$")),
                                                                                       values=c("black","#C11432"))+ facet_wrap(~cond_site)



ggsave(p1,width=10,height=5,filename = "../Documents/acomb_12sites.png")
ggsave(p2,width=10,height=5,filename = "../Documents/bcomb_12sites.png")

# plot estimates spatially for each site ------
# bind data together
est_comball <- bind_rows(est_all %>% mutate("Year_range"=c("1981_2080")),est_allfirst40 %>% mutate("Year_range"=c("1981_2020")),est_alllast40 %>% mutate("Year_range"=c("2041_2080")))
est_comball <- est_comball %>% mutate(Year_range=factor(Year_range,levels=c("1981_2080","1981_2020","2041_2080")))
for (i in 1:12) {
  pa <- tm_shape(est_comball %>% filter(cond_site==names(df_sites)[i]))  + tm_dots(col="a",palette="viridis",n=8,size=0.3,colorNA="aquamarine",title=TeX("$\\alpha$"),textNA = "Conditioning site") + tm_facets(by= c("cond_site","Year_range"))+  tm_layout(legend.outside.size=0.3,asp=0.5,legend.text.size = 1,legend.title.size=1.5)
    tmap_save(pa,filename=paste0("../Documents/acomp_",names(df_sites[i]),".png"),width=10,height=6)
    pb <- tm_shape(est_comball %>% filter(cond_site==names(df_sites)[i]))  + tm_dots(col="b",palette="viridis",n=8,size=0.3,colorNA="aquamarine",title=TeX("$\\beta$"),textNA = "Conditioning site") + tm_facets(by= c("cond_site","Year_range"))+  tm_layout(legend.outside.size=0.3,asp=0.5,legend.text.size = 1,legend.title.size=1.5)
    tmap_save(pb,filename=paste0("../Documents/bcomp_",names(df_sites[i]),".png"),width=10,height=6)
    
    }



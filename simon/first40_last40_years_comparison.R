library(tmap) # spatial map plots
library(sf) # for handling spatial sf objects
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
#spatial_par_est saves parameter estimates as est_all_sf sf object in data_processed folder
q <- 0.98 # quantile threshold
# spatial_par_est(data_Lap = data_mod_Lap,cond_sites = df_sites,dayshift = 0,v=q,Ndays_season = 90,title = paste0("all12sites",q*100))
# spatial_par_est(data_Lap = data_mod_Lap[1:3600,],cond_sites = df_sites,dayshift = 0,v=q,Ndays_season = 90,title = paste0("first40_12sites",q*100))
# spatial_par_est(data_Lap = data_mod_Lap[5401:9000,],cond_sites = df_sites,dayshift = 0,v=q,Ndays_season = 90,title = paste0("last40_12sites",q*100))

# load all three parameter estimates sf objects ----
load(paste0("data_processed/N9000_sequential2_AGG_all12sites",q*100,".RData"))
est_all <- est_all_sf

load(paste0("data_processed/N3600_sequential2_AGG_first40_12sites",q*100,".RData"))
est_allfirst40 <- est_all_sf

load(paste0("data_processed/N3600_sequential2_AGG_last40_12sites",q*100,".RData"))
est_alllast40 <- est_all_sf

# combine the first and the last 40 years
est_comb <- bind_rows(est_allfirst40 %>% mutate("Year_range"=c("1981_2020")),est_alllast40 %>% mutate("Year_range"=c("2041_2080")))

p1 <- ggplot(est_comb) + geom_density(aes(x=a,fill=Year_range),alpha=0.5)  +  scale_fill_manual(name="Year_range",labels=c(TeX("$1980-2020$"),
                                                                                                             TeX("$2041-2080$")),
                                                                                       values=c("black","#C11432")) + facet_wrap(~cond_site)


p2 <- ggplot(est_comb) + geom_density(aes(x=b,fill=Year_range),alpha=0.5)  +  scale_fill_manual(name="Years",labels=c(TeX("$1980-2020$"),
                                                                                                             TeX("$2041-2080$")),
                                                                                       values=c("black","#C11432"))+ facet_wrap(~cond_site)



ggsave(p1,width=10,height=5,filename = paste0("../Documents/acomb_12sites_q",q*100,".png"))
ggsave(p2,width=10,height=5,filename = paste0("../Documents/bcomb_12sitesq",q*100,".png"))

# plot estimates spatially for each site ------
# bind data together
est_comball <- bind_rows(est_all %>% mutate("Year_range"=c("1981_2080")),est_allfirst40 %>% mutate("Year_range"=c("1981_2020")),est_alllast40 %>% mutate("Year_range"=c("2041_2080")))
est_comball <- est_comball %>% mutate(Year_range=factor(Year_range,levels=c("1981_2080","1981_2020","2041_2080")))
#amin <- min(est_comball$a,na.rm=TRUE)
amin <- -0.4
amax <- 1
bmin <- 0
bmax <- 1
#amax <- max(est_comball$a,na.rm=TRUE)
#bmin <- min(est_comball$b,na.rm=TRUE)
#bmax <- max(est_comball$b,na.rm=TRUE)
for (i in 1:12) {
  pa <- tm_shape(est_comball %>% filter(cond_site==names(df_sites)[i]))  + tm_dots(col="a",palette="viridis",breaks=seq(amin,amax,length.out=11),size=0.3,colorNA="aquamarine",title=TeX("$\\alpha$"),textNA = "Conditioning site") + tm_facets(by= c("cond_site","Year_range"))+  tm_layout(legend.outside.size=0.3,asp=0.5,legend.text.size = 1,legend.title.size=1.5)
    tmap_save(pa,filename=paste0("../Documents/acomp_q",q*100,"_",names(df_sites[i]),".png"),width=10,height=6)
    pb <- tm_shape(est_comball %>% filter(cond_site==names(df_sites)[i]))  + tm_dots(col="b",palette="viridis",breaks=seq(bmin,bmax,length.out=11),size=0.3,colorNA="aquamarine",title=TeX("$\\beta$"),textNA = "Conditioning site") + tm_facets(by= c("cond_site","Year_range"))+  tm_layout(legend.outside.size=0.3,asp=0.5,legend.text.size = 1,legend.title.size=1.5)
    tmap_save(pb,filename=paste0("../Documents/bcomp_q",q*100,"_",names(df_sites[i]),".png"),width=10,height=6)
    }

# extend the analysis to a sliding window
# spatial_par_est(data_Lap = data_mod_Lap,cond_sites = df_sites,dayshift = 0,v=q,Ndays_season = 90,title = paste0("all12sites",q*100))
# spatial_par_est(data_Lap = data_mod_Lap[1:3600,],cond_sites = df_sites,dayshift = 0,v=q,Ndays_season = 90,title = paste0("first40_12sites",q*100))
# spatial_par_est(data_Lap = data_mod_Lap[5401:9000,],cond_sites = df_sites,dayshift = 0,v=q,Ndays_season = 90,title = paste0("last40_12sites",q*100))

time_bracket <- 900 # 10 years period is 900 summer days
block_size <- 3600 # 40 years is block size
# there are 7 40 year periods in 100 years
# run at first only for Birmingham to test
for (i in 1:((9000-(block_size-time_bracket))/time_bracket)) {
spatial_par_est(data_Lap = data_mod_Lap[((i-1)*time_bracket+1):(block_size+(i-1)*time_bracket),],cond_sites = df_sites,dayshift = 0,v=q,Ndays_season = 90,title = paste0("12sites_40years_period",i))
}
# load estimates
load(paste0("data_processed/N3600_sequential2_AGG_12sites_40years_period",1,".RData"))
est_all <- est_all_sf %>% mutate(period = as.character(1))
for (p in 2:((9000-(block_size-time_bracket))/time_bracket)) {
load(paste0("data_processed/N3600_sequential2_AGG_12sites_40years_period",p,".RData"))
est_all <- rbind(est_all,est_all_sf%>% mutate(period = as.character(p)))
  }
# plot alphas
tm_shape(est_all) + tm_dots(fill="a",size=1) + tm_facets(by=c("period"),nrow=1)
est_all$a[est_all$a<0] <- 0
point_size <- 0.3
lims <- c(0,1)
misscol <- "aquamarine"
panel_labels <- sapply(1:((9000-(block_size-time_bracket))/time_bracket),FUN=function(i){paste0((1980+(i-1)*time_bracket/90),"-",(1980+(block_size+(i-1)*time_bracket)/90))})
for (i in 1:12) {
#  pa <- tm_shape(uk_tmp1) + tm_dots(fill="a",fill.scale = tm_scale_continuous(limits=lims,values="viridis",value.na=misscol,label.na = "Conditioning site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\alpha$"))) + tm_facets(by=facet_var,nrow = nrow_facet) +  tm_layout(panel.labels = facet_label,legend.outside.size=legend_outside_size,asp=aspect,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE) + tm_title(text=title_map) 
  pa <- tm_shape(est_all %>% filter(cond_site==names(df_sites)[i]))  + tm_dots(fill="a",fill.scale = tm_scale_continuous(limits=lims,values="viridis",value.na=misscol,label.na = "Conditioning site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\alpha$"))) + tm_facets(by= c("period"),nrow = 1)+  tm_layout(panel.labels = panel_labels, legend.outside.size=0.3,asp=0.5,legend.text.size = 1,legend.title.size=1.5,legend.reverse = TRUE,legend.position = tm_pos_out("right","center",pos.h="left",pos.v="top")) + tm_title(names(df_sites)[i])
  tmap_save(pa,filename=paste0("../Documents/dependence_stationary/acomp_q",q*100,"_",names(df_sites[i]),".png"),width=12,height=4)
  pb <- tm_shape(est_all %>% filter(cond_site==names(df_sites)[i]))  + tm_dots(fill="b",fill.scale = tm_scale_continuous(limits=c(0,0.6),values="viridis",value.na=misscol,label.na = "Conditioning site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\beta$"))) + tm_facets(by= c("period"),nrow=1)+  tm_layout(panel.labels = panel_labels, legend.outside.size=0.3,asp=0.5,legend.text.size = 1,legend.title.size=1.5,legend.reverse = TRUE,legend.position = tm_pos_out("right","center",pos.h="left",pos.v="top")) + tm_title(names(df_sites)[i])
  tmap_save(pb,filename=paste0("../Documents/dependence_stationary/bcomp_q",q*100,"_",names(df_sites[i]),".png"),width=12,height=4)
# plot also difference between first and last panel
  
  }

names(panel_labels) <- c(as.character(1:7))
for (i in 1:12) {
pa <- est_all %>% filter(cond_site==names(df_sites)[i]) %>% ggplot() + geom_point(aes(x=period,y=a)) + geom_line(aes(x=period,y=a,group=res)) + ylab(TeX("$\\alpha$")) + xlab("40 year period T") + ggtitle(names(df_sites)[i])
ggsave(pa, filename = paste0("../Documents/dependence_stationary/acomp_plot_q",q*100,"_",names(df_sites[i]),".png"),width=10,height=10)
pa <- est_all %>% filter(cond_site==names(df_sites)[i]) %>% ggplot() + geom_boxplot(aes(y=a)) + facet_wrap("period",nrow = 1,labeller=as_labeller(panel_labels)) + ylab(TeX("$\\alpha$")) + xlab("40 year period T") + ggtitle(names(df_sites)[i])
ggsave(pa, filename = paste0("../Documents/dependence_stationary/acomp_boxplot_q",q*100,"_",names(df_sites[i]),".png"),width=10,height=6)

pb <- est_all %>% filter(cond_site==names(df_sites)[i]) %>% ggplot() + geom_point(aes(x=period,y=b)) + geom_line(aes(x=period,y=b,group=res)) + ylab(TeX("$\\beta$")) + xlab("40 year period T") + ggtitle(names(df_sites)[i])
ggsave(pb, filename = paste0("../Documents/dependence_stationary/bcomp_plot_q",q*100,"_",names(df_sites[i]),".png"),width=10,height=10)
pb <- est_all %>% filter(cond_site==names(df_sites)[i]) %>% ggplot() + geom_boxplot(aes(y=b)) + facet_wrap("period",nrow = 1,labeller=as_labeller(panel_labels))  + ylab(TeX("$\\beta$")) + xlab("40 year period") + ggtitle(names(df_sites)[i])
ggsave(pb, filename = paste0("../Documents/dependence_stationary/bcomp_boxplot_q",q*100,"_",names(df_sites[i]),".png"),width=10,height=6)
}


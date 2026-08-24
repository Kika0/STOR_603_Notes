library(tidyverse)
library(latex2exp)
library(gridExtra)
library(sf)
library(tmap)
library(evd)

theme_set(theme_bw())
theme_replace(
  panel.spacing = unit(2, "lines"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA) )
load("data_processed/spatial_helper.RData", verbose = TRUE)
load("data_processed/data_mod_Lap.RData",verbose = TRUE)
load("data_processed/temperature_data.RData",verbose = TRUE) # for data_mod and data_mod_Lap
folder_name <- "../Documents/spatial_chapter/"

chi_sites <- function(j,i,u=0.9){
  if (i==j) {return(NA)}
  sum(data_mod[,j][data_mod[,i]>quantile(data_mod[,i],u)]>quantile(data_mod[,i],u))/nrow(data_mod) /(1-u)
}

u <- 0.9
Birm_chi <- sapply(1:ncol(data_mod),FUN=chi_sites,i=df_sites[3,1],u=u)
Gla_chi <- sapply(1:ncol(data_mod),FUN=chi_sites,i=df_sites[3,2],u=u)
Lon_chi <- sapply(1:ncol(data_mod),FUN=chi_sites,i=df_sites[3,3],u=u)

# plot of chi as a map -------------------------------------------------------
tmp <- xyUK20_sf %>% mutate(Birm_chi,Gla_chi,Lon_chi)
cond_site_names <- names(df_sites)[1:3]
title_map <- ""
misscol <- "aquamarine"
legend_text_size <- 0.7
point_size <- 0.6
legend_title_size <- 1.2
lims <- c(0,1)
nrow_facet <- 1
p1 <- tm_shape(tmp) + tm_dots(fill="Birm_chi",fill.scale = tm_scale_continuous(limits=lims,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\chi$"))) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,legend.show=FALSE,frame=FALSE) + tm_title(text="Birmingham") 
p2 <- tm_shape(tmp) + tm_dots(fill="Gla_chi",fill.scale = tm_scale_continuous(limits=lims,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\chi$"))) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,legend.show = FALSE,frame=FALSE) + tm_title(text="Glasgow") 
p3 <- tm_shape(tmp) + tm_dots(fill="Lon_chi",fill.scale = tm_scale_continuous(limits=lims,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\chi$"))) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,legend.show = TRUE,frame=FALSE) + tm_title(text="London") 
tmap_save(tmap_arrange(p3,p1,p2,ncol=3),filename=paste0(folder_name,"chi_selected_sites",u*100,".png"),height=6,width=8)
tmap_save(tmap_arrange(p3,p1,p2,ncol=3),filename=paste0(folder_name,"chi_selected_sites_",u*100,".pdf"),height=6,width=8)


  



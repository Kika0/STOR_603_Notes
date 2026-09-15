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

# 1. plot of chi as a map -----------------------------------------------------
tmp <- xyUK20_sf %>% mutate(Birm_chi,Gla_chi,Lon_chi)
cond_site_names <- names(df_sites)[1:3]
title_map <- ""
misscol <- "aquamarine"
legend_text_size <- 0.7
point_size <- 0.6
legend_title_size <- 1.2
lims <- c(0,1)
nrow_facet <- 1
p1 <- tm_shape(tmp) + tm_dots(fill="Birm_chi",fill.scale = tm_scale_continuous(limits=lims,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\chi_u$"))) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,legend.show=FALSE,frame=FALSE) + tm_title(text="Birmingham") 
p2 <- tm_shape(tmp) + tm_dots(fill="Gla_chi",fill.scale = tm_scale_continuous(limits=lims,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\chi_u$"))) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,legend.show = FALSE,frame=FALSE) + tm_title(text="Glasgow") 
p3 <- tm_shape(tmp) + tm_dots(fill="Lon_chi",fill.scale = tm_scale_continuous(limits=lims,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\chi_u$"))) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,legend.show = TRUE,frame=FALSE) + tm_title(text="London") 
tmap_save(tmap_arrange(p3,p1,p2,ncol=3),filename=paste0(folder_name,"chi_selected_sites",u*100,".png"),height=6,width=8)
tmap_save(tmap_arrange(p3,p1,p2,ncol=3),filename=paste0(folder_name,"chi_selected_sites_",u*100,".pdf"),height=6,width=8)

# 2. plot of chi against distance from the conditioning site ------------------
Birmingham <- as.numeric(unlist(st_distance(xyUK20_sf[df_sites[3,1],],xyUK20_sf) %>%
                                  units::set_units(km)))
Birmingham[df_sites[3,1]] <- NA
Glasgow <- as.numeric(unlist(st_distance(xyUK20_sf[df_sites[3,2],],xyUK20_sf) %>%
                               units::set_units(km)))
Glasgow[df_sites[3,2]] <- NA
London <- as.numeric(unlist(st_distance(xyUK20_sf[df_sites[3,3],],xyUK20_sf) %>%
                              units::set_units(km)))
London[df_sites[3,3]] <- NA

tmp1 <- tmp %>% mutate(Birmingham,Glasgow,London)
tmp2 <- data.frame("cond_site_dist"=rep(NA,length(Birm_chi)*3),"chi"=rep(NA,length(Birm_chi)*3),"cond_site"=rep(NA,length(Birm_chi)*3))
tmp2$chi <- tmp1 %>% pivot_longer(cols=c(Birm_chi,Gla_chi,Lon_chi)) %>% pull(value)
tmp2$cond_site_dist <- tmp1 %>% pivot_longer(cols=c(Birmingham,Glasgow,London)) %>% pull(value)
tmp2$cond_site <- tmp1 %>% pivot_longer(cols=c(Birmingham,Glasgow,London)) %>% pull(name)
tmp2 <- tmp2 %>% mutate("cond_site"=factor(tmp2$cond_site,levels=c("London","Birmingham","Glasgow")))

ui <- c(0.9,0.95,0.99)
tmp3 <- data.frame("Birmingham"=numeric(),"Glasgow"=numeric(),"London"=numeric(),"u"=numeric())
for (u in ui) {
  Birm_chi <- sapply(1:ncol(data_mod),FUN=chi_sites,i=df_sites[3,1],u=u)
  Gla_chi <- sapply(1:ncol(data_mod),FUN=chi_sites,i=df_sites[3,2],u=u)
  Lon_chi <- sapply(1:ncol(data_mod),FUN=chi_sites,i=df_sites[3,3],u=u)
  tmp3 <- rbind(tmp3,data.frame("Birmingham"=Birm_chi,"Glasgow"=Gla_chi,"London"=Lon_chi,"u"=u) )
}
tmp4 <- tmp3 %>% pivot_longer(c(Birmingham,Glasgow,London),names_to = "cond_site",values_to="chi")
tmp5 <- tmp4 %>% mutate("cond_site_dist"=rep(tmp2$cond_site_dist,length(ui)))
tmp5 <- tmp5 %>% mutate("u"=factor(u,levels=c(ui))) %>% mutate(cond_site=factor(cond_site,levels=c("London","Birmingham","Glasgow")))
point_size <- 0.5
p <- ggplot(tmp5) + geom_point(aes(x=cond_site_dist,y=chi,col=u),size=point_size) + facet_wrap(~factor(cond_site)) + ylim(c(0,1)) + labs(x="Distance from conditioning site [km]",y=TeX("$\\chi_u$"),col=TeX("$u$")) + scale_color_manual(values=c("#009ADA","#033EB3","#090262"))
plot_name <- "chi_scatterplot_selected_sites"
ggsave(p,filename=paste0(folder_name,plot_name,".png"),width=10,height=3)
ggsave(p,filename=paste0(folder_name,plot_name,".pdf"),width=10,height=3)

# 3. chi outliers for Birmingham ----------------------------------------------
# try as a function of longitude and latitude
coastal_point <- function(grid) {
  sapply(1:nrow(grid),FUN = function(i) {sum(as.vector(st_distance(grid[i,],grid))<20500)<5 & grid$lat[i]+2*grid$lon[i]>48.5})
}

cp <- coastal_point(grid = xyUK20_sf) 
cp[df_sites[3,5]] <- FALSE
t3<- tm_shape(cbind(xyUK20_sf,data.frame(cp))) + tm_dots(fill="cp",fill.scale = tm_scale_categorical(values=c("FALSE"="black","TRUE"="#C11432")),size=point_size, fill.legend = tm_legend(title="")) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,legend.show=FALSE,frame=FALSE) + tm_title(text="East Coast")
tmp6 <- data.frame("cp"=cp,"Birm_chi"=Birm_chi,"Birm_dist"=Birmingham)
p <- ggplot(tmp6) + geom_point(aes(x=Birm_dist,y=Birm_chi,col=cp)) + scale_color_manual(values=c("black","#C11432"),labels=c("East coast", "Not East coast"),name="") + labs(x="Distance from conditioning site [km]",y=TeX("$\\chi_u$"))
p1 <- grid.arrange(t3,p,ncol=2)
plot_name <- "chi_scatterplot_selected_sites"
ggsave(p,filename=paste0(folder_name,plot_name,".png"),width=10,height=3)
ggsave(p,filename=paste0(folder_name,plot_name,".pdf"),width=10,height=3)


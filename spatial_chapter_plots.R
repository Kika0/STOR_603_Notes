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

c13 <- c(
  "#009ADA", "#C11432", # red
           "green4",
           "#6A3D9A", # purple
           "#FF7F00", # orange
           "black", "gold1",
           "#FB9A99", # lt pink
           "gray70", 
           "darkturquoise", "green1", 
           "darkorange4","#F6A7B8"
)

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
point_size <- 0.3
cp <- coastal_point(grid = xyUK20_sf) 
cp[df_sites[3,5]] <- FALSE
t3<- tm_shape(cbind(xyUK20_sf,data.frame(cp))) + tm_dots(fill="cp",fill.scale = tm_scale_categorical(values=c("FALSE"="black","TRUE"="#C11432")),size=point_size, fill.legend = tm_legend(title="")) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,legend.show=FALSE,frame=FALSE) + tm_title(text="East coast")
tmp6 <- data.frame("cp"=cp,"Birm_chi"=Birm_chi,"Birm_dist"=Birmingham)
p <- ggplot(tmp6) + geom_point(aes(x=Birm_dist,y=Birm_chi,col=cp)) + scale_color_manual(values=c("black","#C11432"),labels=c( "Not East coast","East coast"),name="") + labs(x="Distance from conditioning site [km]",y=TeX("$\\chi_u$"))
p1 <- grid.arrange(tmap_grob(t3),p,ncol=2,nrow=1,widths=c(1.2,2))
plot_name <- "chi_scatterplot_Birmingham_east_coast"
ggsave(p1,filename=paste0(folder_name,plot_name,".png"),width=10,height=4)
ggsave(p1,filename=paste0(folder_name,plot_name,".pdf"),width=10,height=4)

# 4. illustrate maps of conditioning sites
uk_diag <- xyUK20_sf %>% mutate(siteID=as.numeric(1:nrow(xyUK20_sf)))
sites_index <- as.numeric(df_sites[3,1:12])
sites_name<- names(df_sites)[1:12]
# plot these points on a map
uk_diag <- uk_diag %>% mutate(sites12=NA)
uk_diag$sites12[sites_index] <- sites_name
legend_text_size <- 0.72
Birm_point <- xyUK20_sf[sites_index[1],]
t1 <- tm_shape(uk_diag) + tm_dots("sites12",size=0.5,fill.scale = tm_scale_categorical(label.na="",values=c13[1:12]),fill.legend = tm_legend(title="",frame=FALSE) ) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,legend.show=TRUE,frame=FALSE) + tm_title(text="12 conditioning sites") + tm_shape(Birm_point) + tm_dots(shape=1,size=1)

# identify diagonal sites Birmingham-London
site_start <- (df_sites %>% dplyr::select(Birmingham))[3,1]
site_end <- (df_sites %>% dplyr::select(London))[3,1]
uk_diag <- xyUK20_sf %>% mutate(siteID=as.numeric(1:nrow(xyUK20_sf)))
sites_index_diagonal <- c(192,174,156,157,137,138,118,100) # first is Birmingham and last is London
site_name_diagonal <- c("Birmingham", paste0("diagonal",1:(length(sites_index_diagonal)-2)),"London")
# plot these points on a map
uk_diag <- uk_diag %>% mutate(sites_diagonal=factor(case_match(siteID,c(site_start) ~ "Birmingham",c(site_end)~"London",sites_index_diagonal[2:(length(site_name_diagonal)-1)]~"diagonal_sites")))
t2 <- tm_shape(uk_diag) + tm_dots("sites_diagonal",size=0.5,fill.scale = tm_scale_categorical(values=c("Birmingham"="#C11432","London" = "#009ADA", "diagonal_sites" = "#FDD10A"))) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,legend.show=FALSE,frame=FALSE) + tm_title(text="Birmingham to London") 

# identify diagonal sites Birmingham-Cromer
site_start <- (df_sites %>% dplyr::select(Birmingham))[3,1]
site_end <- (df_sites %>% dplyr::select(Cromer))[3,1]
uk_diag <- xyUK20_sf %>% mutate(siteID=as.numeric(1:nrow(xyUK20_sf)))
sites_index_diagonal <- c(192,193,194,195,196,197,218,219,220,221,242,263) # first is Birmingham and last is Cromer
site_name_diagonal <- c("Birmingham", paste0("diagonal",1:(length(sites_index_diagonal)-2)),"Cromer")
# plot these points on a map
uk_diag <- uk_diag %>% mutate(sites_diagonal=factor(case_match(siteID,c(site_start) ~ "Birmingham",c(site_end)~"Cromer",sites_index_diagonal[2:(length(site_name_diagonal)-1)]~"diagonal_sites")))
t3 <- tm_shape(uk_diag) + tm_dots("sites_diagonal",size=0.5,fill.scale = tm_scale_categorical(values=c("Birmingham"="#C11432","Cromer" = "#009ADA", "diagonal_sites" = "#FDD10A"))) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,legend.show=FALSE,frame=FALSE) + tm_title(text="Birmingham to Cromer") 

plot_name <- "cond_sites_illustration"
tmap_save(tmap_arrange(t1,t2,t3,ncol=3),filename=paste0(folder_name,plot_name,".png"),height=6,width=8)
tmap_save(tmap_arrange(t1,t2,t3,ncol=3),filename=paste0(folder_name,plot_name,".pdf"),height=6,width=8)

# 5. plot of illustration of marginal transform -------------------------------
Birm_temp <- data_mod_temp[,(df_sites %>% dplyr::select(Birmingham))[3,1]]
Gla_temp <- data_mod_temp[,(df_sites %>% dplyr::select(Glasgow))[3,1]]
tmp <- data.frame(Birm_temp,Gla_temp) %>% mutate(year5=factor(rep(seq(1980,2075,by=5),each=90*5)),obsmod="CPM_data")
# include also observed data
Birm_temp <- data_obs_all[,(df_sites %>% dplyr::select(Birmingham))[3,1]]
Gla_temp <- data_obs_all[,(df_sites %>% dplyr::select(Glasgow))[3,1]]
tmp1 <- data.frame(Birm_temp,Gla_temp) %>% mutate(year5=factor(rep(seq(1960,2020,by=5),each=92*5)[1:length(Birm_temp)]),obsmod="observed")
tmp2 <- rbind(tmp,tmp1) %>% mutate(year5=factor(year5,levels=seq(1960,2075,by=5)))
names(tmp2)[1:2] <- c("Birmingham","Glasgow")
tmp3 <- tmp2 %>% pivot_longer(c(Birmingham,Glasgow),names_to="site",values_to = "temp")
p <- ggplot(tmp3) + geom_boxplot(aes(x=year5,y=temp,fill=obsmod)) + facet_wrap(~site) + scale_fill_manual(values=c("#C11432","black"),labels=c("CPM data","Observed data")) + labs(fill="",x="",y=TeX("Temperature ($^\\circ C$)"))
plot_name <- "CPM_observed_Birmingham_Glasgow"
ggsave(p,filename=paste0(folder_name,plot_name,".png"),width=16,height=4)
ggsave(p,filename=paste0(folder_name,plot_name,".pdf"),width=16,height=4)

Birm_temp <- data_mod_Lap[,(df_sites %>% dplyr::select(Birmingham))[3,1]]
Gla_temp <- data_mod_Lap[,(df_sites %>% dplyr::select(Glasgow))[3,1]]
tmp <- data.frame(Birm_temp,Gla_temp) %>% mutate(year5=factor(rep(seq(1980,2075,by=5),each=90*5)),obsmod="Laplace")
# include also observed data
Birm_temp <- data_mod_Lap_star[,(df_sites %>% dplyr::select(Birmingham))[3,1]]
Gla_temp <- data_mod_Lap_star[,(df_sites %>% dplyr::select(Glasgow))[3,1]]
tmp1 <- data.frame(Birm_temp,Gla_temp) %>% mutate(year5=factor(rep(seq(1980,2075,by=5),each=90*5)),obsmod="double_Laplace")
tmp2 <- rbind(tmp,tmp1) %>% mutate(year5=factor(year5,levels=seq(1960,2075,by=5)))
names(tmp2)[1:2] <- c("Birmingham","Glasgow")
tmp3 <- tmp2 %>% pivot_longer(c(Birmingham,Glasgow),names_to="site",values_to = "temp")
tmp3 <- tmp3 %>% mutate(obsmod=factor(obsmod,levels=c("Laplace", "double_Laplace")))
p <- ggplot(tmp3) + geom_boxplot(aes(x=year5,y=temp,fill=obsmod)) + facet_wrap(~site) + scale_fill_manual(values=c("#009ADA","#66A64F"),labels=c("Laplace","Double Laplace"))  + labs(fill="",x="",y="Temperature (Laplace scale)")
plot_name <- "CPM_Laplace_Birmingham_Glasgow"
ggsave(p,filename=paste0(folder_name,plot_name,".png"),width=16,height=4)
ggsave(p,filename=paste0(folder_name,plot_name,".pdf"),width=16,height=4)

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
df_sites <- df_sites[,1:12]

c13 <- c(
  "#C11432", "#74499C", # red
           "#66A64F",
"#0047AB",
           "#FF7F00", # orange
           "black", "#007191",
           "#FB9A99", # lt pink
           "#009ADA",
           "darkturquoise", "green1", 
           "darkorange4","#F6A7B8"
)

chi_sites <- function(j,i,u=0.9){
  if (i==j) {return(NA)}
  sum(data_mod[,j][data_mod[,i]>quantile(data_mod[,i],u)]>quantile(data_mod[,i],u))/nrow(data_mod) /(1-u)
}

chi_bar_sites <- function(j,i,u=0.9){
  if (i==j) {return(NA)}
P_uv <-  sum(data_mod[,j][data_mod[,i]>quantile(data_mod[,i],u)]>quantile(data_mod[,i],u))/nrow(data_mod)
#  chiu <- 2 - log(P_uv)/log(u)
return(2*log(1-u)/log(P_uv)-1)
  
  ## calculate chibar
#  P_uv <- mean((U > u)*(V > u))
#  chibaru <- 2*log(1 - u)/log(P_uv) - 1
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

# repeat for chi bar ---------------------------------------------------------
u <- 0.9
Birm_chi <- sapply(1:ncol(data_mod),FUN=chi_bar_sites,i=df_sites[3,1],u=u)
Gla_chi <- sapply(1:ncol(data_mod),FUN=chi_bar_sites,i=df_sites[3,2],u=u)
Lon_chi <- sapply(1:ncol(data_mod),FUN=chi_bar_sites,i=df_sites[3,3],u=u)

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
  tmap_save(tmap_arrange(p3,p1,p2,ncol=3),filename=paste0(folder_name,"chi_bar_selected_sites",u*100,".png"),height=6,width=8)
  tmap_save(tmap_arrange(p3,p1,p2,ncol=3),filename=paste0(folder_name,"chi_bar_selected_sites_",u*100,".pdf"),height=6,width=8)
  
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
    Birm_chi <- sapply(1:ncol(data_mod),FUN=chi_bar_sites,i=df_sites[3,1],u=u)
    Gla_chi <- sapply(1:ncol(data_mod),FUN=chi_bar_sites,i=df_sites[3,2],u=u)
    Lon_chi <- sapply(1:ncol(data_mod),FUN=chi_bar_sites,i=df_sites[3,3],u=u)
    tmp3 <- rbind(tmp3,data.frame("Birmingham"=Birm_chi,"Glasgow"=Gla_chi,"London"=Lon_chi,"u"=u) )
  }
  tmp4 <- tmp3 %>% pivot_longer(c(Birmingham,Glasgow,London),names_to = "cond_site",values_to="chi")
  tmp5 <- tmp4 %>% mutate("cond_site_dist"=rep(tmp2$cond_site_dist,length(ui)))
  tmp5 <- tmp5 %>% mutate("u"=factor(u,levels=c(ui))) %>% mutate(cond_site=factor(cond_site,levels=c("London","Birmingham","Glasgow")))
  point_size <- 0.5
  p <- ggplot(tmp5) + geom_point(aes(x=cond_site_dist,y=chi,col=u),size=point_size) + facet_wrap(~factor(cond_site)) + ylim(c(0,1)) + labs(x="Distance from conditioning site [km]",y=TeX("$\\bar{\\chi}_u$"),col=TeX("$u$")) + scale_color_manual(values=c("#009ADA","#033EB3","#090262"))
  plot_name <- "chi_bar_scatterplot_selected_sites"
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
plot_name <- "chi_bar_scatterplot_Birmingham_east_coast"
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
t1 <- tm_shape(uk_diag) + tm_dots("sites12",size=0.5,fill.scale = tm_scale_categorical(label.na="",values=c13[1:12]),fill.legend = tm_legend(title="",frame=FALSE) ) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=FALSE,legend.show=TRUE,frame=FALSE) + tm_title(text="12 conditioning sites") + tm_shape(Birm_point) + tm_dots(shape=1,size=1)

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
t3 <- tm_shape(uk_diag) + tm_dots("sites_diagonal",size=0.5,fill.scale = tm_scale_categorical(values=c("Birmingham"="#C11432","Cromer" = "#66A64F", "diagonal_sites" = "#FDD10A"))) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,legend.show=FALSE,frame=FALSE) + tm_title(text="Birmingham to Cromer") 

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
data_mod_Lap_star <- as.data.frame((data_mod_Lap %>% apply(c(2),FUN=row_number))/(nrow(data_mod_Lap)+1)) %>% apply(c(1,2),FUN=unif_laplace_pit) %>% as.data.frame()
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


# 10.1 plot comparing Laplace and double Laplace Birmingham and Glasgow ------
Birm_temp_Lap <- data_mod_Lap[,(df_sites %>% dplyr::select(Birmingham))[3,1]]
Gla_temp_Lap <- data_mod_Lap[,(df_sites %>% dplyr::select(Glasgow))[3,1]]
Birm_temp <- data_mod_Lap_star[,(df_sites %>% dplyr::select(Birmingham))[3,1]]
Gla_temp <- data_mod_Lap_star[,(df_sites %>% dplyr::select(Glasgow))[3,1]]
tmp4 <- rbind(data.frame("site"="Birmingham","data_Lap"=Birm_temp_Lap,"double_Lap"=Birm_temp),data.frame("site"="Glasgow","data_Lap"=Gla_temp_Lap,"double_Lap"=Gla_temp))
p <- ggplot(tmp4) + geom_point(aes(x=double_Lap,y=data_Lap,col=site)) + geom_abline(slope=1,intercept=0,linetype="dashed") + scale_color_manual(values=c("#C11432", "#74499C")) + labs(x="Double Laplace",y="Laplace",col="Conditioning site")
plot_name <- "DoubleLaplace_Laplace_Birmingham_Glasgow"
ggsave(p,filename=paste0(folder_name,plot_name,".png"),width=6,height=4)
ggsave(p,filename=paste0(folder_name,plot_name,".pdf"),width=6,height=4)

# 5.2 try to find keef estimates
load("data_processed/est_all.RData",verbose=TRUE)
load("data_processed/paramest.rda",verbose=TRUE)
load("data_processed/alphaspatial.RData",verbose=TRUE)

# plot alpha and beta against distance
# load estimates with keef constraints
load("data_processed/N9000_sequential2_AGG_all12sites90.RData",verbose=TRUE)
sites_selection <- c("Birmingham","Glasgow","London", "Inverness","Lancaster","Cromer")
#sites_selection <- names(df_sites)[!(names(df_sites) %in% c("Birmingham","Glasgow","London", "Inverness","Lancaster","Cromer"))]
tmp <- data.frame("b_keef"=est_all_sf$b,"given"=est_all_sf$given,"res"=est_all_sf$res,"cond_site"=est_all_sf$cond_site)
tmp <- tmp %>% mutate(cond_site=factor(cond_site)) %>% dplyr::filter(cond_site %in% sites_selection) |> arrange(cond_site)

# repeat estimates no keef constraints
#spatial_par_est(data_Lap = data_mod_Lap,cond_sites = df_sites,dayshift = 0,v=q,Ndays_season = 90,title = paste0("all12sites_nokeef",q*100),keef=0)
# load these estimates
load("data_processed/N9000_sequential2_AGG_all12sites_nokeef90.RData",verbose=TRUE)
df_sites <- df_sites[ , order(names(df_sites))]
# calculate distance from the conditioning site
get_cond_site_distance <- function(i, grid = xyUK20_sf,sites=df_sites) {
  cond_site_index <- get_site_index(i)
  dist_cond_site <- as.numeric(unlist(st_distance(grid[cond_site_index,],grid)%>%
                                        units::set_units(km)))
  return(dist_cond_site)
}
dij <- unlist(lapply(1:ncol(df_sites),FUN=get_cond_site_distance))
est_all_sf <- est_all_sf %>% dplyr::filter(cond_site %in% names(df_sites)) |> arrange(cond_site) %>% mutate("dij"=dij)
tmp1 <- est_all_sf %>% dplyr::filter(cond_site %in% sites_selection) |> arrange(cond_site)

tmp2 <- data.frame("a"=tmp1$a,"cond_site"=tmp1$cond_site,"b"=tmp1$b,"b_keef"=tmp$b_keef,"dij"=tmp1$dij)
point_size <- 0.3
# split by two sets of three
p1 <- ggplot(tmp2 |> dplyr::filter(cond_site %in% sites_selection[1:3])) + geom_point(aes(x=dij,y=a,col=cond_site),size=point_size) + labs(x="Distance from conditioning site [km]",y=TeX("$\\alpha$"),col=c("Conditioning site")) + scale_color_manual(values=c13[1:ncol(df_sites)][names(df_sites) %in% sites_selection[1:3]])
p2 <- ggplot(tmp2|> dplyr::filter(cond_site %in% sites_selection[4:6])) + geom_point(aes(x=dij,y=a,col=cond_site),size=point_size) + labs(x="Distance from conditioning site [km]",y=TeX("$\\alpha$"),col=c("Conditioning site")) + scale_color_manual(values=c13[1:ncol(df_sites)][names(df_sites) %in% sites_selection[4:6]])
margin = theme(plot.margin = unit(c(0.5,1,0.5,0.5), "cm"))
p3 <- grid.arrange(grobs = lapply(list(p1,p2), "+", margin),ncol=2)
p4 <- ggplot(tmp2 |> dplyr::filter(cond_site %in% sites_selection[1:3])) + geom_point(aes(x=dij,y=b,col=cond_site),size=point_size) + labs(x="Distance from conditioning site [km]",y=TeX("$\\beta$"),col=c("Conditioning site")) + scale_color_manual(values=c13[1:ncol(df_sites)][names(df_sites) %in% sites_selection[1:3]]) + ylim(c(0,1))
p5 <- ggplot(tmp2 |> dplyr::filter(cond_site %in% sites_selection[4:6])) + geom_point(aes(x=dij,y=b,col=cond_site),size=point_size) + labs(x="Distance from conditioning site [km]",y=TeX("$\\beta$"),col=c("Conditioning site")) + scale_color_manual(values=c13[1:ncol(df_sites)][names(df_sites) %in% sites_selection[4:6]]) + ylim(c(0,1))
p6 <- grid.arrange(grobs = lapply(list(p4,p5), "+", margin),ncol=2)
plot_name <- "alpha_scatterplot_selected_sites" # change for different site selection
ggsave(p3,filename=paste0(folder_name,plot_name,".png"),width=13,height=4)
ggsave(p3,filename=paste0(folder_name,plot_name,".pdf"),width=13,height=4)
plot_name <- "beta_scatterplot_selected_sites" # change for different site selection
ggsave(p6,filename=paste0(folder_name,plot_name,".png"),width=13,height=4)
ggsave(p6,filename=paste0(folder_name,plot_name,".pdf"),width=13,height=4)
# plot alpha and beta values on a map
pa <- pb <- list()
for (i in 1:length(sites_selection)) {
title_map <- sites_selection[i]
misscol <- "aquamarine"
  legend_text_size <- 0.7
  point_size <- 0.4
  legend_title_size <- 0.9
  limsa <- c(0,1)
  limsb <- c(0,1)
  nrow_facet <- 1
  a <- tmp2 %>% dplyr::filter(cond_site %in% sites_selection[i]) %>% pull(a)
  b <- tmp2 %>% dplyr::filter(cond_site %in% sites_selection[i]) %>% pull(b)
  estsf <- cbind(xyUK20_sf, data.frame("a"=a,"b"=b))

  pa[[i]] <- tm_shape(estsf) + tm_dots(fill="a",fill.scale = tm_scale_continuous(limits=limsa,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\alpha$"))) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE,legend.show=(i==2)) + tm_title(text=title_map) 
  pb[[i]] <- tm_shape(estsf) + tm_dots(fill="b",fill.scale = tm_scale_continuous(limits=limsb,values="brewer.blues",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\beta$"))) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE,legend.show=(i==2)) + tm_title(text=title_map) 

}
pall <- tmap_arrange(pa[[1]],pb[[1]],pa[[2]],pb[[2]],pa[[3]],pb[[3]],pa[[4]],pb[[4]],pa[[5]],pb[[5]],pa[[6]],pb[[6]],ncol=4)
plot_name <- "alpha_beta_selected_sites" # change for different site selection
tmap_save(pall,filename=paste0(folder_name,plot_name,".png"),height=15,width=13)
tmap_save(pall,filename=paste0(folder_name,plot_name,".pdf"),height=15,width=13)

# compare beta without and with Keef constraints using a density plot
tmp3 <- rbind(data.frame(a=na.omit(tmp2$a),b=na.omit(tmp2$b),"method"="aoriginal"),data.frame(a=na.omit(tmp2$a),b=na.omit(tmp2$b_keef),"method"="new")) %>% mutate("iteration"=rep(1:length(na.omit(tmp2$a)),2),"cond_site"=rep(tmp2$cond_site[!is.na(tmp2$a)],2))
plot_ab <- function(tmp) { ggplot(tmp) + 
    geom_density(aes(x=b,fill=method),alpha=0.5)+
    xlab(TeX("${\\beta}$")) +
    ylab("Kernel smoothed density") + 
    scale_fill_manual(values = c("aoriginal" = "#009ADA", "new" = "#C11432"),labels = c("Without constraints","With Keef constraints"))  + theme(axis.text.y = element_text(angle = 90, vjust = 0.5),legend.title=element_blank()) + labs(col="") + facet_wrap(~cond_site,ncol=3)  + xlim(c(0,1))
}
p <- plot_ab(tmp=tmp3)
plot_name <- "plot_density_b_nokeef_keef_selected" # change for different site selection

ggsave(p,filename=paste0(folder_name,plot_name,".png"),width=10,height=6) 
ggsave(p,filename=paste0(folder_name,plot_name,".pdf"),width=10,height=6) 


# plot alpha, beta and beta with Keef constraints on a map hedgehog plot
plot_ab <- function(tmp) { ggplot(tmp) + 
    geom_line(aes(x=a,y=b,group=iteration),linewidth=0.1)+
    geom_point(aes(x=a,y=b,col=method),size=0.3) +
    xlab(TeX("${\\alpha}$")) +
    ylab(TeX("${\\beta}$")) + 
    scale_color_manual(values = c("aoriginal" = "#009ADA", "new" = "#C11432"),labels = c("Without constraints","With Keef constraints")) + coord_fixed() + theme(axis.text.y = element_text(angle = 90, vjust = 0.5)) + labs(col="") + facet_wrap(~cond_site,ncol=3) +ylim(c(0,1)) + xlim(c(0,1)) 
}
p <- plot_ab(tmp=tmp3)
plot_name <- "plot_ab_nokeef_keef_selected" # change for different site selection

ggsave(p,filename=paste0(folder_name,plot_name,".png"),width=10,height=6) 
#ggsave(p,filename=paste0(folder_name,plot_name,".pdf"),width=10,height=6) 
library(png)
# load .png image
img <- readPNG(paste0(folder_name,plot_name,".png"))
# save as pdf
pdf(paste0(folder_name,plot_name,".pdf"),width=10,height=6)
grid::grid.raster(img)
dev.off()

# map beta values with Keef constraints
pb <- list()
for (i in 1:length(sites_selection)) {
  title_map <- sites_selection[i]
  misscol <- "aquamarine"
  legend_text_size <- 0.7
  point_size <- 0.4
  legend_title_size <- 0.9
  limsa <- c(0,1)
  limsb <- c(0,1)
  nrow_facet <- 1
  b <- tmp2 %>% dplyr::filter(cond_site %in% sites_selection[i]) %>% pull(b_keef)
  estsf <- cbind(xyUK20_sf, data.frame("a"=a,"b"=b))
  pb[[i]] <- tm_shape(estsf) + tm_dots(fill="b",fill.scale = tm_scale_continuous(limits=limsb,values="brewer.blues",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\beta$"))) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE,legend.show=(i==length(sites_selection))) + tm_title(text=title_map) 
  
}
pall <- tmap_arrange(pb[[1]],pb[[2]],pb[[3]],pb[[4]],pb[[5]],pb[[6]],ncol=6)
plot_name <- "beta_keef_selected_sites" # change for different site selection
tmap_save(pall,filename=paste0(folder_name,plot_name,".png"),height=6,width=16)
tmap_save(pall,filename=paste0(folder_name,plot_name,".pdf"),height=6,width=16)
# plot also for the rest of the sites by changing site selection above

# estimate parameters at temporal lags
#spatial_par_est(data_Lap = data_mod_Lap,cond_sites = df_sites,dayshift = c(-3:3),v=q,Ndays_season = 90,title = paste0("all12sites_keef_temporal",q*100),keef=c(1,2),res_margin = "Normal")

#load estimates
load("data_processed/N9000_sequential2_Normal_all12sites_keef_temporal90.RData",verbose=TRUE)

# map alpha and beta at different temporal lags for Birmingham ----------------
pa <- pb <- list()
tit_day <- c("days","days","day","days","day","days","days")
for (i in 1:length(unique(est_all_sf$tau))) {
  title_map <- paste0(unique(est_all_sf$tau)[i]," ",tit_day[i])
  misscol <- "aquamarine"
  legend_text_size <- 0.7
  point_size <- 0.4
  legend_title_size <- 0.9
  limsa <- c(0,1)
  limsb <- c(0,1)
  nrow_facet <- 1
  a <- est_all_sf %>% dplyr::filter(cond_site %in% sites_selection[1],tau %in% unique(est_all_sf$tau)[i]) %>% pull(a)
  b <- est_all_sf %>% dplyr::filter(cond_site %in% sites_selection[1],tau %in% unique(est_all_sf$tau)[i]) %>% pull(b)
  estsf <- cbind(xyUK20_sf, data.frame("a"=a,"b"=b))
  pa[[i]] <- tm_shape(estsf) + tm_dots(fill="a",fill.scale = tm_scale_continuous(limits=limsb,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\alpha$"))) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE,legend.show=(i == length(unique(est_all_sf$tau)))) + tm_title(text=title_map) 
  pb[[i]] <- tm_shape(estsf) + tm_dots(fill="b",fill.scale = tm_scale_continuous(limits=limsb,values="brewer.blues",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\beta$"))) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE,legend.show=(i == length(unique(est_all_sf$tau)))) + tm_title(text=title_map) 
  
}
pall <- tmap_arrange(pa[[1]],pa[[2]],pa[[3]],pa[[4]],pa[[5]],pa[[6]],pa[[7]],pb[[1]],pb[[2]],pb[[3]],pb[[4]],pb[[5]],pb[[6]],pb[[7]],ncol=7)
plot_name <- "alpha_beta_keef_temporal_Birmingham"
tmap_save(pall,filename=paste0(folder_name,plot_name,".png"),height=12,width=18)
tmap_save(pall,filename=paste0(folder_name,plot_name,".pdf"),height=12,width=18)

# calculate max tau and corresponding alpha value ----------------------------
est_all_sf <- est_all_sf |> dplyr::filter(cond_site %in% sites_selection) 
alpha_max <- tau_max <- cond_site <- c()
for (i in 1:length(sites_selection)) {
  for (j in 1:nrow(xyUK20_sf)) {
 alpha_max <- append(alpha_max,est_all_sf |> st_drop_geometry() %>% dplyr::filter(cond_site==sites_selection[i],site_index==j) %>% dplyr::select(a) %>%  pull() %>% max())
 tau_max <- append(tau_max,if (is.na(alpha_max[length(alpha_max)]))  {factor(NA,levels=unique(est_all_sf$tau))} else{unique(est_all_sf$tau)[est_all_sf |> st_drop_geometry() %>% filter(cond_site==sites_selection[i],site_index==j) %>% dplyr::select(a) %>%  pull() %>% which.max()]})
 cond_site <- append(cond_site,sites_selection[i])
  }
}
#    est_all_max <- cbind(est_all_max,data.frame(max_map, max_time))
#   names(est_all_max)[c(ncol(est_all_max)-1,ncol(est_all_max))]<- c(condmodel_params[k],paste0(condmodel_params[k],"_tau"))

#  est_all_m[(n_sites*(i-1)+1):(n_sites*i),] <- est_all_max %>% mutate(cond_site=names(df_sites)[i])
tmp4 <- cbind(est_all_sf |> dplyr::filter(tau==0) |> dplyr::select(),data.frame("alpha_max"=alpha_max, "tau_max" = tau_max,"cond_site"=cond_site)) 
# map max alpha and corresponding tau
palpha_max <- ptau_max <- list()
for (i in 1:length(sites_selection)) {
  title_map <- sites_selection[i]
  misscol <- "aquamarine"
  legend_text_size <- 0.7
  point_size <- 0.4
  legend_title_size <- 0.9
  limsa <- c(0,1)
  limsb <- c(0,1)
  nrow_facet <- 1
  alpha_max <- tmp4 %>% dplyr::filter(cond_site %in% sites_selection[i]) %>% pull(alpha_max)
  tau_max <- tmp4 %>% dplyr::filter(cond_site %in% sites_selection[i]) %>% pull(tau_max)
  estsf <- cbind(xyUK20_sf, data.frame("alpha_max"=alpha_max,"tau_max"=tau_max))
  palpha_max[[i]] <- tm_shape(estsf) + tm_dots(fill="alpha_max",fill.scale = tm_scale_continuous(limits=limsb,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\alpha$"))) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE,legend.show=(i == length(sites_selection))) + tm_title(text=title_map) 
  ptau_max[[i]] <- tm_shape(estsf) + tm_dots(fill="tau_max",fill.scale = tm_scale_ordinal(values="-brewer.rd_bu",labels=levels(tmp4$tau),value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\tau$"))) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE,legend.show=(i == length(sites_selection))) + tm_title(text=title_map) 
  
}
pall <- tmap_arrange(palpha_max[[1]],palpha_max[[2]],palpha_max[[3]],palpha_max[[4]],palpha_max[[5]],palpha_max[[6]],ptau_max[[1]],ptau_max[[2]],ptau_max[[3]],ptau_max[[4]],ptau_max[[5]],ptau_max[[6]],ncol=6)
plot_name <- "alpha_max_tau_keef_temporal_Birmingham"
tmap_save(pall,filename=paste0(folder_name,plot_name,".png"),height=12,width=16)
tmap_save(pall,filename=paste0(folder_name,plot_name,".pdf"),height=12,width=16)

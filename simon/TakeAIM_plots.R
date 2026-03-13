library(tidyverse)
library(tmap)
library(tmaptools)
library(sf)
library(lubridate)
library(gridExtra)
library(LaplacesDemon)
theme_set(theme_bw())
theme_replace(
  panel.spacing = unit(2, "lines"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA) )
folder_name <- "../Documents/TakeAIM2026/"

# load observed data
#source("spatial_parameter_estimation.R") # for spatial_par_est function
load("data_processed/temperature_data.RData",verbose = TRUE)
load("data_processed/spatial_helper.RData", verbose = TRUE)
source("simon/P2q_function_helpers.R")

# load to determine index
load("../luna/kristina/P2q/ukgd_cpm85_5k_x84y20_MSp2q.RData", verbose = TRUE)
load("../luna/kristina/MSdata01/ukgd_cpm85_5k_x100y24_MSdata01.RData")
# add column for date
data01_obs <- data01 %>% rowid_to_column() %>% filter(class=="obs") %>% mutate("date_obs"=seq(ymd('1960-01-01'),ymd('2023-12-31'),by='days') )
# filter date
july3_1976 <- data01_obs %>% dplyr::filter(date_obs==lubridate::ymd("1976-07-03")) %>% pull(rowid)
# find future dates that match from model data
which.min(abs(data01$time[data01$class=="mod"]- (data01_obs$time[july3_1976] +50)))
july3_2026 <- which.min(abs(data01$time[data01$class=="mod"]- (data01_obs$time[july3_1976] +50))) + nrow(data01_obs)
july3_2076 <- which.min(abs(data01$time[data01$class=="mod"]- (data01_obs$time[july3_1976] +100))) + nrow(data01_obs)

# compare P2q of the observed 22846 and model 38362 data
# load("../luna/kristina/MSGpdParam/ukgd_cpm85_5k_x84y20.MSGpdParam.2025-02-26.RData", verbose = TRUE)
# gpdpar[c(22846,38362),]
# qgam.p2q.fn[[22845]](data01$u[22845])
# qgam.p2q.fn[[38362]](data01$u[22846])
# data01[c(22846,38362),]


# load below threshold functions ---------------------------------------------
folder_p2q_below94 <- "../luna/kristina/P2q/"
files <- list.files(folder_p2q_below94)
P2q_sites1 <-P2q_sites2 <-P2q_sites3 <- list() #create empty list
# only subset for x in x20 and y in y20
files_subset <- sapply(1:nrow(xyUK20_sf),function(i){paste0("ukgd_cpm85_5k_x",xyUK20_sf$x[i],"y",xyUK20_sf$y[i],"_MSp2q.RData")})
#loop through the files
files_subset1 <- files_subset[files_subset %in% files]
# could take only x and y divisible by 4 to subset and speed up data loading
for (i in 1:length(files_subset1)) {
  print(files_subset1[i])
  load(paste0(folder_p2q_below94, files_subset1[i]))
  P2q_sites1[[i]] <- qgam.p2q.fn[[july3_1976]] #add files to list position
  P2q_sites2[[i]] <- qgam.p2q.fn[[july3_2026]] #add files to list position
  P2q_sites3[[i]] <- qgam.p2q.fn[[july3_2076]] #add files to list position
  
}

# load GPD parameters ---------------------------------------------------------
folder_gpd_above94 <- "../luna/kristina/MSGpdParam/"
files <- list.files(folder_gpd_above94)
gpdpar_sites1 <- gpdpar_sites2 <-gpdpar_sites3 <- as.data.frame(matrix(ncol=3,nrow=nrow(xyUK20_sf))) #create empty dataframe
names(gpdpar_sites) <- c("scale","shape","threshold")
# only subset for x in x20 and y in y20
files_subset <- sapply(1:nrow(xyUK20_sf),function(i){paste0("ukgd_cpm85_5k_x",xyUK20_sf$x[i],"y",xyUK20_sf$y[i],".MSGpdParam.2025-02-26.RData")})
#loop through the files
files_subset1 <- files_subset[files_subset %in% files]
# could take only x and y divisible by 4 to subset and speed up data loading
for (i in 1:length(files_subset1)) {
  print(files_subset1[i])
  load(paste0(folder_gpd_above94, files_subset1[i]))
  gpdpar_sites1[i,] <- gpdpar[july3_1976,] #add files to df position
  gpdpar_sites2[i,] <- gpdpar[july3_2026,] #add files to df position
  gpdpar_sites3[i,] <- gpdpar[july3_2076,] #add files to df position
  
}

# find july 3rd in the observed data
july3_obs <- as.numeric(unlist(data_obs_stationary[(92*(1976-1960)+34),]))

# transform margin to original to check
july3_P2q <- unif_orig_P2q(u=july3_obs,P2q=P2q_sites1,gpdpar = gpdpar_sites1)
# compare with observed data
july3_obs1 <- as.numeric(unlist(data_obs[(92*(1976-1960)+34),]))
# plot check
plot(july3_P2q,july3_obs1)

# repeat transform for model functions
july3_2026t <- unif_orig_P2q(u=july3_obs,P2q=P2q_sites2,gpdpar = gpdpar_sites2)
july3_2076t <- unif_orig_P2q(u=july3_obs,P2q=P2q_sites3,gpdpar = gpdpar_sites3)
tmp <- xyUK20_sf %>% mutate(july3_P2q,july3_2026t,july3_2076t)
lims <- c(16,42)
t1 <- tm_shape(tmp) + tm_dots(fill="july3_P2q",size=0.5,fill.scale = tm_scale_intervals(values="viridis",breaks=c(16,20,24,28,32,36,40,44)),fill.legend = tm_legend(title = "Temperature",reverse = TRUE)) + tm_layout(legend.position=c(0.57,0.95),legend.height = 10,frame=FALSE) + tm_title("July 3, 1976") 
t2 <- tm_shape(tmp) + tm_dots(fill="july3_2026t",size=0.5,fill.scale = tm_scale_intervals(values="viridis",breaks=c(16,20,24,28,32,36,40,44)),fill.legend = tm_legend(title = "Temperature",reverse=TRUE)) + tm_layout(legend.show=FALSE,frame=FALSE) + tm_title("July 3, 2026 (projection)") 
t3 <- tm_shape(tmp) + tm_dots(fill="july3_2076t",size=0.5,fill.scale = tm_scale_intervals(values="viridis",breaks=c(16,20,24,28,32,36,40,44)),fill.legend = tm_legend(title = "Temperature",reverse=TRUE)) + tm_layout(legend.show=FALSE,legend.height = 10,frame=FALSE) + tm_title("July 3, 2076 (projection)") 
t <- tmap_arrange(t1,t2,t3,ncol=3)
tmap_save(t,filename=paste0(folder_name,"heatwave1976_future.png"),height=6,width=8)


# 1. plot illustrating conditioning on a variable ------
tmp <- data_obs_stationary_Lap[,c(Birm_index,London_index,Inverness_index)]
names(tmp) <- names(df_sites)[c(1,3,4)]
q <- 0.9
vL <- quantile(tmp[,1],q)
tf <- tmp$tf <- (tmp$London>vL)
limx <-limy <-  c(-10,10)
p1 <- ggplot(tmp) + geom_point(aes(x=London,y=Birmingham),size=0.5) + scale_color_manual(values = c("black", "#009ADA")) + theme(legend.position="none") + coord_fixed() + xlim(limx) + ylim(limy) 
p2 <- ggplot(tmp) + geom_point(aes(x=London,y=Inverness),size=0.5)  + scale_color_manual(values = c("black", "#009ADA")) + theme(legend.position="none") + coord_fixed() + xlim(limx) + ylim(limy)  
p1n <- ggplot(tmp) + geom_point(aes(x=London,y=Birmingham,col=tf),size=0.5) + scale_color_manual(values = c("black", "#009ADA")) + theme(legend.position="none") + coord_fixed() + xlim(limx) + ylim(limy)  + geom_vline(xintercept=vL,color="#009ADA",linetype="dashed")
p2n <- ggplot(tmp) + geom_point(aes(x=London,y=Inverness,col=tf),size=0.5)  + scale_color_manual(values = c("black", "#009ADA")) + theme(legend.position="none") + coord_fixed() + xlim(limx) + ylim(limy)  + geom_vline(xintercept=vL,color="#009ADA",linetype="dashed") 
p <- grid.arrange(p1,p2,ncol=2)
ggsave(p,filename="../Documents/condmodel_illustration.png",width=8,height=4)
pn <- grid.arrange(p1n,p2n,ncol=2)
ggsave(pn,filename="../Documents/condmodel_illustrationvL.png",width=8,height = 4)

# 2. plot of London/Birmingham and London/Inverness on original margins -------
size_point <- 0.3
Birm_index <- find_site_index(site=Birmingham,grid_uk = xyUK20_sf)
London_index <- find_site_index(site=London,grid_uk = xyUK20_sf)
Inverness_index <- find_site_index(site=Inverness,grid_uk = xyUK20_sf)
tmp <- data_obs_all[,c(Birm_index,London_index,Inverness_index)]
names(tmp) <- names(df_sites)[c(1,3,4)]
tmp$tf <- tf
limx <-limy <-  c(7,40)
p1 <- ggplot(tmp) + geom_point(aes(x=London,y=Birmingham),size=size_point) + coord_fixed() + xlim(limx) + ylim(limy) 
p2 <- ggplot(tmp) + geom_point(aes(x=London,y=Inverness),size=size_point) + coord_fixed() + xlim(limx) + ylim(limy) 
p <- grid.arrange(p1,p2,ncol=2)
ggsave(p,filename= "../Documents/London_Birmingham_Inverness_original.png",width=8,height=4)
tmp <- tmp %>% mutate(tf=factor(tf,levels=c(TRUE,FALSE)))
p1 <- ggplot(tmp) + geom_point(aes(x=London,y=Birmingham,col=tf),size=size_point) + scale_color_manual(values = c( "#009ADA","black")) + theme(legend.position="none") + coord_fixed() + xlim(limx) + ylim(limy) 
p2 <- ggplot(tmp) + geom_point(aes(x=London,y=Inverness,col=tf),size=size_point) + scale_color_manual(values = c("#009ADA","black")) + theme(legend.position="none") + coord_fixed() + xlim(limx) + ylim(limy) 
p <- grid.arrange(p1,p2,ncol=2)
ggsave(p,filename= "../Documents/London_Birmingham_Inverness_original_above.png",width=8,height=4)

# load to determine index
load("../luna/kristina/P2q/ukgd_cpm85_5k_x84y20_MSp2q.RData", verbose = TRUE)
# add column for date
data01_obs <- data01 %>% rowid_to_column() %>% filter(class=="obs") %>% mutate("date_obs"=seq(ymd('1960-01-01'),ymd('2023-12-31'),by='days') )
# filter date
july19_2022 <- data01_obs %>% dplyr::filter(date_obs==lubridate::ymd("2022-07-19")) %>% pull(rowid)
# find future dates that match from model data
which.min(abs(data01$time[data01$class=="mod"]- (data01_obs$time[july19_2022] +(2026-2022))))
july19_2026 <- which.min(abs(data01$time[data01$class=="mod"]- (data01_obs$time[july19_2022] +(2026-2022)))) + nrow(data01_obs)
date3 <- which.min(abs(data01$time[data01$class=="mod"]- (data01_obs$time[july19_2022] +(2026-2022+50)))) + nrow(data01_obs)


# load below threshold functions ---------------------------------------------
folder_p2q_below94 <- "../luna/kristina/P2q/"
files <- list.files(folder_p2q_below94)
P2q_sites1 <-P2q_sites2 <-P2q_sites3 <- list() #create empty list
# only subset for x in x20 and y in y20
files_subset <- sapply(1:nrow(xyUK20_sf),function(i){paste0("ukgd_cpm85_5k_x",xyUK20_sf$x[i],"y",xyUK20_sf$y[i],"_MSp2q.RData")})
#loop through the files
files_subset1 <- files_subset[files_subset %in% files]
# could take only x and y divisible by 4 to subset and speed up data loading
for (i in 1:length(files_subset1)) {
  print(files_subset1[i])
  load(paste0(folder_p2q_below94, files_subset1[i]))
  P2q_sites1[[i]] <- qgam.p2q.fn[[july19_2022]] #add files to list position
  P2q_sites2[[i]] <- qgam.p2q.fn[[july19_2026]] #add files to list position
  P2q_sites3[[i]] <- qgam.p2q.fn[[date3]] #add files to list position
  
}

# load GPD parameters ---------------------------------------------------------
folder_gpd_above94 <- "../luna/kristina/MSGpdParam/"
files <- list.files(folder_gpd_above94)
gpdpar_sites1 <- gpdpar_sites2 <-gpdpar_sites3 <- as.data.frame(matrix(ncol=3,nrow=nrow(xyUK20_sf))) #create empty dataframe
names(gpdpar_sites) <- c("scale","shape","threshold")
# only subset for x in x20 and y in y20
files_subset <- sapply(1:nrow(xyUK20_sf),function(i){paste0("ukgd_cpm85_5k_x",xyUK20_sf$x[i],"y",xyUK20_sf$y[i],".MSGpdParam.2025-02-26.RData")})
#loop through the files
files_subset1 <- files_subset[files_subset %in% files]
# could take only x and y divisible by 4 to subset and speed up data loading
for (i in 1:length(files_subset1)) {
  print(files_subset1[i])
  load(paste0(folder_gpd_above94, files_subset1[i]))
  gpdpar_sites1[i,] <- gpdpar[july19_2022,] #add files to df position
  gpdpar_sites2[i,] <- gpdpar[july19_2026,] #add files to df position
  gpdpar_sites3[i,] <- gpdpar[date3,] #add files to df position
  
}

# find the hot day in the observed data
july19_obs <- as.numeric(unlist(data_obs_stationary[(92*(2022-1960)+50),]))
july19_P2q <- unif_orig_P2q(u=july3_obs,P2q=P2q_sites1,gpdpar = gpdpar_sites1)

# repeat transform for model functions
july19_2026t <- unif_orig_P2q(u=july19_obs,P2q=P2q_sites2,gpdpar = gpdpar_sites2)
date3t <- unif_orig_P2q(u=july19_obs,P2q=P2q_sites3,gpdpar = gpdpar_sites3)
tmp <- xyUK20_sf %>% mutate(july19_P2q,july19_2026t,date3t)
lims <- c(16,42)
t1 <- tm_shape(tmp) + tm_dots(fill="july19_P2q",size=0.5,fill.scale = tm_scale_intervals(values="viridis",breaks=c(16,20,24,28,32,36,40,44)),fill.legend = tm_legend(title = "Temperature",reverse = TRUE)) + tm_layout(legend.position=c(0.57,0.95),legend.height = 10,frame=FALSE) + tm_title("July 19, 2022 (TO DO)") 
t2 <- tm_shape(tmp) + tm_dots(fill="july19_2026t",size=0.5,fill.scale = tm_scale_intervals(values="viridis",breaks=c(16,20,24,28,32,36,40,44)),fill.legend = tm_legend(title = "Temperature",reverse=TRUE)) + tm_layout(legend.position=c(0.57,0.95),legend.height = 10,frame=FALSE) + tm_title("July 19, 2026 (projection)") 
t3 <- tm_shape(tmp) + tm_dots(fill="date3t",size=0.5,fill.scale = tm_scale_intervals(values="viridis",breaks=c(16,20,24,28,32,36,40,44)),fill.legend = tm_legend(title = "Temperature",reverse=TRUE)) + tm_layout(legend.position=c(0.57,0.95),legend.height = 10,frame=FALSE) + tm_title("July 19, 2076 (projection)") 
t <- tmap_arrange(t1,t2,t3,ncol=3)
tmap_save(t,filename=paste0(folder_name,"heatwave2022_future.png"),height=6,width=8)

# simulate 10 fields overall --------------------------------------------------
# get index for london
x <- july3_obs[London_index]
# transform to Laplace
xL <- unif_laplace_pit(x)
# get parameter values
load(paste0("data_processed/N9000_sequential2_AGG_all12sites",q*100,".RData"),verbose = TRUE)
aest <- discard(est_all_sf %>% filter(cond_site=="London") %>% pull(a),is.na)
best <- discard(est_all_sf %>% filter(cond_site=="London") %>% pull(b),is.na)
load("data_processed/iterative_phi0l_phi0u_estimates_London.RData",verbose=TRUE)
pe <- as.data.frame(result_new[[12]] %>% dplyr::select(mu_agg,sigl,sigu,deltal,deltau))
names(pe) <- c("mu","sigl","sigu","deltal","deltau")
# get fields
# reconstruct the fields
y_sim <- apply(random10N,MARGIN=c(2),FUN=function(xk){xL*aest+xL^best*xk})
# add row for the conditioning site
y_sim <- as.data.frame(y_sim)
#y_sim <- y_sim %>% add_row(.before=London_index)
y_sim[London_index,] <- xL
names(y_sim) <- paste0("random",1:10)
# plot on standard Normal scale
tmpsf <- st_as_sf(cbind(y_sim,xyUK20_sf)) %>% pivot_longer(cols=contains("random"))
t <- tm_shape(tmpsf %>% mutate("name"=factor(name, levels=unique(tmpsf$name)))) + tm_dots(fill="value",size=0.5,fill.scale = tm_scale_continuous(values="-brewer.rd_bu",limits=c(-10,30)),fill.legend = tm_legend(title = "")) + tm_facets(by="name",ncol=5)
tmap_save(tm=t, filename=paste0(folder_name,"random10_y_sim_London.png"),width=10,height=8)

# transform all 10 onto 2 different scales
#unif_orig_P2q(u=july3_obs,P2q=P2q_sites2,gpdpar = gpdpar_sites2)
x_sim1 <- apply(y_sim,MARGIN=c(2),FUN=function(xk){unif_orig_P2q(u=plaplace(xk),P2q=P2q_sites2,gpdpar = gpdpar_sites2)}) %>% as.data.frame()
names(x_sim1) <- paste0("random",1:10)
# plot on standard Normal scale
tmpsf <- st_as_sf(cbind(x_sim1,xyUK20_sf)) %>% pivot_longer(cols=contains("random"))
t <- tm_shape(tmpsf %>% mutate("name"=factor(name, levels=unique(tmpsf$name)))) + tm_dots(fill="value",size=0.5,fill.scale = tm_scale_intervals(values="-brewer.rd_bu",breaks=seq(10,50,5)),fill.legend = tm_legend(title = "")) + tm_facets(by="name",ncol=5)
tmap_save(tm=t, filename=paste0(folder_name,"random10_x_sim2026.png"),width=10,height=8)
t <- tm_shape(tmpsf %>% mutate("name"=factor(name, levels=unique(tmpsf$name)))) + tm_dots(fill="value",size=0.5,fill.scale = tm_scale_continuous(values="-brewer.rd_bu",limits=c(10,50)),fill.legend = tm_legend(title = "")) + tm_facets(by="name",ncol=5)
tmap_save(tm=t, filename=paste0(folder_name,"random10_xcont_sim2026_London.png"),width=10,height=8)

x_sim2 <- apply(y_sim,MARGIN=c(2),FUN=function(xk){unif_orig_P2q(u=plaplace(xk),P2q=P2q_sites3,gpdpar = gpdpar_sites3)}) %>% as.data.frame()
names(x_sim2) <- paste0("random",1:10)
# plot on standard Normal scale
tmpsf <- st_as_sf(cbind(x_sim2,xyUK20_sf)) %>% pivot_longer(cols=contains("random"))
t <- tm_shape(tmpsf %>% mutate("name"=factor(name, levels=unique(tmpsf$name)))) + tm_dots(fill="value",size=0.5,fill.scale = tm_scale_intervals(values="-brewer.rd_bu",breaks=seq(10,50,5)),fill.legend = tm_legend(title = "")) + tm_facets(by="name",ncol=5)
tmap_save(tm=t, filename=paste0(folder_name,"random10_x_sim2076.png"),width=10,height=8)
t <- tm_shape(tmpsf %>% mutate("name"=factor(name, levels=unique(tmpsf$name)))) + tm_dots(fill="value",size=0.5,fill.scale = tm_scale_continuous(values="-brewer.rd_bu",limits=c(10,50)),fill.legend = tm_legend(title = "")) + tm_facets(by="name",ncol=5)
tmap_save(tm=t, filename=paste0(folder_name,"random10_xcont_sim2076_London.png"),width=10,height=8)

# plot selected 3,4,10
x1sub <- x_sim1[,c(3,4,10)]
x2sub <- x_sim2[,c(3,4,10)]
tmpsf <- st_as_sf(cbind(x1sub,xyUK20_sf)) %>% pivot_longer(cols=contains("random"))
lims <-c(min(min(x1sub),min(x2sub)),max(x1sub,x2sub))
leg_height <- 10
t1 <- tm_shape(tmpsf %>% mutate("name"=factor(name, levels=unique(tmpsf$name))) %>% dplyr::filter(name=="random3")) + tm_dots(fill="value",size=0.5,fill.scale = tm_scale_continuous(values="-brewer.rd_bu",limits=lims),fill.legend = tm_legend(title = "Temperature")) + tm_layout(legend.position=c(0.57,0.95),legend.height = leg_height,frame=FALSE) 
t2 <- tm_shape(tmpsf %>% mutate("name"=factor(name, levels=unique(tmpsf$name))) %>% dplyr::filter(name=="random4")) + tm_dots(fill="value",size=0.5,fill.scale = tm_scale_continuous(values="-brewer.rd_bu",limits=lims),fill.legend = tm_legend(title = "Temperature")) + tm_layout(legend.show=FALSE,legend.height = leg_height,frame=FALSE) 
t3 <- tm_shape(tmpsf %>% mutate("name"=factor(name, levels=unique(tmpsf$name))) %>% dplyr::filter(name=="random10")) + tm_dots(fill="value",size=0.5,fill.scale = tm_scale_continuous(values="-brewer.rd_bu",limits=lims),fill.legend = tm_legend(title = "Temperature")) + tm_layout(legend.show=FALSE,legend.height = leg_height,frame=FALSE) 


tmpsf <- st_as_sf(cbind(x2sub,xyUK20_sf)) %>% pivot_longer(cols=contains("random"))
t4 <- tm_shape(tmpsf %>% mutate("name"=factor(name, levels=unique(tmpsf$name)))%>% filter(name=="random3")) + tm_dots(fill="value",size=0.5,fill.scale = tm_scale_continuous(values="-brewer.rd_bu",limits=lims),fill.legend = tm_legend(title = "Temperature")) + tm_layout(legend.show=FALSE,legend.height = leg_height,frame=FALSE) 
t5 <- tm_shape(tmpsf %>% mutate("name"=factor(name, levels=unique(tmpsf$name)))%>% filter(name=="random4")) + tm_dots(fill="value",size=0.5,fill.scale = tm_scale_continuous(values="-brewer.rd_bu",limits=lims),fill.legend = tm_legend(title = "Temperature")) + tm_layout(legend.show=FALSE,legend.height = leg_height,frame=FALSE) 
t6 <- tm_shape(tmpsf %>% mutate("name"=factor(name, levels=unique(tmpsf$name)))%>% filter(name=="random10")) + tm_dots(fill="value",size=0.5,fill.scale = tm_scale_continuous(values="-brewer.rd_bu",limits=lims),fill.legend = tm_legend(title = "Temperature")) + tm_layout(legend.show=FALSE,legend.height = leg_height,frame=FALSE) 

t <- tmap_arrange(t1,t2,t3,t4,t5,t6,nrow=2)
tmap_save(tm=t, filename=paste0(folder_name,"random3_xcont_sim_London.png"),width=10,height=8)

# 4. AGG density illustration -------------------------------------------------
Lanc_index <- find_site_index(Lancaster,grid_uk = xyUK20_sf)
Z1 <- Z %>% add_column(.before=London_index)
# create a subset of residuals
tmp <- Z1[,c(Lanc_index,Birm_index,Inverness_index)] 
names(tmp) <- c("Birmingham","Lancaster","Inverness")
p <- ggplot(tmp%>% pivot_longer(cols=everything())) + geom_histogram(aes(x=value,fill=name)) + labs(fill="Residual site",x="Residual value",y="Count") + scale_fill_manual(values = c("#C11432","#66A64F", "#009ADA")) 
ggsave(p,filename=paste0(folder_name,"London_residual_site_illustrate.png"),width=6,height=4)
# change the plot to density
# get the parameters
pe3 <- (pe %>% add_row(.before=London_index))[c(Lanc_index,Birm_index,Inverness_index),]
tr1 <- data.frame(y=as.numeric(),x=as.numeric(),Method=as.character())
x <- seq(min(Z1),max(Z1),by=0.1)
for (i in 1:3) {
  sitepar <- as.numeric(pe3[i,])
  tr1 <- rbind(tr1,data.frame(y=AGG_density(x=x,theta = sitepar),x=x,Method=names(tmp)[i]))
}
pl1 <- ggplot() + geom_line(data=tr1 %>% mutate(Method_optim=factor(Method)),aes(x=x,y=y,col=Method)) + labs(col="Residual site",x="Residual value",y="") + scale_fill_manual(values = c("#C11432","#66A64F", "#009ADA")) 
ggsave(pl1,filename=paste0(folder_name,"London_residual_site__density_illustrate.png"),width=6,height=4)


# 5. univariate plot gpd ------------------------------------------------------
xo <- as.numeric(unlist(data_obs_all[,London_index]))
ux <- quantile(xo,0.94)
isabove <- xo>ux
p <- ggplot(data.frame(xo=xo,y=0,tf=isabove)) + geom_density(aes(x=xo)) + geom_point(aes(x=xo,y=y,col=tf),size=0.5) + scale_color_manual(values = c("black", "#009ADA")) + xlab("London summer temperature") + ylab("") + geom_vline(xintercept = ux,linetype="dashed",col="#009ADA") + theme(legend.position="none") +
  annotate("text",x=31,y=0.1,label=">28.34 °C",col="#009ADA")
ggsave(p,filename=paste0(folder_name,"London_illust.png"),width=6,height=4)

# 6. conditional probabilities ------------------------------------------------
# gpd parameters for July 3, 2026
gpd2026 <- as.numeric(unlist(gpdpar_sites2[cond_index,]))
gpd2076 <- as.numeric(unlist(gpdpar_sites3[cond_index,]))
# simulate 900 fields
n_sim <- 900
set.seed(1)
random900 <- as.data.frame(  spam::rmvnorm(n=n_sim,Sigma = Zcov)  )
random900N <- sapply(1:ncol(random900),FUN=function(k) {Normal_AGG_PIT(z = random900[,k],theta=c(pe$mu[k],pe$sigl[k],pe$sigu[k],pe$deltal[1],pe$deltau[1]))}) %>% t %>%  as.data.frame()

# select threshold
u <- 32
p_repeat <- function(u) {
# simulate x
v <- 0.94 + 0.06*evd::pgpd(u,loc = gpd2026[3], scale = gpd2026[1], shape = gpd2026[2])
xL <- qlaplace(v) +  rexp(n = n_sim, rate = 1)
# reconstruct the fields
y_sim <- sapply(1:n_sim,FUN=function(i){xL[i]*aest+xL[i]^best*random900N[,i]})
# add row for the conditioning site
y_sim <- as.data.frame(y_sim)
#y_sim <- y_sim %>% add_row(.before=London_index)
y_sim[London_index,] <- xL
names(y_sim) <- paste0("randomY",1:n_sim)
x_sim1 <- apply(y_sim,MARGIN=c(2),FUN=function(xk){unif_orig_P2q(u=plaplace(xk),P2q=P2q_sites2,gpdpar = gpdpar_sites2)}) %>% as.data.frame()
names(x_sim1) <- paste0("Xstar",1:n_sim)

# count 3 sites exceeding
p1 <- sum(apply(x_sim1,MARGIN=c(2),FUN=function(xk){xk[Lanc_index]>u & xk[Birm_index]>u}))/n_sim
p1i <- sum(apply(x_sim1,MARGIN=c(2),FUN=function(xk){xk[Lanc_index]>u & xk[Birm_index]>u & xk[Inverness_index]>u}))
return(data.frame(u=u,p1=p1,p1i=p1i))
}

tmp <- data.frame(u=numeric(),p1=numeric(),p1i=numeric())
ui <- seq(32,40,1)
for (u in ui) {
 tmp <- rbind(tmp,p_repeat(u))
}
p <- ggplot(tmp) + geom_point(aes(x=u,y=p1)) + labs(x="Temperature at London exceeding",y="Probability of joint exceedance at Birmingham and Lancaster")
ggsave(p,filename=paste0(folder_name,"London_p_exceedance.png"),width=6,height=4)
# repeat for 2076
# simulate x
v <- 0.94 + 0.06*evd::pgpd(u,loc = gpd2076[3], scale = gpd2076[1], shape = gpd2076[2])
xL <- qlaplace(v) +  rexp(n = n_sim, rate = 1)
# reconstruct the fields
y_sim <- sapply(1:n_sim,FUN=function(i){xL[i]*aest+xL[i]^best*random900N[,i]})
# add row for the conditioning site
y_sim <- as.data.frame(y_sim)
#y_sim <- y_sim %>% add_row(.before=London_index)
y_sim[London_index,] <- xL
names(y_sim) <- paste0("randomY",1:n_sim)
x_sim1 <- apply(y_sim,MARGIN=c(2),FUN=function(xk){unif_orig_P2q(u=plaplace(xk),P2q=P2q_sites3,gpdpar = gpdpar_sites3)}) %>% as.data.frame()
names(x_sim1) <- paste0("Xstar",1:n_sim)

# count 3 sites exceeding
p1 <- sum(apply(x_sim1,MARGIN=c(2),FUN=function(xk){xk[Lanc_index]>u & xk[Birm_index]>u}))/n_sim
p1i <- sum(apply(x_sim1,MARGIN=c(2),FUN=function(xk){xk[Lanc_index]>u & xk[Birm_index]>u & xk[Inverness_index]>u}))
    
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
# add column for date
data01_obs <- data01 %>% rowid_to_column() %>% filter(class=="obs") %>% mutate("date_obs"=seq(ymd('1960-01-01'),ymd('2023-12-31'),by='days') )
# filter date
july3_1976 <- data01_obs %>% dplyr::filter(date_obs==lubridate::ymd("1976-07-03")) %>% pull(rowid)
# find future dates that match from model data
which.min(abs(data01$time[data01$class=="mod"]- (data01_obs$time[july3_1976] +50)))
july3_2026 <- which.min(abs(data01$time[data01$class=="mod"]- (data01_obs$time[july3_1976] +50))) + nrow(data01_obs)
july3_2076 <- which.min(abs(data01$time[data01$class=="mod"]- (data01_obs$time[july3_1976] +100))) + nrow(data01_obs)

# compare P2q of the observed 22846 and model 38362 data
load("../luna/kristina/MSGpdParam/ukgd_cpm85_5k_x84y20.MSGpdParam.2025-02-26.RData", verbose = TRUE)
gpdpar[c(22846,38362),]
qgam.p2q.fn[[22845]](data01$u[22845])
qgam.p2q.fn[[38362]](data01$u[22846])
data01[c(22846,38362),]


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
t2 <- tm_shape(tmp) + tm_dots(fill="july3_2026t",size=0.5,fill.scale = tm_scale_intervals(values="viridis",breaks=c(16,20,24,28,32,36,40,44)),fill.legend = tm_legend(title = "Temperature",reverse=TRUE)) + tm_layout(legend.position=c(0.57,0.95),legend.height = 10,frame=FALSE) + tm_title("July 3, 2026 (projection)") 
t3 <- tm_shape(tmp) + tm_dots(fill="july3_2076t",size=0.5,fill.scale = tm_scale_intervals(values="viridis",breaks=c(16,20,24,28,32,36,40,44)),fill.legend = tm_legend(title = "Temperature",reverse=TRUE)) + tm_layout(legend.position=c(0.57,0.95),legend.height = 10,frame=FALSE) + tm_title("July 3, 2076 (projection)") 
t <- tmap_arrange(t1,t2,t3,ncol=3)
tmap_save(t,filename=paste0(folder_name,"heatwave1976_future.png"),height=6,width=8)

# 1. plot of London/Birmingham and London/Inverness on original margins -------
Birm_index <- find_site_index(site=Birmingham,grid_uk = xyUK20_sf)
London_index <- find_site_index(site=London,grid_uk = xyUK20_sf)
Inverness_index <- find_site_index(site=Inverness,grid_uk = xyUK20_sf)
tmp <- data_obs[,c(Birm_index,London_index,Inverness_index)]
names(tmp) <- names(df_sites)[c(1,3,4)]
p1 <- ggplot(tmp) + geom_point(aes(x=London,y=Birmingham),size=0.5)
p2 <- ggplot(tmp) + geom_point(aes(x=London,y=Inverness),size=0.5)
p <- grid.arrange(p1,p2,ncol=2)
ggsave(p,filename= "../Documents/London_Birmingham_Inverness_original.png",width=8,height=4)

# 2. plot illustrating conditioning on a variable ------
tmp <- data_obs[,c(Birm_index,London_index,Inverness_index)]
names(tmp) <- names(df_sites)[c(1,3,4)]
q <- 0.9
vL <- quantile(tmp[,1],q)
tmp$tf <- (tmp$London>vL)
limx <-limy <-  c(7,36)
p1 <- ggplot(tmp) + geom_point(aes(x=London,y=Birmingham),size=0.5) + xlab(TeX("$Y_1$")) + ylab(TeX("$Y_2$")) + scale_color_manual(values = c("black", "#009ADA")) + theme(legend.position="none") + coord_fixed() + xlim(limx) + ylim(limy) 
p2 <- ggplot(tmp) + geom_point(aes(x=London,y=Inverness),size=0.5) + xlab(TeX("$Y_1$")) + ylab(TeX("$Y_3$")) + scale_color_manual(values = c("black", "#009ADA")) + theme(legend.position="none") + coord_fixed() + xlim(limx) + ylim(limy)  
p1n <- ggplot(tmp) + geom_point(aes(x=London,y=Birmingham,col=tf),size=0.5) + scale_color_manual(values = c("black", "#009ADA")) + theme(legend.position="none") + coord_fixed() + xlim(limx) + ylim(limy)  + geom_vline(xintercept=vL,color="#009ADA",linetype="dashed")
p2n <- ggplot(tmp) + geom_point(aes(x=London,y=Inverness,col=tf),size=0.5)  + scale_color_manual(values = c("black", "#009ADA")) + theme(legend.position="none") + coord_fixed() + xlim(limx) + ylim(limy)  + geom_vline(xintercept=vL,color="#009ADA",linetype="dashed") 
p <- grid.arrange(p1,p2,ncol=2)
ggsave(p,filename="../Documents/condmodel_illustration.png",width=8,height=4)
pn <- grid.arrange(p1n,p2n,ncol=2)
ggsave(pn,filename="../Documents/condmodel_illustrationvL.png",width=8,height = 4)


library(tidyverse)
library(tmap)
library(lubridate)
theme_set(theme_bw())
theme_replace(
  panel.spacing = unit(2, "lines"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA) )

# load observed data
#source("spatial_parameter_estimation.R") # for spatial_par_est function
load("data_processed/temperature_data.RData",verbose = TRUE)
load("data_processed/spatial_helper.RData", verbose = TRUE)

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
list_of_files <- list() #create empty list
# only subset for x in x20 and y in y20
files_subset <- sapply(1:nrow(xyUK20_sf),function(i){paste0("ukgd_cpm85_5k_x",xyUK20_sf$x[i],"y",xyUK20_sf$y[i],"_MSp2q.RData")})
#loop through the files
files_subset1 <- files_subset[files_subset %in% files]
# plot the missing files
t <- tm_shape(xyUK20_sf) + tm_dots() + tm_shape(xyUK20_sf %>% filter(x %in% 96,!(y %in% c(120,124,128,132,136,140,144)))) + tm_dots("red")
tmap_save(t,filename = "../Documents/missing_p2q_files.png")
# could take only x and y divisible by 4 to subset and speed up data loading
for (i in 1:length(files_subset1)) {
  print(files_subset1[i])
  load(paste0(folder_p2q_below94, files_subset1[i]))
  list_of_files[[i]] <- qgam.p2q.fn[[july3_1976]] #add files to list position
}

# load GPD parameters ---------------------------------------------------------
folder_p2q_below94 <- "../luna/kristina/P2q/"
files <- list.files(folder_p2q_below94)
list_of_files <- list() #create empty dataframe
# only subset for x in x20 and y in y20
files_subset <- sapply(1:nrow(xyUK20_sf),function(i){paste0("ukgd_cpm85_5k_x",xyUK20_sf$x[i],"y",xyUK20_sf$y[i],"_MSp2q.RData")})
#loop through the files
files_subset1 <- files_subset[files_subset %in% files]
# plot the missing files
t <- tm_shape(xyUK20_sf) + tm_dots() + tm_shape(xyUK20_sf %>% filter(x %in% 96,!(y %in% c(120,124,128,132,136,140,144)))) + tm_dots("red")
tmap_save(t,filename = "../Documents/missing_p2q_files.png")
# could take only x and y divisible by 4 to subset and speed up data loading
for (i in 1:length(files_subset1)) {
  print(files_subset1[i])
  load(paste0(folder_p2q_below94, files_subset1[i]))
  list_of_files[[i]] <- qgam.p2q.fn[[july3_1976]] #add files to list position
}


# find july 3rd in the observed data
july3_obs <- as.numeric(unlist(data_obs[(92*(1976-1960)+34),]))

# join the summer data together could try 1960-1999 for non-stationary data -----
# June 1 is 152 doy, August 31 is 243 doy (92 days per year)
# create a dataframe
xyUK20 <- xyUK20_sf %>% dplyr::select(-temp) %>% st_drop_geometry() %>% cbind(as.data.frame(matrix(data=numeric(),ncol=92*40,nrow=nrow(xyUK20_sf))))
names(xyUK20)[5:ncol(xyUK20)] <- paste0(rep(152:243,40),"_",rep(1960:1999,each=92)) 
for (i in 1: length(list_of_files)) {
  xyUK20[i,5:ncol(xyUK20)] <- list_of_files[[i]] %>% mutate(year=floor(time)) %>% filter(class=="obs",doy>=152,doy<=243, year<=1999) %>% pull(x) 
}

# explore issue with deltas Birmingham to Cromer
summary(est_all_sf)
est_all_sf %>% 
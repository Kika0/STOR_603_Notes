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
x=seq.Date(from=as.Date("1960-01-01"), to=as.Date("2023-12-31"), by="day")
data01_obs <- data01 %>% rowid_to_column() %>% filter(class=="obs") %>% mutate("date"=as.Date(x))
# filter date
july3_1976 <- data01_obs %>% filter("date"==as.Date("1976-07-19")) %>% pull(rowid)
data01 %>% rowid_to_column() %>% filter(class=="mod",time<2022.55,time>2022.54)
# compare P2q of the observed 22846 and model 38362 data
load("../luna/kristina/MSGpdParam/ukgd_cpm85_5k_x84y20.MSGpdParam.2025-02-26.RData", verbose = TRUE)
gpdpar[c(22846,38362),]
qgam.p2q.fn[[22845]](data01$u[22845])
qgam.p2q.fn[[38362]](data01$u[22846])
data01[c(22846,38362),]



#
files <- list.files("../luna/kristina/MSdata01/")
list_of_files <- list() #create empty list
# only subset for x in x20 and y in y20
files_subset <- sapply(1:nrow(xyUK20_sf),function(i){paste0("ukgd_cpm85_5k_x",xyUK20_sf$x[i],"y",xyUK20_sf$y[i],"_MSdata01.RData")})
#loop through the files
files_subset1 <- files_subset[files_subset %in% files]
# could take only x and y divisible by 4 to subset and speed up data loading
for (i in 1:length(files_subset1)) {
  print(files_subset1[i])
  load(paste0("../luna/kristina/MSdata01/", files_subset1[i]))
  list_of_files[[i]] <- data01 #add files to list position
}
xyUK20_sf <- xyUK20_sf[files_subset %in% files_subset1,]

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
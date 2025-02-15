library(tidyverse)
library(ncdf4)
library(fields)
# load the data
load("../luna/kristina/MSdata_CPM5km_member1/ukgd_cpm85_5k_x100y23_MSdata01.RData")
data01.std.param # not much info here
data01 <- data01 %>% mutate(year=floor(time))
dim(data01)
summary(data01)
glimpse(data01)
data01 %>% group_by(year) %>% summarize(mean_year=mean(x)) %>% ggplot(aes(x=year,y=mean_year)) + geom_point()
ggplot(data01 %>% filter(doy %in% c(1:20))) + geom_point(aes(x=time,y=x))
# to explore gpd functions for the tails
str(gpdpar)
summary(chosen.MSgpd)

# load the Gpd parameters
load("../luna/kristina/MSGpdParam/ukgd_cpm85_5k_x100y23.MSGpdParam.2024-10-22-064747.RData")
glimpse(gpdpar) # dataframe for columns for scale, shape and threshold

###############################################################################
# 5km
#
# NOTE the obs are on a bigger grid than the CPM grid
###############################################################################
obs_cpm_offset_y <- 33  # the difference in size of the y dimension

xcoord_m         <- 146 # London x146 y44 on the CPM grid from Laura's city file.
ycoord_m         <- 44
xcoord_o         <- xcoord_m
ycoord_o         <- ycoord_m + obs_cpm_offset_y

### OBS 5km
obs_example  <- '../luna/kristina/UKgrid5km/tasmax_rcp85_land-cpm_uk_5km_01_ann_206012-208011.nc'
nc1      <- nc_open(obs_example)
vlist    <- nc1$var
shape.o  <- vlist$tasmax$size
tas5.o   <- ncvar_get(nc1, "tasmax")
lon5.o   <- ncvar_get(nc1, "longitude")
lat5.o   <- ncvar_get(nc1, "latitude")
nc_close(nc1)

### CPM 5km
cpm2k_example <- "//tasmax_rcp85_land-cpm_uk_5km_01_day_19901201-20001130.nc"
nc1      <- nc_open(cpm2k_example)
vlist    <- nc1$var
shape.m  <- vlist$tasmax$size
tas5.m   <- ncvar_get(nc1, "tasmax")
lon5.m   <- ncvar_get(nc1, "longitude")
lat5.m   <- ncvar_get(nc1, "latitude")
nc_close(nc1)

par(mfcol=c(2,2))

quilt.plot(as.vector(lon5.o), as.vector(lat5.o), tas5.o[,,10],  nx=shape.o[1], ny=shape.o[2], main='OBS 5km')
quilt.plot(as.vector(lon5.m), as.vector(lat5.m), tas5.m[,,180], nx=shape.m[1], ny=shape.m[2], main='CPM 5km')

image.plot(tas5.o[,,10],  main='OBS 5km')
image.plot(tas5.m[,,180], main='CPM 5km')

#  str(tas5.o) num [1:180, 1:290, 1:31] NA NA NA NA NA NA NA NA NA NA ...
#  str(lon5.o) num [1:180, 1:290] -9.99 -9.92 -9.86 -9.79 -9.72 ...
#  str(lat5.o) num [1:180, 1:290] 47.8 47.8 47.9 47.9 47.9 ...
#
#  str(tas5.m) num [1:180, 1:244, 1:3600] 10.1 10.2 10.3 10.3 10.2 ...
#  str(lon5.m) num [1:180, 1:244] -10.23 -10.16 -10.09 -10.02 -9.95 ...
#  str(lat5.m) num [1:180, 1:244] 49.3 49.3 49.3 49.3 49.3 ...

save(file="grid-info-5km.RData", lon5.o, lat5.o, lon5.m, lat5.m, obs_example, cpm2k_example)

# load all non-stationary data
file.sources = list.files(pattern="*helpers.R")
sapply(file.sources,source,.GlobalEnv)

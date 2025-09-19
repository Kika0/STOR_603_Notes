library(tidyverse)
library(ncdf4)
library(fields)
library(sf)
library(tmap)
load("data_processed/spatial_helper.RData")
# load the data
load("../luna/kristina/MSdata01/ukgd_cpm85_5k_x100y24_MSdata01.RData")
data01 <- data01 %>% mutate(year=floor(time))

#data01 %>% group_by(year) %>% summarize(mean_year=mean(x)) %>% ggplot(aes(x=year,y=mean_year)) + geom_point()

# load the Gpd parameters
#load("../luna/kristina/scratch/hadsx/heatwave/HadUKGrid/dur-clim/CPM5km/v2kristina/UK/01/MSGpdParam/ukgd_cpm85_5k_x76y148.MSGpdParam.2025-02-26.RData")
#glimpse(gpdpar) # dataframe for columns for scale, shape and threshold

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
obs_example  <- '../luna/kristina_old/UKgrid5km/tasmax_rcp85_land-cpm_uk_5km_01_ann_206012-208011.nc'
nc1      <- nc_open(obs_example)
vlist    <- nc1$var
shape.o  <- vlist$tasmax$size
tas5.o   <- ncvar_get(nc1, "tasmax")
lon5.o   <- ncvar_get(nc1, "longitude")
lat5.o   <- ncvar_get(nc1, "latitude")
nc_close(nc1)
par(mfcol=c(1,2))

quilt.plot(as.vector(lon5.o), as.vector(lat5.o), tas5.o[,,10],  nx=shape.o[1], ny=shape.o[2], main='OBS 5km')

image.plot(tas5.o[,,10],  main='OBS 5km')

#  str(tas5.o) num [1:180, 1:290, 1:31] NA NA NA NA NA NA NA NA NA NA ...
#  str(lon5.o) num [1:180, 1:290] -9.99 -9.92 -9.86 -9.79 -9.72 ...
#  str(lat5.o) num [1:180, 1:290] 47.8 47.8 47.9 47.9 47.9 ...
#
#  str(tas5.m) num [1:180, 1:244, 1:3600] 10.1 10.2 10.3 10.3 10.2 ...
#  str(lon5.m) num [1:180, 1:244] -10.23 -10.16 -10.09 -10.02 -9.95 ...
#  str(lat5.m) num [1:180, 1:244] 49.3 49.3 49.3 49.3 49.3 ...

#save(file="grid-info-5km.RData", lon5.o, lat5.o, lon5.m, lat5.m, obs_example, cpm2k_example)

# NOTE: only need to work with files that overlap with mainland UK
# find a subset of x and y that overlap with mainland UK
summary(as.vector(lon5.o))
summary(as.vector(lat5.o))
# create a dataframe with also x (row) and y (column) indeces
# as.vector does column by column
xy_df <- data.frame("lon"=as.vector(lon5.o),"lat"=as.vector(lat5.o),"x"= rep(1:dim(lon5.o)[1],n=dim(lon5.o)[2]), "y" = rep(1:dim(lon5.o)[2],each=dim(lon5.o)[1]), "temp"=as.vector(tas5.o[,,10]))
xy_sf <- xy_df %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POINT") %>% st_make_valid()
xy_sf <- cbind(xy_sf,xy_df)
# map
tmap_mode("view")
tm_shape(xy_sf) + tm_dots(col="temp")
# looks resonable, now subset over mainland UK
lad <- st_read("data/LAD_boundary_UK/LAD_MAY_2022_UK_BFE_V3.shp")
# UK is comprised of many polygons (islands), simplify to only
# take mainland UK (Great Britain)
uk_notsimplified <- (lad %>% st_union() %>% st_cast( "MULTIPOLYGON" ) %>% st_cast("POLYGON"))[1531]
uk <- st_simplify(uk_notsimplified,dTolerance = 2000) %>% st_transform(crs = 4326)  
# check
# tm_shape(uk_notsimplified) + tm_polygons()
# great, now subset
xyUK_sf <- st_filter(xy_sf,uk)
tm_shape(xyUK_sf) + tm_dots(col="temp")
# take only every fourth dot in each x and y
x20 <- seq(from=4,to=dim(lon5.o)[1],by=4)
y20 <- seq(from=4,to=dim(lon5.o)[2],by=4)
xyUK20_sf <-xyUK_sf %>% filter(x %in% x20,y %in% y20)
# save sf objects used for spatial analysis
save(uk,uk_notsimplified,xyUK20_sf,files_subset1,file="data_processed/spatial_helper.RData")


# start here to get data---
tm_shape(xyUK20_sf) + tm_dots(col="temp")
# great, now load only files that overlap this grid or perhaps delete all other files?
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
# setup for par_est
data_obs <- xyUK20 %>% dplyr::select(-all_of(1:4)) %>% t() %>% as.data.frame()
# ordered alphabetically so Y1 Birmingham, Y2 Glasgow and Y3 is London
colnames(data_obs) <- paste0("Y",1:ncol(data_obs))
# transform to Laplace margins
data_obs_Lap <- as.data.frame((data_obs %>% apply(c(2),FUN=row_number))/(nrow(data_obs)+1)) %>% apply(c(1,2),FUN=unif_laplace_pit) %>% as.data.frame()

# create a dataframe for stationary data ----
xyUK20 <- xyUK20_sf %>% dplyr::select(-temp) %>% st_drop_geometry() %>% cbind(as.data.frame(matrix(data=numeric(),ncol=90*100,nrow=nrow(xyUK20_sf))))
names(xyUK20)[5:ncol(xyUK20)] <- paste0(rep(91:180,40),"_",rep(1981:2080,each=90)) 

for (i in 1: length(list_of_files)) {
  xyUK20[i,5:ncol(xyUK20)] <- list_of_files[[i]] %>% mutate(year=floor(time)) %>% filter(class=="mod") %>% pull(x) 
}
# setup for par_est
data_mod <- xyUK20 %>% dplyr::select(-all_of(1:4)) %>% t() %>% as.data.frame()
# ordered alphabetically so Y1 Birmingham, Y2 Glasgow and Y3 is London
colnames(data_mod) <- paste0("Y",1:ncol(data_mod))
# transform to Laplace margins
data_mod_Lap <- as.data.frame((data_mod %>% apply(c(2),FUN=row_number))/(nrow(data_mod)+1)) %>% apply(c(1,2),FUN=unif_laplace_pit) %>% as.data.frame()

# save as R object for further analysis
save(data_obs,data_obs_Lap,data_mod,data_mod_Lap,file = "data_processed/temperature_data.RData")

# explore other datasets
files <- list.files("../luna/kristina/MSGpdParam_CPM5km_member1_20x20/scratch/hadsx/heatwave/HadUKGrid/dur-clim/CPM5km/v2kristina/UK/01/MakeStationary/MSGpdParam/")
list_of_files <- list() #create empty list
# only subset for x in x20 and y in y20
files_subset <- sapply(1:nrow(xyUK20_sf),function(i){paste0("ukgd_cpm85_5k_x",xyUK20_sf$x[i],"y",xyUK20_sf$y[i],".MSGpdParam.2025-02-26.RData")})
#loop through the files
files_subset1 <- files_subset[files_subset %in% files]
# could take only x and y divisible by 4 to subset and speed up data loading
for (i in 1:length(files_subset1)) {
  print(files_subset1[i])
  load(paste0("../luna/kristina/MSGpdParam_CPM5km_member1_20x20/scratch/hadsx/heatwave/HadUKGrid/dur-clim/CPM5km/v2kristina/UK/01/MakeStationary/MSGpdParam/", files_subset1[i]))
  list_of_files[[i]] <- gpdpar #add files to list position
}

i    <- which(trunc(data01$time)==2020 & data01$doy==205 & data01$class=='obs')
k    <- which(trunc(data01$time)==2079 & data01$doy==205 & data01$class=='mod')

glimpse(list_of_files[[1]]$threshold[i])
# setup for par_est
thresCPM <- thresobs <- c()
for (j in 1: length(list_of_files)) {
thresCPM[j] <-   list_of_files[[j]]$threshold[k]
thresobs[j] <-   list_of_files[[j]]$threshold[i]
}
tmap_mode("plot")
tm_thres <- xyUK20_sf %>% mutate("CPM" =thresCPM, "observed"=thresobs) %>% pivot_longer(cols=c("CPM","observed"),values_to = "temperature", names_to = "data_source")
t <- tm_shape(tm_thres) + tm_dots(fill="temperature",size=0.8,fill.scale =tm_scale_continuous(values="-matplotlib.rd_yl_bu")) + tm_facets(by = c("data_source")) + tm_layout(legend.position=c("right","top"),legend.height = 12)
# save map
tmap_save(t,filename=paste0("../Documents/threshold_explore.png"),width=8,height=6)

# add this to object of analysis data to be potentially used as a covariate

library(tidyverse)
library(ncdf4)
library(fields)
library(sf)
library(tmap)
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
gridlon5.o
summary(as.vector(lon5.o))
summary(as.vector(lat5.o))
# create a dataframe with also x (row) and y (column) indeces
# as.vector does column by column
xy_df <- data.frame("lon"=as.vector(lon5.o),"lat"=as.vector(lat5.o),"x"= rep(1:dim(lon5.o)[1],n=dim(lon5.o)[2]), "y" = rep(1:dim(lon5.o)[2],each=dim(lon5.o)[1]))
xy_sf <- xy_df %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POINT") %>% st_make_valid()
xy_sf <- cbind(xy_sf,xy_df)
# map
tmap_mode("view")
tm_shape(xy_sf) + tm_dots()
# looks resonable, now subset over mainland UK
lad <- st_read("data/LAD_boundary_UK/LAD_MAY_2022_UK_BFE_V3.shp")
# UK is comprised of many polygons (islands), simplify to only
# take mainland UK (Great Britain)
uk <- (lad %>% st_union() %>% st_cast( "MULTIPOLYGON" ) %>% st_cast("POLYGON"))[1531]
uk <- st_simplify(uk,dTolerance = 2000) %>% st_transform(crs = 4326)  
# check
# tm_shape(uk) + tm_polygons()
# great, now subset
xyUK_sf <- st_filter(xy_sf,uk)
tm_shape(xyUK_sf) + tm_dots()
# take only every fourth dot in each x and y
x20 <- seq(from=4,to=dim(lon5.o)[1],by=4)
y20 <- seq(from=4,to=dim(lon5.o)[2],by=4)
xyUK20_sf <-xyUK_sf %>% filter(x %in% x20,y %in% y20)
tm_shape(xyUK20_sf) + tm_dots()
# great, now load only files that overlap this grid or perhaps delete all other files?
files <- list.files("../luna/kristina/MSdata_CPM5km_member1/")
list_of_files <- list() #create empty list
# only subset for x in x20 and y in y20
files_subset <- sapply(1:nrow(xyUK20_sf),function(i){paste0("ukgd_cpm85_5k_x",xyUK20_sf$x[i],"y",xyUK20_sf$y[i],"_MSdata01.RData")})
#loop through the files
files_subset1 <- files[files %in% files_subset]
# could take only x and y divisible by 4 to subset and speed up data loading
for (i in 1:length(files_subset1)) {
  print(files_subset1[i])
  load(paste0("../luna/kristina/MSdata_CPM5km_member1/", files_subset1[i]))
  list_of_files[[i]] <- data01 #add files to list position
}
xyUK20_sf <- xyUK20_sf[files_subset %in% files_subset1,]
tm_shape(xyUK20_sf) + tm_dots()

# join the summer data together could try 1960-1999 for non-stationary data 
# June 1 is 152 doy, August 31 is 243 doy (92 days per year)
# create a dataframe
xyUK20 <- xyUK20_sf %>% st_drop_geometry() %>% cbind(as.data.frame(matrix(data=numeric(),ncol=92*40,nrow=nrow(xyUK20_sf))))
names(xyUK20)[5:ncol(xyUK20)] <- paste0(rep(152:243,40),"_",rep(1960:1999,each=92)) 

for (i in 1: length(list_of_files)) {
  xyUK20[i,5:ncol(xyUK20)] <- list_of_files[[i]] %>% mutate(year=floor(time)) %>% filter(class=="obs",doy>=152,doy<=243, year<=1999) %>% pull(x) 

}

# save as R object for further analysis

# explore NA values
# p <- data01 %>% mutate(year=floor(time)) %>% filter(is.na(uqgam)) %>% group_by(year,class) %>% summarize(n=n()) %>% arrange(-n) %>% ggplot() + geom_line(aes(x=year,y=n,col=class)) + ggtitle("Counts of NA values of uqgam for each year for observation and model data")
# ggsave(p,filename="../Documents/newdata/uqgamNA.png",width=10,height=5)

# add data

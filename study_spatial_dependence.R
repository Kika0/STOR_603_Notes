library(ncdf4)
library(fields)
library(maps)
library(tmap)
library(mapdata)
library(PCICt)
library(rgdal)
library(viridis) # for palettes
library(tidyverse)
library(sf)
library(units)
#library(gridExtra)
source("rotate_unrotate_coordinates.R")

### subset to show only UK mainland overlapping grid points ----
lad <- st_read("data/LAD_boundary_UK/LAD_MAY_2022_UK_BFE_V3.shp")
# UK is comprised of many polygons (islands), simplify to only
# take mainland UK (Great Britain)
uk <- (lad %>% st_union() %>% st_cast( "MULTIPOLYGON" ) %>% st_cast("POLYGON"))[1019]
uk <- st_simplify(uk,dTolerance = 7000)
# rotate coordinates of uk polygon
long <- st_coordinates(uk %>% st_transform(4326))[,1]
lat <- st_coordinates(uk %>% st_transform(4326))[,2]
# loop to replace coordinates with rotated
lat_rot <- c()
long_rot <- c()
for (i in 1:length(long)){
  r.latlon  <- pp.ll.to.rg(lat=lat[i],long=long[i], gr_npole_lat, gr_npole_lon)
  lat_rot[i] <- r.latlon[1]
  long_rot[i] <- r.latlon[2]
}

df <- data.frame("lon"=long_rot,"lat"=lat_rot)
# convert back to sf polygon
uk_rot <-  df %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON") %>% st_make_valid()
rm(df)

### read one year of data ----
flist  <- "data/tasmax_rcp85_land-cpm_uk_2.2km_01_day_19991201-20001130.nc"
### read netcdf file
nc1 <- nc_open(flist[1])
# time coordinate
time1    <- ncvar_get(nc1,'time')
time1att <- ncatt_get(nc1,'time')
# read all 360 days
# time dimension is last in the start argument
v1    <- ncvar_get(nc1,'tasmax',start=c(1,1,1,1), count=c(-1,-1,360,1))
# lat long of grid
rlon     <- apply(ncvar_get(nc1,"grid_longitude_bnds"),2,mean)
rlat     <- apply(ncvar_get(nc1,"grid_latitude_bnds"),2,mean)
cpm.pole <- ncatt_get(nc1,'rotated_latitude_longitude')
nc_close(nc1)
# this climate model has a 360 day year so needs special treatment
# PCICt does this for us
origin   <- tail(unlist(strsplit(time1att$units,'hours since ')),1)
pcictime <- as.PCICt( time1*(60*60), cal=time1att$calendar ,origin=origin)

# the CPM grid has the pole in a different place to usual.  therefore need to transform the map coordinates
gr_npole_lat <- cpm.pole$grid_north_pole_latitude
gr_npole_lon <- cpm.pole$grid_north_pole_longitude


### plot a single field -----
longitude <- rlon
latitude <- rlat
### subset the grid
lon_subset <- seq(from=10,to=(length(longitude)-length(longitude)%%10),by=10)
lat_subset <- seq(from=10,to=(length(latitude)-length(latitude)%%10),by=10)
longitude <- longitude[lon_subset]
latitude <- latitude[lat_subset]
v1_sub <- v1[lon_subset,lat_subset,]
#fields::image.plot( x=longitude, y=latitude, z=v1_sub[,,1], zlim=c(0,20), main=pcictime[2])

### subset rotated grid points overlapping with UK -----
# create grid points for long lat
lon_lat <- expand.grid(longitude,latitude) 
colnames(lon_lat) <- c("lon","lat")
# extract first day (first out of 360 in dim3:time)
Temperature <- c(v1_sub[,,1])
lon_lat_temp <- cbind(lon_lat,Temperature)
for (i in 2:(dim(v1_sub)[3]/1)) {
  lon_lat_temp <- cbind(lon_lat_temp,c(v1_sub[,,i]))
  colnames(lon_lat_temp)[length(lon_lat_temp)] <- as.character(i)
}
# convert to sf points object
temp_sf <- st_as_sf(lon_lat_temp,coords = c("lon","lat"),crs=4326)
tm_shape(temp_sf %>% select(Temperature)) + tm_dots(col="Temperature",style="cont",palette="viridis")
tm_shape(temp_sf %>% select(Temperature)) + tm_dots(col="Temperature",style="cont",size=0.05,palette="-RdYlBu",legend.col.reverse=TRUE)+
  tm_layout(legend.position = c(0.7, 0.6)) + tm_grid(alpha=0.3) + tm_xlab("Longitude") + tm_ylab("Latitude")
# subset to only include grid points in mainland UK
uk_temp_sf <- st_filter(temp_sf,(uk_rot %>% st_cast("MULTIPOLYGON")))
tm_shape(uk_temp_sf) + tm_dots(col="Temperature",style="cont",size=0.01,palette="viridis")

# try to plot for 30 days ----
# subset to only include grid points in Lancaster
uk_temp_sf_long <- uk_temp_sf %>% select(Temperature) %>%
  mutate("day"=as.factor(rep(1,nrow(uk_temp_sf))))
#  mutate("date"=rep(pcictime[1],nrow(uk_temp_sf)))
# December subset 
#pcictime[1:30] # last day is missing

# for loop for all other half-hour intervals
for (j in (2:(dim(v1_sub)[3]/1))) {
  to_bind <- uk_temp_sf  %>% select(all_of(j)) %>% 
    mutate("day"=as.factor(rep(j,nrow(uk_temp_sf)))) 
  #  mutate("date"=rep(pcictime[j],nrow(uk_temp_sf)))
  names(to_bind)[1] <- "Temperature"
  uk_temp_sf_long <- rbind(uk_temp_sf_long,to_bind)
  
}

tm_shape(uk_temp_sf_long) + 
  tm_dots(col="Temperature",style="cont",size=0.05,palette="viridis") +
  tm_facets(by="day",as.layers=TRUE,ncol=30,nrow=12) +
  tm_layout(panel.show = FALSE,between.margin = 2)
# facet by timestamp
# uk_temp_sf_long$date <- factor(uk_temp_sf_long$date,      # Reordering group factor levels
#                                  levels = (uk_temp_sf_long$date %>% unique()))
# tm_shape(uk_temp_sf_long) + 
#   tm_dots(col="Temperature",style="cont",size=0.9,palette="viridis",title = "Temperature") +
#   tm_facets(by="date",as.layers=TRUE,ncol=30,nrow=12)  + 
#   theme(
#     strip.background = element_blank(),
#     strip.text.x = element_blank()
#   )

### examine 1999 data ----
uk_temp_sf_long$Temperature %>% summary()

max_1999 <- uk_temp_sf_long[uk_temp_sf_long$Temperature==max(uk_temp_sf_long$Temperature),]
long <- st_coordinates(max_1999 %>% st_transform(4326))[,1]
lat <- st_coordinates(max_1999 %>% st_transform(4326))[,2]

conv <- CnvRttPol(latlon = data.frame(long,lat),spol_coor = c(gr_npole_lon, gr_npole_lat))
max_point <- data.frame(lon=conv$lon,lat=conv$lat) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  summarise(geometry = st_combine(geometry))

tmap_mode("view")
tm_shape(max_point) + tm_dots()

min_1999 <- uk_temp_sf_long[uk_temp_sf_long$Temperature==min(uk_temp_sf_long$Temperature),]
long <- st_coordinates(min_1999 %>% st_transform(4326))[,1]
lat <- st_coordinates(min_1999 %>% st_transform(4326))[,2]

conv <- CnvRttPol(latlon = data.frame(long,lat),spol_coor = c(gr_npole_lon, gr_npole_lat))
min_point <- data.frame(lon=conv$lon,lat=conv$lat) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  summarise(geometry = st_combine(geometry))

tmap_mode("view")
tm_shape(min_point) + tm_dots()
tmap_mode("plot")

# try unrotating all
long <- st_coordinates(uk_temp_sf %>% st_transform(4326))[,1]
lat <- st_coordinates(uk_temp_sf %>% st_transform(4326))[,2]

conv <- CnvRttPol(latlon = data.frame(long,lat),spol_coor = c(gr_npole_lon, gr_npole_lat))
uk_unrotated <- data.frame(lon=conv$lon,lat=conv$lat) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  summarise(geometry = st_combine(geometry))

tmap_mode("view")
tm_shape(uk_unrotated) + tm_dots()
tmap_mode("plot")


### create London dataset ----
# coord for London centre (Alan Turing Institute)
london_lat <- 51.529972
london_lon <- -0.127676
# plot to check
# london_sf <- data.frame(lon=london_lon,lat=london_lat) %>%
#   st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
#   summarise(geometry = st_combine(geometry))
# tmap_mode("view")
# tm_shape(london_sf) + tm_dots()
# tmap_mode("plot")
# rotate coordinates
london_lon_rot <- pp.ll.to.rg(lat=london_lat,long = london_lon,pole.lat =  gr_npole_lat,pole.long =  gr_npole_lon)[2]
london_lat_rot <- pp.ll.to.rg(lat=london_lat,long = london_lon, gr_npole_lat, gr_npole_lon)[1]
london_rot_sf <- data.frame(lon=london_lon_rot,lat=london_lat_rot) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  summarise(geometry = st_combine(geometry))
# find nearest dot
# st_distance(uk_temp_sf,london_rot_sf)[,1] %>% which.min()
# uk_london_dist <- (uk_temp_sf %>% mutate(d=st_distance(uk_temp_sf,london_rot_sf)[,1]))
# tm_shape(uk_london_dist) +
  # tm_dots("d",style="cont") +
  # tm_shape(london_rot_sf) + tm_dots() +
  # tm_shape(uk_london_dist[st_distance(uk_temp_sf,london_rot_sf)[,1] %>% which.min(),]) + tm_dots()
# calculate distances from London grid point
grid_london <- uk_temp_sf[st_distance(uk_temp_sf,london_rot_sf)[,1] %>% which.min(),]
# st_distance(uk_temp_sf,grid_london)[,1]
# check
uk_1999_sf <-  uk_temp_sf %>% mutate("dist_london"=st_distance(uk_temp_sf,grid_london)[,1])  
# check  
# tm_shape(uk_1999_sf) + tm_dots("dist",style="cont")
is_location <- rep("no",dim(uk_1999_sf)[1])
is_location[uk_1999_sf$dist_london==set_units(0,m)] <- "london"
long <- st_coordinates(uk_1999_sf %>% st_transform(4326))[,1]
lat <- st_coordinates(uk_1999_sf %>% st_transform(4326))[,2]
uk_1999 <-( uk_1999_sf %>% as.data.frame() %>%
                     select(-geometry) %>%
                     mutate(is_location=is_location)%>% 
              mutate(Longitude=long) %>% 
              mutate(Latitude=lat) %>% 
              relocate(dist_london) %>% 
                     relocate(c(Longitude,Latitude,is_location),.before = dist_london))

### create birmingham dataset ----
# coord for birmingham centre (Birmingham Art Gallery)
birmingham_lat <- 52.4806
birmingham_lon <- -1.9032
# plot to check
# birmingham_sf <- data.frame(lon=birmingham_lon,lat=birmingham_lat) %>%
  # st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  # summarise(geometry = st_combine(geometry))
# tm_shape(birmingham_sf) + tm_dots()
# rotate coordinates
birmingham_lon_rot <- pp.ll.to.rg(lat=birmingham_lat,long = birmingham_lon,pole.lat =  gr_npole_lat,pole.long =  gr_npole_lon)[2]
birmingham_lat_rot <- pp.ll.to.rg(lat=birmingham_lat,long = birmingham_lon, gr_npole_lat, gr_npole_lon)[1]
birmingham_rot_sf <- data.frame(lon=birmingham_lon_rot,lat=birmingham_lat_rot) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  summarise(geometry = st_combine(geometry))
# find nearest dot
# st_distance(uk_temp_sf,birmingham_rot_sf)[,1] %>% which.min()
# uk_birmingham_dist <- (uk_temp_sf %>% mutate(d=st_distance(uk_temp_sf,birmingham_rot_sf)[,1]))
# tm_shape(uk_birmingham_dist) +
#   tm_dots("d",style="cont") +
#   tm_shape(birmingham_rot_sf) + tm_dots() +
#   tm_shape(uk_birmingham_dist[st_distance(uk_temp_sf,birmingham_rot_sf)[,1] %>% which.min(),]) + tm_dots()
# calculate distances from birmingham grid point
grid_birmingham <- uk_temp_sf[st_distance(uk_temp_sf,birmingham_rot_sf)[,1] %>% which.min(),]
# st_distance(uk_temp_sf,grid_birmingham)[,1]
# check
uk_1999_sf <-  uk_temp_sf %>% mutate("dist_birmingham"=st_distance(uk_temp_sf,grid_birmingham)[,1])  
# check  
# tm_shape(uk_1999_sf) + tm_dots("dist_birmingham",style="cont")
uk_1999$is_location[uk_1999_sf$dist_birmingham==set_units(0,m)] <- "birmingham"
uk_1999 <- uk_1999 %>%
  mutate("dist_birmingham"=st_distance(uk_temp_sf,grid_birmingham)[,1]) %>% 
  relocate(dist_birmingham, .before = Temperature)


### create glasgow dataset ----
# coord for glasgow centre (Glasgow Central)
glasgow_lat <- 55.859112
glasgow_lon <- -4.258109
# plot to check
#  glasgow_sf <- data.frame(lon=glasgow_lon,lat=glasgow_lat) %>%
# st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
#   summarise(geometry = st_combine(geometry))
# tm_shape(glasgow_sf) + tm_dots()
# rotate coordinates
glasgow_lon_rot <- pp.ll.to.rg(lat=glasgow_lat,long = glasgow_lon,pole.lat =  gr_npole_lat,pole.long =  gr_npole_lon)[2]
glasgow_lat_rot <- pp.ll.to.rg(lat=glasgow_lat,long = glasgow_lon, gr_npole_lat, gr_npole_lon)[1]
glasgow_rot_sf <- data.frame(lon=glasgow_lon_rot,lat=glasgow_lat_rot) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  summarise(geometry = st_combine(geometry))
# find nearest dot
# st_distance(uk_temp_sf,glasgow_rot_sf)[,1] %>% which.min()
# uk_glasgow_dist <- (uk_temp_sf %>% mutate(d=st_distance(uk_temp_sf,glasgow_rot_sf)[,1]))
# tm_shape(uk_glasgow_dist) +
#   tm_dots("d",style="cont") +
#   tm_shape(glasgow_rot_sf) + tm_dots() +
#   tm_shape(uk_glasgow_dist[st_distance(uk_temp_sf,glasgow_rot_sf)[,1] %>% which.min(),]) + tm_dots()
# calculate distances from glasgow grid point
grid_glasgow <- uk_temp_sf[st_distance(uk_temp_sf,glasgow_rot_sf)[,1] %>% which.min(),]
# st_distance(uk_temp_sf,grid_glasgow)[,1]
# check
uk_1999_sf <-  uk_temp_sf %>% mutate("dist_glasgow"=st_distance(uk_temp_sf,grid_glasgow)[,1])  
# check  
# tm_shape(uk_1999_sf) + tm_dots("dist_glasgow",style="cont")
uk_1999$is_location[uk_1999_sf$dist_glasgow==set_units(0,m)] <- "glasgow"
uk_1999 <- uk_1999 %>%
  mutate("dist_glasgow"=st_distance(uk_temp_sf,grid_glasgow)[,1]) %>% 
  relocate(dist_glasgow, .before = Temperature)

### do for the other years -----
# start with 1999 year
no_col <- ncol(uk_1999)
for (i in (no_col-359):no_col) {
names(uk_1999)[i] <- paste0(names(uk_1999)[i],"_",1990)
}
uk_1999_2018 <- uk_1999
uk_1999_2018_winter <- uk_1999_2018[,c(1:6,7:(90+6))]
uk_1999_2018_spring <- uk_1999_2018[,c(1:6,(7+90):(180+6))]
uk_1999_2018_summer <- uk_1999_2018[,c(1:6,(7+180):(270+6))]
uk_1999_2018_autumn <- uk_1999_2018[,c(1:6,(7+270):(360+6))]

for (j in 2000:2018) {
# load data
  flist  <- paste0("data/tasmax_rcp85_land-cpm_uk_2.2km_01_day_",j,"1201-",j+1,"1130.nc")
  ### read netcdf file
  nc1 <- nc_open(flist[1])
  # time coordinate
  time1    <- ncvar_get(nc1,'time')
  time1att <- ncatt_get(nc1,'time')
  # read all 360 days
  # time dimension is last in the start argument
  v1    <- ncvar_get(nc1,'tasmax',start=c(1,1,1,1), count=c(-1,-1,360,1))
  # lat long of grid
  rlon     <- apply(ncvar_get(nc1,"grid_longitude_bnds"),2,mean)
  rlat     <- apply(ncvar_get(nc1,"grid_latitude_bnds"),2,mean)
  cpm.pole <- ncatt_get(nc1,'rotated_latitude_longitude')
  nc_close(nc1)
  # this climate model has a 360 day year so needs special treatment
  # PCICt does this for us
  origin   <- tail(unlist(strsplit(time1att$units,'hours since ')),1)
  pcictime <- as.PCICt( time1*(60*60), cal=time1att$calendar ,origin=origin)
  
  # the CPM grid has the pole in a different place to usual.  therefore need to transform the map coordinates
  gr_npole_lat <- cpm.pole$grid_north_pole_latitude
  gr_npole_lon <- cpm.pole$grid_north_pole_longitude
  
# sparse grid
  longitude <- rlon
  latitude <- rlat
  ### subset the grid
  lon_subset <- seq(from=10,to=(length(longitude)-length(longitude)%%10),by=10)
  lat_subset <- seq(from=10,to=(length(latitude)-length(latitude)%%10),by=10)
  longitude <- longitude[lon_subset]
  latitude <- latitude[lat_subset]
  v1_sub <- v1[lon_subset,lat_subset,]
  #fields::image.plot( x=longitude, y=latitude, z=v1_sub[,,1], zlim=c(0,20), main=pcictime[2])
  
  ### subset rotated grid points overlapping with UK -----
  # create grid points for long lat
  lon_lat <- expand.grid(longitude,latitude) 
  colnames(lon_lat) <- c("lon","lat")
  # extract first day (first out of 360 in dim3:time)
  Temperature <- c(v1_sub[,,1])
  lon_lat_temp <- cbind(lon_lat,Temperature)
  for (i in 2:(dim(v1_sub)[3]/1)) {
    lon_lat_temp <- cbind(lon_lat_temp,c(v1_sub[,,i]))
    colnames(lon_lat_temp)[length(lon_lat_temp)] <- as.character(i)
  }
  # convert to sf points object
  temp_sf <- st_as_sf(lon_lat_temp,coords = c("lon","lat"),crs=4326)
  # subset to only include grid points in mainland UK
  uk_temp_sf <- st_filter(temp_sf,(uk_rot %>% st_cast("MULTIPOLYGON")))

# combine with data set
temp_to_cbind <- uk_temp_sf  %>% as.data.frame()%>% 
  select(-geometry) %>% select(all_of(1:90))
for (i in 1:90) {
  names(temp_to_cbind)[i] <- paste0(names(temp_to_cbind)[i],"_",j)
}
uk_1999_2018_winter <- cbind(uk_1999_2018_winter,temp_to_cbind)

temp_to_cbind <- uk_temp_sf  %>% as.data.frame()%>% 
  select(-geometry) %>% select(all_of(91:180))
for (i in 91:180) {
  names(temp_to_cbind)[(i-90)] <- paste0(names(temp_to_cbind)[(i-90)],"_",j)
}
uk_1999_2018_spring <- cbind(uk_1999_2018_spring,temp_to_cbind)

temp_to_cbind <- uk_temp_sf  %>% as.data.frame()%>% 
  select(-geometry) %>% select(all_of(181:270))
for (i in 1:90) {
  names(temp_to_cbind)[i] <- paste0(names(temp_to_cbind)[i],"_",j)
}
uk_1999_2018_summer <- cbind(uk_1999_2018_summer,temp_to_cbind)

temp_to_cbind <- uk_temp_sf  %>% as.data.frame()%>% 
  select(-geometry) %>% select(all_of(271:360))
for (i in 1:90) {
  names(temp_to_cbind)[i] <- paste0(names(temp_to_cbind)[i],"_",j)
}
uk_1999_2018_autumn <- cbind(uk_1999_2018_autumn,temp_to_cbind)
}
# uk_1999_2018_winter %>% names()
# plot(x=6:ncol(uk_1999_2018_winter),y=uk_1999_2018_winter[1,6:ncol(uk_1999_2018_winter)])
# save as R object
 # saveRDS(uk_1999_2018_winter, "data/uk_1999_2018_winter.RDS")
 # saveRDS(uk_1999_2018_spring, "data/uk_1999_2018_spring.RDS")
 # saveRDS(uk_1999_2018_summer, "data/uk_1999_2018_summer.RDS")
 # saveRDS(uk_1999_2018_autumn, "data/uk_1999_2018_autumn.RDS")
# readRDS("data/uk_1999_2018_winter.RDS")

# start with 1999 year and make a function of end year
add_temperature_years <- function(uk_1999,last_year_start) {
no_col <- ncol(uk_1999)
tmp <- uk_1999
for (i in (no_col-359):no_col) {
  names(tmp)[i] <- paste0(names(tmp)[i],"_",1990)
}

uk_1999_end[[1]] <- tmp[,c(1:6,7:(90+6))]
uk_1999_end[[2]] <- tmp[,c(1:6,(7+90):(180+6))]
uk_1999_end[[3]] <- tmp[,c(1:6,(7+180):(270+6))]
uk_1999_end[[4]] <- tmp[,c(1:6,(7+270):(360+6))]

for (j in 2000:last_year_start) {
  # load data
  flist  <- paste0("data/tasmax_rcp85_land-cpm_uk_2.2km_01_day_",j,"1201-",j+1,"1130.nc")
  ### read netcdf file
  nc1 <- nc_open(flist[1])
  # time coordinate
  time1    <- ncvar_get(nc1,'time')
  time1att <- ncatt_get(nc1,'time')
  # read all 360 days
  # time dimension is last in the start argument
  v1    <- ncvar_get(nc1,'tasmax',start=c(1,1,1,1), count=c(-1,-1,360,1))
  # lat long of grid
  rlon     <- apply(ncvar_get(nc1,"grid_longitude_bnds"),2,mean)
  rlat     <- apply(ncvar_get(nc1,"grid_latitude_bnds"),2,mean)
  cpm.pole <- ncatt_get(nc1,'rotated_latitude_longitude')
  nc_close(nc1)
  # this climate model has a 360 day year so needs special treatment
  # PCICt does this for us
  origin   <- tail(unlist(strsplit(time1att$units,'hours since ')),1)
  pcictime <- as.PCICt( time1*(60*60), cal=time1att$calendar ,origin=origin)
  
  # the CPM grid has the pole in a different place to usual.  therefore need to transform the map coordinates
  gr_npole_lat <- cpm.pole$grid_north_pole_latitude
  gr_npole_lon <- cpm.pole$grid_north_pole_longitude
  
  # sparse grid
  longitude <- rlon
  latitude <- rlat
  ### subset the grid
  lon_subset <- seq(from=10,to=(length(longitude)-length(longitude)%%10),by=10)
  lat_subset <- seq(from=10,to=(length(latitude)-length(latitude)%%10),by=10)
  longitude <- longitude[lon_subset]
  latitude <- latitude[lat_subset]
  v1_sub <- v1[lon_subset,lat_subset,]
  #fields::image.plot( x=longitude, y=latitude, z=v1_sub[,,1], zlim=c(0,20), main=pcictime[2])
  
  ### subset rotated grid points overlapping with UK -----
  # create grid points for long lat
  lon_lat <- expand.grid(longitude,latitude) 
  colnames(lon_lat) <- c("lon","lat")
  # extract first day (first out of 360 in dim3:time)
  Temperature <- c(v1_sub[,,1])
  lon_lat_temp <- cbind(lon_lat,Temperature)
  for (i in 2:(dim(v1_sub)[3]/1)) {
    lon_lat_temp <- cbind(lon_lat_temp,c(v1_sub[,,i]))
    colnames(lon_lat_temp)[length(lon_lat_temp)] <- as.character(i)
  }
  # convert to sf points object
  temp_sf <- st_as_sf(lon_lat_temp,coords = c("lon","lat"),crs=4326)
  # subset to only include grid points in mainland UK
  uk_temp_sf <- st_filter(temp_sf,(uk_rot %>% st_cast("MULTIPOLYGON")))
  
  # combine with data set
  temp_to_cbind <- uk_temp_sf  %>% as.data.frame()%>% 
    select(-geometry) %>% select(all_of(1:90))
  for (i in 1:90) {
    names(temp_to_cbind)[i] <- paste0(names(temp_to_cbind)[i],"_",j)
  }
  uk_1999_end[[1]] <- cbind(uk_1999_end[[1]],temp_to_cbind)
  
  temp_to_cbind <- uk_temp_sf  %>% as.data.frame()%>% 
    select(-geometry) %>% select(all_of(91:180))
  for (i in 91:180) {
    names(temp_to_cbind)[(i-90)] <- paste0(names(temp_to_cbind)[(i-90)],"_",j)
  }
  uk_1999_end[[2]] <- cbind(uk_1999_end[[2]],temp_to_cbind)
  
  temp_to_cbind <- uk_temp_sf  %>% as.data.frame()%>% 
    select(-geometry) %>% select(all_of(181:270))
  for (i in 1:90) {
    names(temp_to_cbind)[i] <- paste0(names(temp_to_cbind)[i],"_",j)
  }
  uk_1999_end[[3]] <- cbind(uk_1999_end[[3]],temp_to_cbind)
  
  temp_to_cbind <- uk_temp_sf  %>% as.data.frame()%>% 
    select(-geometry) %>% select(all_of(271:360))
  for (i in 1:90) {
    names(temp_to_cbind)[i] <- paste0(names(temp_to_cbind)[i],"_",j)
  }
  uk_1999_end[[4]] <- cbind(uk_1999_end[[4]],temp_to_cbind)
}
return(uk_1999_end)
}
t <- add_temperature_years(uk_1999,last_year_start = 2079)
# t[[1]] %>% names()
# plot(x=7:ncol(t[[3]]),y=t[[3]][1,7:ncol(t[[3]])])
d <- t[[3]][1,7:ncol(t[[3]])]
# save as R object
# saveRDS(t[[1]], "data/uk_1999_2079_winter.RDS")
# saveRDS(t[[2]], "data/uk_1999_2079_spring.RDS")
# saveRDS(t[[3]], "data/uk_1999_2079_summer.RDS")
# saveRDS(t[[4]], "data/uk_1999_2079_autumn.RDS")
# readRDS("data/uk_1999_2079_winter.RDS")


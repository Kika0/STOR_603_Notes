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

fields::image.plot( x=longitude, y=latitude, z=v1_sub[,,1], zlim=c(0,20), main=pcictime[2])
points(x=long,y=lat)
### subset to show only Lancaster district ----
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

### subset rotated grid points overlapping with UK -----
# create grid points for long lat
lon_lat <- expand.grid(longitude,latitude) 
colnames(lon_lat) <- c("lon","lat")
# extract first half-hour of the day (first out of 48 in dim3:time)
Temperature <- c(v1_sub[,,1])
lon_lat_temp <- cbind(lon_lat,Temperature)
for (i in 2:(dim(v1_sub)[3]/1)) {
  lon_lat_temp <- cbind(lon_lat_temp,c(v1_sub[,,i]))
  colnames(lon_lat_temp)[length(lon_lat_temp)] <- as.character(i)
}
# convert to sf points object
temp_sf <- st_as_sf(lon_lat_temp,coords = c("lon","lat"),crs=4326)
tm_shape(temp_sf) + tm_dots(col="Temperature",style="cont",palette="viridis")
tm_shape(temp_sf) + tm_dots(col="Temperature",style="cont",size=0.05,palette="-RdYlBu")+
  tm_layout(legend.position = c(0.7, 0.55))
# subset to only include grid points in mainland UK
uk_temp_sf <- st_filter(temp_sf,(uk_rot %>% st_cast("MULTIPOLYGON")))
tm_shape(uk_temp_sf) + tm_dots(col="Temperature",style="cont",size=0.1)

# try to plot for 30 days ----
# subset to only include grid points in Lancaster
uk_temp_sf_long <- uk_temp_sf %>% select(Temperature) %>%
  mutate("day"=as.factor(rep(1,nrow(uk_temp_sf)))) %>% 
  mutate("date"=rep(pcictime[1],nrow(uk_temp_sf)))
# December subset 
pcictime[1:30]

# for loop for all other half-hour intervals
for (j in (2:(dim(v1_sub)[3]/1))) {
  to_bind <- uk_temp_sf  %>% select(all_of(j)) %>% 
    mutate("day"=as.factor(rep(j,nrow(uk_temp_sf)))) %>% 
    mutate("date"=rep(pcictime[j],nrow(uk_temp_sf)))
  names(to_bind)[1] <- "Temperature"
  uk_temp_sf_long <- rbind(uk_temp_sf_long,to_bind)
  
}

tm_shape(uk_temp_sf_long) + 
  tm_dots(col="Temperature",style="cont",size=0.05,palette="viridis") +
  tm_facets(by="day",as.layers=TRUE,ncol=30,nrow=12) +
  tm_layout(panel.show = FALSE)
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
london_sf <- data.frame(lon=london_lon,lat=london_lat) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  summarise(geometry = st_combine(geometry))
tm_shape(london_sf) + tm_dots()
# rotate coordinates
london_lon_rot <- pp.ll.to.rg(lat=london_lat,long = london_lon,pole.lat =  gr_npole_lat,pole.long =  gr_npole_lon)[2]
london_lat_rot <- pp.ll.to.rg(lat=london_lat,long = london_lon, gr_npole_lat, gr_npole_lon)[1]
london_rot_sf <- data.frame(lon=london_lon_rot,lat=london_lat_rot) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  summarise(geometry = st_combine(geometry))
# find nearest dot
st_distance(uk_temp_sf,london_rot_sf)[,1] %>% which.min()
uk_london_dist <- (uk_temp_sf %>% mutate(d=st_distance(uk_temp_sf,london_rot_sf)[,1]))
tm_shape(uk_london_dist) +
  tm_dots("d",style="cont") +
  tm_shape(london_rot_sf) + tm_dots() +
  tm_shape(uk_london_dist[st_distance(uk_temp_sf,london_rot_sf)[,1] %>% which.min(),]) + tm_dots()
# calculate distances from London grid point
grid_london <- uk_temp_sf[st_distance(uk_temp_sf,london_rot_sf)[,1] %>% which.min(),]
# st_distance(uk_temp_sf,grid_london)[,1]
# check
uk_1999_sf <-  uk_temp_sf %>% mutate("dist"=st_distance(uk_temp_sf,grid_london)[,1])  
# check  
# tm_shape(uk_1999_sf) + tm_dots("dist",style="cont")
is_london <- rep("not_london",dim(uk_1999_sf)[1])
is_london[uk_1999_sf$dist==set_units(0,m)] <- "london"
uk_1999_winter <-( uk_1999_sf %>% as.data.frame() %>%
                     select(geometry,dist,everything()) %>%
                     mutate(is_london=is_london))[,1:92]

### create birmingham dataset ----
# coord for birmingham centre (Alan Turing Institute)
birmingham_lat <- 52.4806
birmingham_lon <- -1.9032
# plot to check
birmingham_sf <- data.frame(lon=birmingham_lon,lat=birmingham_lat) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  summarise(geometry = st_combine(geometry))
tm_shape(birmingham_sf) + tm_dots()
# rotate coordinates
birmingham_lon_rot <- pp.ll.to.rg(lat=birmingham_lat,long = birmingham_lon,pole.lat =  gr_npole_lat,pole.long =  gr_npole_lon)[2]
birmingham_lat_rot <- pp.ll.to.rg(lat=birmingham_lat,long = birmingham_lon, gr_npole_lat, gr_npole_lon)[1]
birmingham_rot_sf <- data.frame(lon=birmingham_lon_rot,lat=birmingham_lat_rot) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  summarise(geometry = st_combine(geometry))
# find nearest dot
st_distance(uk_temp_sf,birmingham_rot_sf)[,1] %>% which.min()
uk_birmingham_dist <- (uk_temp_sf %>% mutate(d=st_distance(uk_temp_sf,birmingham_rot_sf)[,1]))
tm_shape(uk_birmingham_dist) +
  tm_dots("d",style="cont") +
  tm_shape(birmingham_rot_sf) + tm_dots() +
  tm_shape(uk_birmingham_dist[st_distance(uk_temp_sf,birmingham_rot_sf)[,1] %>% which.min(),]) + tm_dots()
# calculate distances from birmingham grid point
grid_birmingham <- uk_temp_sf[st_distance(uk_temp_sf,birmingham_rot_sf)[,1] %>% which.min(),]
# st_distance(uk_temp_sf,grid_birmingham)[,1]
# check
uk_1999_sf <-  uk_temp_sf %>% mutate("dist"=st_distance(uk_temp_sf,grid_birmingham)[,1])  
# check  
# tm_shape(uk_1999_sf) + tm_dots("dist",style="cont")
is_london <- rep("not_london",dim(uk_1999_sf)[1])
is_london[uk_1999_sf$dist==set_units(0,m)] <- "birmingham"
uk_1999_winter <-( uk_1999_sf %>% as.data.frame() %>%
                     select(geometry,dist,everything()) %>%
                     mutate(is_london=is_london))[,1:92]

# calculate dependence between X(London) and Y(some other location)
X <- uk_1999_winter[is_london=="london",3:ncol(uk_1999_winter)] 
Y <- uk_1999_winter[1,3:ncol(uk_1999_winter)] 




library(ncdf4)
library(fields)
library(maps)
library(tmap)
library(mapdata)
library(PCICt)
library(rgdal)
# for palettes
library(viridis)
library(tidyverse)
library(sf)
#library(gridExtra)

# function to translate from normal pole coordinates to rotated pole coordinates
pp.ll.to.rg <- function(lat,long,pole.lat,pole.long) {
  while(pole.long>180) pole.long<-pole.long-360
  l0  <- pole.long+180
  dtr <- pi/180
  sin.pole.lat <- sin(pole.lat*dtr)
  cos.pole.lat <- cos(pole.lat*dtr)
  if(pole.lat < 0) {
    sin.pole.lat <- -sin.pole.lat
    cos.pole.lat <- -cos.pole.lat
  }
  long <- long-l0
  while(long >  180) long <- long-360
  while(long < -180) long <- long+360
  
  lat.rotated <- asin(max(-1,min(1,-cos.pole.lat*
                                   cos(long*dtr)*
                                   cos(lat*dtr)+
                                   sin.pole.lat*
                                   sin(lat*dtr))))
  
  long.rotated <- 0
  if(cos(lat.rotated) > 1.0e-6) {
    long.rotated <- acos(max(-1,min(1,(cos.pole.lat*
                                         sin(lat*dtr)+
                                         sin.pole.lat*
                                         cos(long*dtr)*
                                         cos(lat*dtr))/
                                      cos(lat.rotated))))
  }
  long.rotated <- long.rotated*sign(long)
  lat.rotated  <- lat.rotated/dtr
  long.rotated <- long.rotated/dtr
  while(long.rotated >  180) long.rotated <- long.rotated-360
  while(long.rotated < -180) long.rotated <- long.rotated+360
  return(c(lat.rotated,long.rotated))
}
cr <- '\n'

### read one year of data
flist  <- "data/tasmax_rcp85_land-cpm_uk_2.2km_01_day_19991201-20001130.nc"

### read netcdf file
cat("loading",basename(flist[1]),cr) 
nc1 <- nc_open(flist[1])
# all variables and metadata in file
vlist <- nc1$var
# time coordinate
time1    <- ncvar_get(nc1,'time')
time1att <- ncatt_get(nc1,'time')
# rainfall rate
# read all 360 days
# time dimension is last in the start argument
v1    <- ncvar_get(nc1,'tasmax',start=c(1,1,1,1), count=c(-1,-1,360,1))
# names(vlist) contains list of variables 
v1att <- ncatt_get(nc1,'tasmax')
# lat long of grid
rlon     <- apply(ncvar_get(nc1,"grid_longitude_bnds"),2,mean)
rlat     <- apply(ncvar_get(nc1,"grid_latitude_bnds"),2,mean)
cpm.pole <- ncatt_get(nc1,'rotated_latitude_longitude')
nc_close(nc1)
cat("Variable:",v1att$standard_name,cr)
cat("Units:",v1att$units,cr)

# this climate model has a 360 day year so needs special treatment
# PCICt does this for us
origin   <- tail(unlist(strsplit(time1att$units,'hours since ')),1)
pcictime <- as.PCICt( time1*(60*60), cal=time1att$calendar ,origin=origin)

# the CPM grid has the pole in a different place to usual.  therefore need to transform the map coordinates
gr_npole_lat <- cpm.pole$grid_north_pole_latitude
gr_npole_lon <- cpm.pole$grid_north_pole_longitude
n.map <- maps::map('worldHires', xlim=c(-10,10),ylim=c(30,70),interior=F,plot=FALSE)
r.map <- n.map
for (i in 1:length(r.map$x)){
  if(!is.na(n.map$y[i]) & !is.na(n.map$x[i])) r.latlon   <- pp.ll.to.rg(n.map$y[i],n.map$x[i], gr_npole_lat, gr_npole_lon)
  r.map$x[i] <- r.latlon[2]
  r.map$y[i] <- r.latlon[1]
}


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
points(r.map$x+360,r.map$y,pch=46,col='grey80', lwd=0)

### subset to show only Lancaster district ----
#m <- st_read("cnty/infuse_cnty_lyr_2011.shp")
lad <- st_read("data/LAD_boundary_UK/LAD_MAY_2022_UK_BFE_V3.shp")
# UK is comprised of many polygons (islands), simplify to only
# take mainland UK (Great Britain)
uk <- (lad %>% st_union() %>% st_cast( "MULTIPOLYGON" ) %>% st_cast("POLYGON"))[1019]
uk <- st_simplify(uk,dTolerance = 7000)
# rotate coordinates of uk polygon
long <- st_coordinates(uk %>% st_transform(4326))[,1]
lat <- st_coordinates(uk %>% st_transform(4326))[,2]

#tmp <- (lad %>% st_union() %>% st_cast( "MULTIPOLYGON" ) %>% st_cast("POLYGON"))
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
# check by mapping
tm_shape(uk_rot) + tm_polygons()

# study point rotation
# df1 <- data.frame(X=long_rot,Y=lat_rot,type=rep("Rotated",length(long_rot)))
# df2 <- data.frame(X=long,Y=lat,type=rep("Original",length(long)))
# #df <- rbind(df1,df2)
# ggplot() + geom_point(data=df1,mapping=aes(x=X,y=Y))+ geom_point(data=df2,mapping=aes(x=X,y=Y)) 
# can play further to only include UK shapefile

### subset rotated grid points overlapping with UK -----
# create grid points for long lat
lon_lat <- expand.grid(longitude,latitude) 
colnames(lon_lat) <- c("lon","lat")
# extract first half-hour of the day (first out of 48 in dim3:time)
Temperature <- c(v1_sub[,,1])
lon_lat_temp <- cbind(lon_lat,Temperature)
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
for (j in (2:30)) {
  to_bind <- uk_temp_sf  %>% select(all_of(j)) %>% 
    mutate("day"=as.factor(rep(j,nrow(uk_temp_sf)))) %>% 
    mutate("date"=rep(pcictime[j],nrow(uk_temp_sf)))
  names(to_bind)[1] <- "Temperature"
  uk_temp_sf_long <- rbind(uk_temp_sf_long,to_bind)
  
}

lanc_temp_sf_long <- rbind(lanc_temp_sf_long,to_bind)
tm_shape(lanc_temp_sf_long) + 
  tm_dots(col="Temperature",style="cont",size=0.9,palette="viridis") +
  tm_facets(by="day",as.layers=TRUE,ncol=8,nrow=6)
# facet by timestamp
lanc_temp_sf_long$date <- factor(uk_temp_sf_long$date,      # Reordering group factor levels
                                 levels = (uk_temp_sf_long$time %>% unique()))
tm_shape(uk_temp_sf_long) + 
  tm_dots(col="Temperature",style="cont",size=0.9,palette="viridis",title = "Temperature") +
  tm_facets(by="time",as.layers=TRUE,ncol=6,nrow=5) 

# to plot mean temperature in Lancaster
t <- c()
for (i in (1:48)) {
  t[i] <- mean(lanc_temp_sf_long$temp[(132*(i-1)+1):(132*i)])
}
time <- 1:48
plot(time,t,type="l")




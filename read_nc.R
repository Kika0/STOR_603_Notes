# source("read_nc.R")

library(ncdf4)
library(fields)
library(maps)
library(tmap)
library(mapdata)
library(PCICt)
library(rgdal)
library(viridis) # for palettes
library(ggthemes)
library(tidyverse)
library(sf)
library(gridExtra)
source("rotate_unrotate_coordinates.R")
#flist  <- "/data/users/hadsx/model_data/cpm/halfhourly/leeds_msc/r001i1p00000_19991201-19991230_pr.nc"
# update file path
flist  <- "data/tasmax_rcp85_land-cpm_uk_2.2km_01_day_19991201-20001130.nc"
cr <- '\n'

### read netcdf file
cat("loading",basename(flist[1]),cr) 
nc1 <- nc_open(flist[1])
    # all variables and metadata in file
    vlist <- nc1$var
    # time coordinate
    time1    <- ncvar_get(nc1,'time')
    time1att <- ncatt_get(nc1,'time')
    # rainfall rate
    # read one days worth of 1/2 hourly data
    v1    <- ncvar_get(nc1,'tasmax',start=c(1,1,1,1), count=c(-1,-1,48,1))
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
fields::image.plot( x=longitude, y=latitude, z=v1[,,2], zlim=c(0,20), main=pcictime[2])
points(r.map$x+360,r.map$y,pch=46,col='grey80', lwd=0.1)

## try to just plot a subset
image.plot( x=rlon, y=rlat[400:606], z=v1[,,1][,400:606], zlim=c(0,20), main=pcictime[1])
points(r.map$x+360,r.map$y,pch=46,col='grey80', lwd=0.1)


### data structure -----
# > str(v1)
#  num [1:532, 1:654, 1:48] 0.172 0.164 0.151 0.133 0.116 ...
# dim1 is x/longitude
# dim2 is y/latitude
# dim3 is time

### subset to show only Lancaster district ----
#m <- st_read("cnty/infuse_cnty_lyr_2011.shp")
lad <- st_read("data/LAD_boundary_UK/LAD_MAY_2022_UK_BFE_V3.shp")
lanc <- lad %>% filter(LAD22NM=="Lancaster")


# rotate coordinates of lanc
long <- st_coordinates(lanc %>% st_transform(4326))[,1]
lat <- st_coordinates(lanc %>% st_transform(4326))[,2]

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
lanc_rot <-  df %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")
# check by mapping
tm_shape(lanc_rot) + tm_polygons()



# compare results with random sample
# n.map <- maps::map('worldHires', xlim=c(-10,10),ylim=c(30,70),interior=F,plot=FALSE)
# n.map$x <- sample(n.map$x,20)
# n.map$y <- sample(n.map$y,20)
# r.map <- n.map
# for (i in 1:length(r.map$x)){
#   if(!is.na(n.map$y[i]) & !is.na(n.map$x[i])) r.latlon   <- pp.ll.to.rg(n.map$y[i],n.map$x[i], gr_npole_lat, gr_npole_lon)
#   r.map$x[i] <- r.latlon[2]
#   r.map$y[i] <- r.latlon[1]
# }
n.map <- maps::map('worldHires', xlim=c(-10,3),ylim=c(48,70),interior=F,plot=FALSE)
r.map <- n.map
for (i in 1:length(r.map$x)){
  if(!is.na(n.map$y[i]) & !is.na(n.map$x[i])) r.latlon   <- pp.ll.to.rg(n.map$y[i],n.map$x[i], gr_npole_lat, gr_npole_lon)
  r.map$x[i] <- r.latlon[2]
  r.map$y[i] <- r.latlon[1]
}
# study point rotation
df1 <- data.frame(X=long_rot,Y=lat_rot,type=rep("Rotated",length(long_rot)))
df2 <- data.frame(X=long,Y=lat,type=rep("Original",length(long)))
#df <- rbind(df1,df2)
ggplot() + geom_point(data=df1,mapping=aes(x=X,y=Y))+ geom_point(data=df2,mapping=aes(x=X,y=Y)) 

# study whole UK since Lancaster district only is squished
df1 <- data.frame(X=n.map$x,Y=n.map$y,Type=rep("Original",length(n.map$x)))
df2 <- data.frame(X=r.map$x,Y=r.map$y,Type=rep("Rotated",length(r.map$x)))
#df3 <- cbind(df1,df2)
df <- rbind(df1,df2)
ggplot() + geom_point(data=df1,mapping=aes(x=X,y=Y))
ggplot(df) + geom_point(aes(x=X,y=Y,col=Type),lwd=0.1) + coord_fixed()
rm(df)
# can play further to only include UK shapefile

### subset rotated grid points overlapping with Lancaster -----
# create grid points for long lat
lon_lat <- expand.grid(rlon,rlat) 
colnames(lon_lat) <- c("lon","lat")
# extract first half-hour of the day (first out of 48 in dim3:time)
Temperature <- c(v1[,,1])
lon_lat_temp <- cbind(lon_lat,Temperature)
# convert to sf points object
temp_sf <- st_as_sf(lon_lat_temp,coords = c("lon","lat"),crs=4326)
tm_shape(temp_sf) + tm_dots(col="Temperature",style="cont",palette="viridis")
tm_shape(temp_sf) + tm_dots(col="Temperature",style="cont",size=0.05,palette="-RdYlBu")+
  tm_layout(legend.position = c(0.7, 0.55))
# subset to only include grid points in Lancaster
lanc_temp_sf <- st_filter(temp_sf,lanc_rot)
tm_shape(lanc_temp_sf) + tm_dots(col="Temperature",style="cont",size=2)

# try to plot for the whole day ----
temp <- c(v1[,,1])
lon_lat_temp <- cbind(lon_lat,temp)
for (i in 2:dim(v1)[3]) {
  lon_lat_temp <- cbind(lon_lat_temp,c(v1[,,i]))
  colnames(lon_lat_temp)[length(lon_lat_temp)] <- as.character(i)
}
# convert to sf points object
temp_sf <- st_as_sf(lon_lat_temp,coords = c("lon","lat"),crs=4326)
tm_shape(temp_sf) + tm_dots(col="temp",style="cont")

# subset to only include grid points in Lancaster
lanc_temp_sf <- st_filter(temp_sf,lanc_rot)

# try to plot for the whole day ----
# subset to only include grid points in Lancaster
lanc_temp_sf <- st_filter(temp_sf,lanc_rot)

lanc_temp_sf_long <- lanc_temp_sf %>% select(Temperature) %>%
  mutate("time_of_day"=as.factor(rep(1,nrow(lanc_temp_sf)))) %>% 
  mutate("time"=rep(strftime(seq(from = as.POSIXct("1999-12-01 12:00"), 
                        to = as.POSIXct("1999-12-02 11:30"), by = "30 min")[1],format="%H:%M"),nrow(lanc_temp_sf)))
# subset only time
strftime(seq(from=as.POSIXct("1999-12-01 12:00"), 
         to = as.POSIXct("1999-12-02 11:30"), by = "30 min")[1], format="%H:%M")

# for loop for all other half-hour intervals
for (j in (2:48)) {
  to_bind <- lanc_temp_sf  %>% select(all_of(j)) %>% 
    mutate("time_of_day"=as.factor(rep(j,nrow(lanc_temp_sf)))) %>% 
    mutate("time"=rep(strftime(seq(from = as.POSIXct("1999-12-01 12:00"), 
        to = as.POSIXct("1999-12-02 11:30"), by = "30 min")[j],format="%H:%M"),nrow(lanc_temp_sf)))
  names(to_bind)[1] <- "temp"
  lanc_temp_sf_long <- rbind(lanc_temp_sf_long,to_bind)
  
}

lanc_temp_sf_long <- rbind(lanc_temp_sf_long,to_bind)
tm_shape(lanc_temp_sf_long) + 
  tm_dots(col="temp",style="cont",size=0.9,palette="viridis") +
  tm_facets(by="time_of_day",as.layers=TRUE,ncol=8,nrow=6)
# facet by timestamp
lanc_temp_sf_long$time <- factor(lanc_temp_sf_long$time,      # Reordering group factor levels
                         levels = (lanc_temp_sf_long$time %>% unique()))
tm_shape(lanc_temp_sf_long) + 
  tm_dots(col="temp",style="cont",size=0.9,palette="viridis",title = "Temperature") +
  tm_facets(by="time",as.layers=TRUE,ncol=8,nrow=6) 

# to plot mean temperature in Lancaster
t <- c()
for (i in (1:48)) {
  t[i] <- mean(lanc_temp_sf_long$temp[(132*(i-1)+1):(132*i)])
}
time <- 1:48
plot(time,t,type="l")

#----------------------------------------------------------------
### additional code unused
## Get land mask based off of UKCP high-res shapefile -----
UK_shp = readOGR(dsn="data/LAD_boundary_UK/LAD_MAY_2022_UK_BFE_V3.shp")

latlong = "+init=epsg:4326"
UK_shp_ll = spTransform(UK_shp, CRS(latlong))
rm(UK_shp)
UK_shp_plt = as(UK_shp_ll[(UK_shp_ll$LAD22NM == "Lancaster"),],"SpatialPolygons")
tmp = lonlat_df
#lonlat_df (north pole orientated coordinates) is rotated boundary polygon
colnames(tmp) = c("lon", "lat")#, "ind")
coordinates(tmp) = ~lon+lat
tmp@proj4string = CRS(latlong)
over_land = over(tmp,as(UK_shp_plt,"SpatialPolygons"))
tmp = tmp[!is.na(over_land),]
rm(UK_shp_ll, tmp)


library(tmap)
tm_shape(m) + tm_polygons()
tm_shape(lad) + tm_polygons()
tm_shape(lad %>% filter(LAD22NM=="Lancaster")) + tm_polygons()

# source("read_nc.R")

library(ncdf4)
library(fields)
library(maps)
library(mapdata)
library(PCICt)
library(rgdal)
library(tidyverse)
library(sf)
library(gridExtra)

flist  <- "/data/users/hadsx/model_data/cpm/halfhourly/leeds_msc/r001i1p00000_19991201-19991230_pr.nc"
# update file path
flist  <- "data/tasmax_rcp85_land-cpm_uk_2.2km_01_day_19991201-20001130.nc"


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
fields::image.plot( x=rlon, y=rlat, z=v1[,,2], zlim=c(0,20), main=pcictime[2])
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

# subset to show only Lancaster district
m <- st_read("cnty/infuse_cnty_lyr_2011.shp")
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



# compare results with random sample
n.map <- maps::map('worldHires', xlim=c(-10,10),ylim=c(30,70),interior=F,plot=FALSE)
n.map$x <- sample(n.map$x,20)
n.map$y <- sample(n.map$y,20)
r.map <- n.map
for (i in 1:length(r.map$x)){
  if(!is.na(n.map$y[i]) & !is.na(n.map$x[i])) r.latlon   <- pp.ll.to.rg(n.map$y[i],n.map$x[i], gr_npole_lat, gr_npole_lon)
  r.map$x[i] <- r.latlon[2]
  r.map$y[i] <- r.latlon[1]
}
n.map <- maps::map('worldHires', xlim=c(-10,10),ylim=c(30,70),interior=F,plot=FALSE)
r.map <- n.map
for (i in 1:length(r.map$x)){
  if(!is.na(n.map$y[i]) & !is.na(n.map$x[i])) r.latlon   <- pp.ll.to.rg(n.map$y[i],n.map$x[i], gr_npole_lat, gr_npole_lon)
  r.map$x[i] <- r.latlon[2]
  r.map$y[i] <- r.latlon[1]
}
# study point rotation
df1 <- data.frame(X=long_rot,Y=lat_rot,type=rep("Rotated",length(long_rot)))
df2 <- data.frame(X=long,Y=lat,type=rep("Original",length(long)))
df <- rbind(df1,df2)
ggplot() + geom_point(data=df1,mapping=aes(x=X,y=Y))+ geom_point(data=df2,mapping=aes(x=X,y=Y)) 

# study by linking points
df1 <- data.frame(X=n.map$x,Y=n.map$y,type=rep("Rotated",length(n.map$x)))
df2 <- data.frame(X=r.map$x,Y=r.map$y,type=rep("Rotated",length(r.map$x)))
df3 <- cbind(df1,df2)
ggplot() + geom_point(data=df1,mapping=aes(x=X,y=Y))
ggplot(df) + geom_point(aes(x=X,y=Y,col=type))
rm(df)
# can play further to only include UK shapefile

## Get land mask based off of UKCP high-res shapefile

[(UK_shp_ll$geo_region == UK_shp_ll$geo_region[1]),]
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

CnvRttPol = function(latlon, spol_coor, option=1){
  if(option==1){
    lon = latlon[,1]
    lat = latlon[,2]
    
    lon = (lon*pi)/180
    lat = (lat*pi)/180
    
    spol_lon = spol_coor[1]
    spol_lat = spol_coor[2]
    
    theta = 90+spol_lat# Rotation around y-axis
    phi = spol_lon # Rotation around z-axis
    
    phi = (phi*pi)/180 # Convert degrees to radians
    theta = (theta*pi)/180
    
    x = cos(lon)*cos(lat) # Convert from spherical to cartesian coordinates
    y = sin(lon)*cos(lat)
    z = sin(lat)
    
    phi = -phi
    theta = -theta
    
    x_new = cos(theta)*cos(phi)*x + sin(phi)*y + sin(theta)*cos(phi)*z;
    y_new = -cos(theta)*sin(phi)*x + cos(phi)*y - sin(theta)*sin(phi)*z;
    z_new = -sin(theta)*x + cos(theta)*z;
    
    lon_new = atan2(y_new, x_new) # Convert cartesian back to spherical coordinates
    lat_new = asin(z_new)
    
    lon_new = (lon_new*180)/pi
    lat_new = (lat_new*180)/pi
    
    lonlat_df = data.frame(cbind(lon_new, lat_new, 1:length(lon_new)))
    colnames(lonlat_df) = c("lon", "lat", "ind")
    
    return(lonlat_df)
  }
}

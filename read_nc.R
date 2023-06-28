# source("read_nc.R")

library(ncdf4)
library(fields)
library(maps)
library(mapdata)
library(PCICt)

flist  <- "/data/users/hadsx/model_data/cpm/halfhourly/leeds_msc/r001i1p00000_19991201-19991230_pr.nc"

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
    v1    <- ncvar_get(nc1,'precipitation_flux',start=c(1,1,1), count=c(-1,-1,48))
    v1att <- ncatt_get(nc1,'precipitation_flux')
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
n.map <- map('worldHires', xlim=c(-10,10),ylim=c(30,70),interior=F,plot=FALSE)
r.map <- n.map
for (i in 1:length(r.map$x)){
        if(!is.na(n.map$y[i]) & !is.na(n.map$x[i])) r.latlon   <- pp.ll.to.rg(n.map$y[i],n.map$x[i], gr_npole_lat, gr_npole_lon)
        r.map$x[i] <- r.latlon[2]
        r.map$y[i] <- r.latlon[1]
}


### plot a single field
image.plot( x=rlon, y=rlat, z=v1[,,1], zlim=c(0,20), main=pcictime[1])
points(r.map$x+360,r.map$y,pch=46,col='grey80', lwd=0.1)


### data structure
# > str(v1)
#  num [1:532, 1:654, 1:48] 0.172 0.164 0.151 0.133 0.116 ...
# dim1 is x/longitude
# dim2 is y/latitude
# dim3 is time

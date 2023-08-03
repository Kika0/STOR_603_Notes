# function to translate from normal pole coordinates ----
# to rotated pole coordinates
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

CnvRttPol <-  function(latlon, spol_coor, option=1){
  if(option==1){
    lon = latlon[,1]
    lat = latlon[,2]
    
    lon = (lon*pi)/180
    lat = (lat*pi)/180
    
    spol_lon = spol_coor[1]+180
    spol_lat = -spol_coor[2]
    
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
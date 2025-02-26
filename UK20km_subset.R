# libraries
library(ncdf4)
library(sf)
# load sf object of UK mainland polygon


### OBS 5km ----- please change file path to any /UKgrid5km/*.nc file
obs_example  <- '../luna/kristina/UKgrid5km/tasmax_rcp85_land-cpm_uk_5km_01_ann_206012-208011.nc'
nc1      <- nc_open(obs_example)
tas5.o   <- ncvar_get(nc1, "tasmax")
lon5.o   <- ncvar_get(nc1, "longitude")
lat5.o   <- ncvar_get(nc1, "latitude")
nc_close(nc1)

# NOTE: only need to work with files that overlap with mainland UK
# find a subset of x and y that overlap with mainland UK
# create a dataframe with also x (row) and y (column) indeces
# as.vector() does column by column
xy_df <- data.frame("lon"=as.vector(lon5.o),"lat"=as.vector(lat5.o),"x"= rep(1:dim(lon5.o)[1],n=dim(lon5.o)[2]), "y" = rep(1:dim(lon5.o)[2],each=dim(lon5.o)[1]), "temp"=as.vector(tas5.o[,,10]))
xy_sf <- xy_df %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POINT") %>% st_make_valid()
xy_sf <- cbind(xy_sf,xy_df)
# map to check
tmap_mode("view")
tm_shape(xy_sf) + tm_dots(col="temp")

# looks resonable, now subset over mainland UK
lad <- st_read("data/LAD_boundary_UK/LAD_MAY_2022_UK_BFE_V3.shp")
# UK is comprised of many polygons (islands), simplify to only
# take mainland UK (Great Britain)
uk <- (lad %>% st_union() %>% st_cast( "MULTIPOLYGON" ) %>% st_cast("POLYGON"))[1531]
uk <- st_simplify(uk,dTolerance = 2000) %>% st_transform(crs = 4326)  
save(uk,"uk.RData")

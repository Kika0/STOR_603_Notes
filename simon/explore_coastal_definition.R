library(tmap)
library(tidyverse)
library(sf)

# load("data_processed/UK_all_coastal_points.RData")
# summary(coor.coast)
# # plot
# ggplot(coor.coast) + geom_point(aes(x=ilon,y=ilat))

# load points
load("data_processed/spatial_helper.RData")

coastal_point <- function(grid) {
  sapply(1:nrow(grid),FUN = function(i) {sum(as.vector(st_distance(grid[i,],grid))<20500)<5})
}


cp <- coastal_point(grid = xyUK20_sf) 
t1 <- tm_shape(cbind(xyUK20_sf,data.frame(cp))) + tm_dots(fill="cp",size=0.5)+ tm_layout(legend.position=c("right","top"),legend.height = 12)

coast_within <- function(grid,uk_polygon,x) {
  # add buffer for calculating distance from the coast (distance from within the polygon sometimes causes error)
  uk_buffered <- st_buffer(uk_polygon, 50000)   
  # polygon of only the buffer
  buffer_only <- st_difference(uk_buffered, uk_polygon) 
  y <- sapply(1:nrow(grid),FUN = function(i) {as.vector(st_distance(grid[i,],buffer_only))})
  y1 <- rep(FALSE,nrow(grid))
  for (i in length(x):1) {
    y1[y<x[i]*1000] <- as.character(x[i])
  }
  return(y1)
}

cw <- coast_within(grid = xyUK20_sf,uk_polygon = uk,x=c(1,5,10,20))

t2 <- tm_shape(uk) + tm_polygons() + tm_shape(cbind(xyUK20_sf,data.frame(cw))) + tm_dots(fill="cw",size=0.5,fill.scale = tm_scale_categorical(values=c("20" = "#C11432","10" = "#66A64F","5"="#009ADA","1"="#FDD10A", "FALSE" = "black")))+ tm_layout(legend.position=c("right","top"),legend.height = 12)
t <- tmap_arrange(t1,t2)
tmap_save(t,filename = "../Documents/coastal_definition.png",width=6,height=6)

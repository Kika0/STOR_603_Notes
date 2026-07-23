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

# explore east coast definition
legend_text_size <- 0.7
point_size <- 0.3
legend_title_size <- 0.9
coastal_point <- function(grid) {
  sapply(1:nrow(grid),FUN = function(i) {sum(as.vector(st_distance(grid[i,],grid))<20500)<5 & grid$lon[i]>-1})
}

cp <- coastal_point(grid = xyUK20_sf) 
t1 <- tm_shape(cbind(xyUK20_sf,data.frame(cp))) + tm_dots(fill="cp",fill.scale = tm_scale_categorical(values=c("FALSE"="black","TRUE"="#C11432")),size=point_size, fill.legend = tm_legend(title="")) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,legend.show=FALSE,frame=FALSE) + tm_title(text="Longitude > -1")

# try as a function of lon and lat
coastal_point <- function(grid) {
  sapply(1:nrow(grid),FUN = function(i) {sum(as.vector(st_distance(grid[i,],grid))<20500)<5 & grid$lon[i]>-2})
}

cp <- coastal_point(grid = xyUK20_sf) 
t2 <- tm_shape(cbind(xyUK20_sf,data.frame(cp))) + tm_dots(fill="cp",fill.scale = tm_scale_categorical(values=c("FALSE"="black","TRUE"="#C11432")),size=point_size, fill.legend = tm_legend(title="")) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,legend.show=FALSE,frame=FALSE) + tm_title(text="Longitude > -2")

# try as a function of longitude and latitude
coastal_point <- function(grid) {
  sapply(1:nrow(grid),FUN = function(i) {sum(as.vector(st_distance(grid[i,],grid))<20500)<5 & grid$lat[i]+2*grid$lon[i]>48.5})
}

cp <- coastal_point(grid = xyUK20_sf) 
cp[df_sites[3,5]] <- FALSE
t3<- tm_shape(cbind(xyUK20_sf,data.frame(cp))) + tm_dots(fill="cp",fill.scale = tm_scale_categorical(values=c("FALSE"="black","TRUE"="#C11432")),size=point_size, fill.legend = tm_legend(title="")) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,legend.show=FALSE,frame=FALSE) + tm_title(text="East Coast?")

# save as one picture
tmap_save(tmap_arrange(t1,t2,t3,ncol=3),filename="../Documents/east_coast.png",height=6,width=9)

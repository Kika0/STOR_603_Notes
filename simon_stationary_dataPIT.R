library(tidyverse)
library(ncdf4)
library(fields)
library(sf)
library(tmap)
load("data_processed/spatial_helper.RData")
# load the data
load("../luna/kristina/MSdata_CPM5km_member1/ukgd_cpm85_5k_x100y23_MSdata01.RData")
data01.std.param # not much info here
data01 <- data01 %>% mutate(year=floor(time))
# NA in uqgam may mean high value of u
data01 %>% filter(is.na(uqgam)) %>% view()
summary(data01 %>% filter(is.na(uqgam)) )
# is it right to replace these with low temperature? maybe replace with an average of neighbouring temperatures?

# when u is NA, uqgam is also NA
data01 %>% filter(is.na(u)) %>% view()
summary(data01 %>% filter(is.na(u)))

dim(data01)
summary(data01)
glimpse(data01)
data01 %>% group_by(year) %>% summarize(mean_year=mean(x)) %>% ggplot(aes(x=year,y=mean_year)) + geom_point()
ggplot(data01 %>% filter(doy %in% c(1:20))) + geom_point(aes(x=time,y=x))
# to explore gpd functions for the tails
str(gpdpar)
summary(chosen.MSgpd)

# load the Gpd parameters
load("../luna/kristina/MSGpdParam/ukgd_cpm85_5k_x100y23.MSGpdParam.2024-10-22-064747.RData")
glimpse(gpdpar) # dataframe for columns for scale, shape and threshold

###############################################################################
# 5km
#
# NOTE the obs are on a bigger grid than the CPM grid
###############################################################################
obs_cpm_offset_y <- 33  # the difference in size of the y dimension

xcoord_m         <- 146 # London x146 y44 on the CPM grid from Laura's city file.
ycoord_m         <- 44
xcoord_o         <- xcoord_m
ycoord_o         <- ycoord_m + obs_cpm_offset_y

### OBS 5km
obs_example  <- '../luna/kristina/UKgrid5km/tasmax_rcp85_land-cpm_uk_5km_01_ann_206012-208011.nc'
nc1      <- nc_open(obs_example)
vlist    <- nc1$var
shape.o  <- vlist$tasmax$size
tas5.o   <- ncvar_get(nc1, "tasmax")
lon5.o   <- ncvar_get(nc1, "longitude")
lat5.o   <- ncvar_get(nc1, "latitude")
nc_close(nc1)
par(mfcol=c(1,2))

quilt.plot(as.vector(lon5.o), as.vector(lat5.o), tas5.o[,,10],  nx=shape.o[1], ny=shape.o[2], main='OBS 5km')

image.plot(tas5.o[,,10],  main='OBS 5km')

#  str(tas5.o) num [1:180, 1:290, 1:31] NA NA NA NA NA NA NA NA NA NA ...
#  str(lon5.o) num [1:180, 1:290] -9.99 -9.92 -9.86 -9.79 -9.72 ...
#  str(lat5.o) num [1:180, 1:290] 47.8 47.8 47.9 47.9 47.9 ...
#
#  str(tas5.m) num [1:180, 1:244, 1:3600] 10.1 10.2 10.3 10.3 10.2 ...
#  str(lon5.m) num [1:180, 1:244] -10.23 -10.16 -10.09 -10.02 -9.95 ...
#  str(lat5.m) num [1:180, 1:244] 49.3 49.3 49.3 49.3 49.3 ...

#save(file="grid-info-5km.RData", lon5.o, lat5.o, lon5.m, lat5.m, obs_example, cpm2k_example)

# NOTE: only need to work with files that overlap with mainland UK
# find a subset of x and y that overlap with mainland UK
summary(as.vector(lon5.o))
summary(as.vector(lat5.o))
# create a dataframe with also x (row) and y (column) indeces
# as.vector does column by column
xy_df <- data.frame("lon"=as.vector(lon5.o),"lat"=as.vector(lat5.o),"x"= rep(1:dim(lon5.o)[1],n=dim(lon5.o)[2]), "y" = rep(1:dim(lon5.o)[2],each=dim(lon5.o)[1]), "temp"=as.vector(tas5.o[,,10]))
xy_sf <- xy_df %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POINT") %>% st_make_valid()
xy_sf <- cbind(xy_sf,xy_df)
# map
tmap_mode("view")
tm_shape(xy_sf) + tm_dots(col="temp")
# looks resonable, now subset over mainland UK
lad <- st_read("data/LAD_boundary_UK/LAD_MAY_2022_UK_BFE_V3.shp")
# UK is comprised of many polygons (islands), simplify to only
# take mainland UK (Great Britain)
uk_notsimplified <- (lad %>% st_union() %>% st_cast( "MULTIPOLYGON" ) %>% st_cast("POLYGON"))[1531]
uk <- st_simplify(uk_notsimplified,dTolerance = 2000) %>% st_transform(crs = 4326)  
# check
# tm_shape(uk_notsimplified) + tm_polygons()
# great, now subset
xyUK_sf <- st_filter(xy_sf,uk)
tm_shape(xyUK_sf) + tm_dots(col="temp")
# take only every fourth dot in each x and y
x20 <- seq(from=4,to=dim(lon5.o)[1],by=4)
y20 <- seq(from=4,to=dim(lon5.o)[2],by=4)
xyUK20_sf <-xyUK_sf %>% filter(x %in% x20,y %in% y20)
# save sf objects used for spatial analysis
save(uk,uk_notsimplified,xyUK20_sf,files_subset1,file="data_processed/spatial_helper.RData")

tm_shape(xyUK20_sf) + tm_dots(col="temp")
# great, now load only files that overlap this grid or perhaps delete all other files?
files <- list.files("../luna/kristina/MSdata01/")
list_of_files <- list() #create empty list
# only subset for x in x20 and y in y20
files_subset <- sapply(1:nrow(xyUK20_sf),function(i){paste0("ukgd_cpm85_5k_x",xyUK20_sf$x[i],"y",xyUK20_sf$y[i],"_MSdata01.RData")})
# explore an example of a missing file
head(sort(files_subset),n=50)
head(sort(files_subset1),n=50)

# what is the one missing files?
files_missing <- files_subset[!(files_subset %in% files)]
files_missing[1] # it is over Heysham in the sea, a mistake due to a simplified UK polygon

#loop through the files
files_subset1 <- files_subset[files_subset %in% files]
# could take only x and y divisible by 4 to subset and speed up data loading
for (i in 1:length(files_subset1)) {
  print(files_subset1[i])
  load(paste0("../luna/kristina/MSdata01/", files_subset1[i]))
  list_of_files[[i]] <- data01 #add files to list position
}
xyUK20_sf <- xyUK20_sf[files_subset %in% files_subset1,]

# join the summer data together could try 1960-1999 for non-stationary data 
# June 1 is 152 doy, August 31 is 243 doy (92 days per year)
# create a dataframe
xyUK20 <- xyUK20_sf %>% dplyr::select(-temp) %>% st_drop_geometry() %>% cbind(as.data.frame(matrix(data=numeric(),ncol=92*40,nrow=nrow(xyUK20_sf))))
names(xyUK20)[5:ncol(xyUK20)] <- paste0(rep(152:243,40),"_",rep(1960:1999,each=92)) 

for (i in 1: length(list_of_files)) {
  xyUK20[i,5:ncol(xyUK20)] <- list_of_files[[i]] %>% mutate(year=floor(time)) %>% filter(class=="obs",doy>=152,doy<=243, year<=1999) %>% pull(x) 
}
# link back to spatial
tm_shape(cbind(xyUK20_sf,xyUK20)) + tm_dots(col="X155_1990")
# save as R object for further analysis

# explore NA values
# p <- data01 %>% mutate(year=floor(time)) %>% filter(is.na(uqgam)) %>% group_by(year,class) %>% summarize(n=n()) %>% arrange(-n) %>% ggplot() + geom_line(aes(x=year,y=n,col=class)) + ggtitle("Counts of NA values of uqgam for each year for observation and model data")
# ggsave(p,filename="../Documents/newdata/uqgamNA.png",width=10,height=5)

# setup for par_est
sims <- xyUK20 %>% dplyr::select(-all_of(1:4)) %>% t() %>% as.data.frame()
# ordered alphabetically so Y1 Birmingham, Y2 Glasgow and Y3 is London
colnames(sims) <- paste0("Y",1:ncol(sims))
# transform to Laplace margins
sims <- as.data.frame((sims %>% apply(c(2),FUN=row_number))/(nrow(sims)+1)) %>% apply(c(1,2),FUN=unif_laplace_pit) %>% as.data.frame()
v <- 0.9 # set threshold

# create a dataframe for stationary data
# create a dataframe
xyUK20 <- xyUK20_sf %>% dplyr::select(-temp) %>% st_drop_geometry() %>% cbind(as.data.frame(matrix(data=numeric(),ncol=90*100,nrow=nrow(xyUK20_sf))))
names(xyUK20)[5:ncol(xyUK20)] <- paste0(rep(91:180,40),"_",rep(1981:2080,each=90)) 

for (i in 1: length(list_of_files)) {
  xyUK20[i,5:ncol(xyUK20)] <- list_of_files[[i]] %>% mutate(year=floor(time)) %>% filter(class=="mod") %>% pull(x) 
}
# link back to spatial
tm_shape(cbind(xyUK20_sf,xyUK20)) + tm_dots(col="X155_1990")
# save as R object for further analysis

# explore NA values
# p <- data01 %>% mutate(year=floor(time)) %>% filter(is.na(uqgam)) %>% group_by(year,class) %>% summarize(n=n()) %>% arrange(-n) %>% ggplot() + geom_line(aes(x=year,y=n,col=class)) + ggtitle("Counts of NA values of uqgam for each year for observation and model data")
# ggsave(p,filename="../Documents/newdata/uqgamNA.png",width=10,height=5)

# setup for par_est
sims <- xyUK20 %>% dplyr::select(-all_of(1:4)) %>% t() %>% as.data.frame()
# ordered alphabetically so Y1 Birmingham, Y2 Glasgow and Y3 is London
colnames(sims) <- paste0("Y",1:ncol(sims))
# transform to Laplace margins
sims <- as.data.frame((sims %>% apply(c(2),FUN=row_number))/(nrow(sims)+1)) %>% apply(c(1,2),FUN=unif_laplace_pit) %>% as.data.frame()
v <- 0.9 # set threshold

sum(is.na(sims))
# run for Birmingham and tau=0 to check
cond_site <- find_site_index(as.numeric(df_sites[,1]),grid_uk = xyUK20_sf)
pe <- par_est(df=sims,v=v,given=cond_site,margin = "AGG", method="sequentialGG",keef_constraints = c(1,2))
i <- 1
pet <- pe %>% add_row(.before=cond_site) %>%  mutate(cond_site=names(df_sites[i]),tau=as.character(0))

# plot parameters
petsf <- est_join_spatial(tmp_est=pet,grid_uk=xyUK20_sf)
p1 <- map_param(tmp_est=petsf %>% filter(cond_site==names(df_sites)[i]),method="AGG",facet_var = "cond_site",title_map = "With a bug",grid_uk = xyUK20_sf)

# estimate residual margins
# compare parameter estimates with AGG fitted to the observed residuals
obsr <- observed_residuals(df=sims,given=cond_site,v = v,a=pe$a,b=pe$b)
# fit AGG for each column
dum1 <- data.frame("mu_agg"=numeric(),"sigl"=numeric(),"sigu" = numeric(),"deltal"=numeric(),"deltau" = numeric())
opt2 <- list()
for (i in 1:ncol(obsr)) {
  Z2 <- as.numeric(unlist(obsr[,i]))
  opt2[[i]] <- optim(fn=NLL_AGG,x=Z2,par=c(mean(Z2),sd(Z2),sd(Z2),1.2,1.8),control=list(maxit=2000),method = "Nelder-Mead")
  dum1[nrow(dum1)+1,] <- optim(fn=NLL_AGG,x=Z2,par=c(mean(Z2),sd(Z2),sd(Z2),1.2,1.8),control=list(maxit=2000),method = "Nelder-Mead")$par
}

# compare estimates
plot(pe$mu_agg,dum1$mu_agg) # different as expected

# compare estimates on a map
i <- 1
pet1 <- cbind(pe[,c(5:7,9,13,16,17)],dum1) %>% add_row(.before=cond_site) %>%  mutate(cond_site=names(df_sites)[i])
petsf1 <- est_join_spatial(tmp_est=pet1,grid_uk=xyUK20_sf)
p2 <- map_param(tmp_est=petsf1 %>% filter(cond_site==names(df_sites)[i]),method="AGG",facet_var = "cond_site",title_map = "Bug fixed",grid_uk = xyUK20_sf)
x <- c("mu_agg","sigl","sigu","sigdiff","deltal","deltau","deltadiff")
for (k in c(5:11)) {
tmap_save(tmap_arrange(p1[[k]],p2[[k]]),filename=paste0("../Documents/AGG_troubleshoot/Birmingham_",x[k-4],".png"),width=8,height=8)
}


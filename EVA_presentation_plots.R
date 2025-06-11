#library(tmap) # spatial map plots
library(gridExtra)
library(sf) # for handling spatial sf objects
library(tidyverse)
library(latex2exp) # latex expressions for plot labels
theme_set(theme_bw())
theme_replace(
  panel.spacing = unit(2, "lines"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA) )
file.sources = list.files(pattern="*helpers.R")
sapply(file.sources,source,.GlobalEnv)
source("spatial_parameter_estimation.R") # for spatial_par_est function
load("data_processed/temperature_data.RData")
load("data_processed/spatial_helper.RData")
Birmingham <- c(-1.9032,52.4806)
Glasgow <- c(-4.258109,55.859112)
London <- c(-0.127676,51.529972)
Inverness <- c(-4.22498,57.48065) # Inverness bus station
Lancaster <- c(-2.78440,54.00871) # PSC building Lancaster University
Newcastle <- c(-1.61682,54.96902) # Newcastle railway station
Cromer <- c(1.28486,53.05349)
Hull <- c(-0.335827,53.767750)
Lowestoft <- c(1.72431,52.48435)
Truro <- c(-5.05125342465549,50.263075821232704)
Dolgellau <- c(-3.8844362867080897,52.74213275545185)
Bournemouth <- c(-1.8650607066137428,50.72173094587856)
df_sites <- data.frame(Birmingham,Glasgow,London,Inverness,Lancaster,Newcastle,Cromer,Hull,Lowestoft,Truro,Dolgellau,Bournemouth)

# 1. plot of Birmingham and Glasgow on original and on Laplace margins
Birm_index <- find_site_index(site=Birmingham,grid_uk = xyUK20_sf)
Glasgow_index <- find_site_index(site=Glasgow,grid_uk = xyUK20_sf)

tmp <- data_mod[,c(Birm_index,Glasgow_index)]
names(tmp) <- names(df_sites)[1:2]
p1 <- ggplot(tmp) + geom_point(aes(x=Birmingham,y=Glasgow),size=0.5)
tmpL <- data_mod_Lap[,c(Birm_index,Glasgow_index)]
names(tmpL) <- names(df_sites)[1:2]
p2 <- ggplot(tmpL) + geom_point(aes(x=Birmingham,y=Glasgow),size=0.5)
p <- grid.arrange(p1,p2,ncol=2)
ggsave(p,filename= "../Documents/Birmingham_Glasgow_original_Laplace.png",width=8,height=4)

# 2. plot illustrating conditioning on a variable
London_index <- find_site_index(site=London,grid_uk = xyUK20_sf)
Inverness_index <- find_site_index(site=Inverness,grid_uk = xyUK20_sf)
tmp <- data_mod_Lap[,c(Birm_index,Glasgow_index,Inverness_index)]
names(tmp) <- names(df_sites)[c(1,3,4)]
p1 <- ggplot(tmp) + geom_point(aes(x=London,y=Birmingham),size=0.5) + xlab(TeX("$Y_1$")) + ylab(TeX("$Y_2$")) + scale_color_manual(values = c("black", "#009ADA")) + theme(legend.position="none") + coord_fixed() + xlim(-10,10) + ylim(-10,10)

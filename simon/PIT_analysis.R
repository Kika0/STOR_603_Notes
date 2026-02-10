#library(tmap)
#library(sf)
library(tidyverse)
library(latex2exp)
library(gridExtra)
library(LaplacesDemon)

theme_set(theme_bw())
theme_replace(
  panel.spacing = unit(2, "lines"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA) )

file.sources = list.files(pattern="*helpers.R")
sapply(file.sources,source,.GlobalEnv)
#source("spatial_parameter_estimation.R") # for spatial_par_est function
load("data_processed/temperature_data.RData",verbose=TRUE) # for data_mod_Lap
load("data_processed/spatial_helper.RData",verbose = TRUE) # for xyUK20_sf
q <- 0.9

# check maxima for each site
hist(apply(data_mod_Lap,MARGIN=c(2),max))
dev.print(pdf, '../Documents/histogram_site_maxima_laplace.pdf')
dev.off()
apply(data_mod_Lap,MARGIN=c(2),max)[c(192,260,321)]

# transform onto new Laplace margins


# estimate alpha and beta

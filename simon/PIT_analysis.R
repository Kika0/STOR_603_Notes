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
data_mod_Lap_star <- as.data.frame((data_mod_Lap %>% apply(c(2),FUN=row_number))/(nrow(data_mod_Lap)+1)) %>% apply(c(1,2),FUN=unif_laplace_pit) %>% as.data.frame()
# check a couple of sites
head(data_mod_Lap[,1:5])
head(data_mod_Lap_star[,1:5])
# estimate alpha and beta
pe <- par_est(df=data_mod_Lap_star,v=q,given=192,keef_constraints = c(1,2),margin="Normal",method="sequential2")
# compare with existing estimates
summary(pe$b)
summary(est_all_sf %>% dplyr::filter(cond_site=="Birmingham") %>% pull(b))

a <- est_all_sf %>% dplyr::filter(cond_site=="Birmingham") %>% pull(a)
b <- est_all_sf %>% dplyr::filter(cond_site=="Birmingham") %>% pull(b)
a_star <- pe$a
b_star <- pe$b

# histogram comparison
tmpa <- rbind(data.frame("alpha"=a,"method"="original_estimate"),data.frame("alpha"=a_star,"method"="new_estimate"))
tmpb <- rbind(data.frame("beta"=b,"method"="original_estimate"),data.frame("beta"=b_star,"method"="new_estimate"))
p1 <- ggplot(tmpa) + geom_density(aes(x=alpha,fill=method),alpha=0.5) + scale_fill_manual(values=c("original_estimate"="black","new_estimate"="#C11432"))

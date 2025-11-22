library(tmap)
library(sf)
library(tidyverse)
library(latex2exp)
library(gridExtra)

theme_set(theme_bw())
theme_replace(
  panel.spacing = unit(2, "lines"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA) )

load("data_processed/spatial_helper.RData") # xyUK20_sf 20km grid
file.sources = list.files(pattern="*helpers.R")
sapply(file.sources,source,.GlobalEnv)

# set folder name to save plots to
folder_name <- "../Documents/GLEN_workshop_plots/"

# 1. plot alphas temporally and spatially ------------------------------
# load estimates

# 2. plot of sigma values comparison -----------------------------------
# load estimates
load("data_processed/iterative_abmu_fixed_delta.RData")
c12 <- c(
  "#009ADA", "#C11432", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "#FB9A99", # lt pink
  "gray70", 
  "darkturquoise", "green1", 
  "darkorange4"
)
# identify diagonal sites
# find indeces of start and end sites
Birmingham <- c(-1.9032,52.4806)
Cromer <- c(1.28486,53.05349)
site_start <- find_site_index(site=Birmingham,grid_uk = xyUK20_sf)
site_end <- find_site_index(site=Cromer,grid_uk = xyUK20_sf)
sites_index_diagonal <- c(192,193,194,195,196,197,218,219,220,221,242,263) # first is Birmingham and last is Cromer
site_name_diagonal <- c("Birmingham", paste0("diagonal",1:(length(sites_index_diagonal)-2)),"Cromer")

# plot sigmas against distance
get_sigma_distance <- function(i, grid = xyUK20_sf,site_names=site_name_diagonal,tmp = tmp_fixed_deltas,sig_which="sigu", cond_site_index = sites_index_diagonal) {
  sig <- tmp[[i]][[12]] %>% dplyr::select(contains(sig_which)) %>% pull()
  sig <- append(sig,NA,after=(cond_site_index[i]-1))
  dist_cond_site <- as.numeric(unlist(st_distance(grid[cond_site_index[i],],grid) %>%
                                        units::set_units(mi)))
  sigd <- data.frame(sig=sig,dist=dist_cond_site) 
  return(sigd %>% mutate(cond_site = site_name_diagonal[i]))
}

lims <- c(0,2.3)
point_size <- 0.7
tmp_sigl <- do.call(rbind,lapply(1:length(site_name_diagonal),FUN=get_sigma_distance,sig_which="sigl")) %>% mutate(cond_site=factor(cond_site,levels=site_name_diagonal))
p1 <- ggplot(tmp_sigl) + geom_point(aes(y=sig,x=dist,col=cond_site),size=point_size) + scale_color_manual(values = c12) + xlab("Distance [mi]") + ylab(TeX("$\\sigma_l$")) + theme(legend.position="none") + ylim(lims)
tmp_sigu <- do.call(rbind,lapply(1:length(site_name_diagonal),FUN=get_sigma_distance,sig_which="sigu")) %>% mutate(cond_site=factor(cond_site,levels=site_name_diagonal))
p2 <- ggplot(tmp_sigu) + geom_point(aes(y=sig,x=dist,col=cond_site),size=point_size) + scale_color_manual(values = c12) + xlab("Distance [mi]") + ylab(TeX("$\\sigma_u$")) + theme(legend.position="none") + ylim(lims)
tmp_sigmas <- data.frame("sigl"=tmp_sigl$sig, "sigu" = tmp_sigu$sig, "cond_site" = tmp_sigl$cond_site)
p3 <- ggplot(tmp_sigmas) + geom_point(aes(y=sigu,x=sigl,col=cond_site),size=point_size) + scale_color_manual(values = c12) + xlab(TeX("$\\sigma_l$")) + ylab(TeX("$\\sigma_u$")) + xlim(lims) + ylim(lims) + coord_fixed() + labs(col="Conditioning site") + guides(col=guide_legend(override.aes=list(size=4))) + geom_abline(slope=1,linetype="dashed")
p <- grid.arrange(p1,p2,p3,ncol=3)
# save
ggsave(p,filename=paste0(folder_name,"sigmas_distance.png"),width=15,height=3.5)
# remove
rm(tmp_sigl,tmp_sigu,tmp_sigmas,lims,point_size,p1,p2,p3,p)

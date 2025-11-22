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

# 1. plot alphas for a couple of sites ------------------------------
# load estimates
q <- 0.9
load(paste0("data_processed/N9000_sequential2_AGG_all12sites",q*100,".RData"))
#est_all <- as.data.frame(est_all_sf)
# plot alpha values for Lancaster, Birmingham, Cromer
cond_site_names <- c("Lancaster","Birmingham","Cromer")
est_sites <- est_all_sf %>% filter(cond_site %in% cond_site_names) %>% mutate(cond_site=factor(cond_site,levels = cond_site_names))
title_map <- ""
misscol <- "aquamarine"
legend_text_size <- 0.7
point_size <- 0.6
legend_title_size <- 1.2
lims <- c(0,1)
nrow_facet <- 1
p <- tm_shape(est_sites) + tm_dots(fill="a",fill.scale = tm_scale_continuous(limits=lims,values="brewer.blues",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\alpha$"))) + tm_facets(by="cond_site",nrow = nrow_facet) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE) + tm_title(text=title_map) 
tmap_save(p,filename=paste0(folder_name,"alpha_selected_sites.png"))
lims <- c(0,0.5)
p <- tm_shape(est_sites) + tm_dots(fill="b",fill.scale = tm_scale_continuous(limits=lims,values="brewer.blues",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\beta$"))) + tm_facets(by="cond_site",nrow = nrow_facet) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE) + tm_title(text=title_map) 
tmap_save(p,filename=paste0(folder_name,"beta_selected_sites.png"))
# remove objects
rm(est_sites,p)

# calculate estimates from Birmingham to London


# 3. plot of sigma values comparison -----------------------------------
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

library(tmap)
library(sf)
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

a <- est_all_sf %>% dplyr::filter(cond_site=="Birmingham") %>% pull(a) %>% na.omit()
b <- est_all_sf %>% dplyr::filter(cond_site=="Birmingham") %>% pull(b) %>% na.omit()
a_star <- pe$a
b_star <- pe$b

# density comparison
tmpa <- rbind(data.frame("alpha"=a,"method"="original_estimate"),data.frame("alpha"=a_star,"method"="new_estimate"))
tmpb <- rbind(data.frame("beta"=b,"method"="original_estimate"),data.frame("beta"=b_star,"method"="new_estimate"))
p1 <- ggplot(tmpa) + geom_density(aes(x=alpha,fill=method),alpha=0.5) + scale_fill_manual(values=c("original_estimate"="black","new_estimate"="#C11432"))
p2 <- ggplot(tmpb) + geom_density(aes(x=beta,fill=method),alpha=0.5) + scale_fill_manual(values=c("original_estimate"="black","new_estimate"="#C11432"))
p <- grid.arrange(p1,p2,ncol=2)
ggsave(p,filename="../Documents/ab_compare_PIT.png",width=10,height=5)

# map these estimates
tmpab <- cbind(data.frame(a,a_star,b,b_star) %>% add_row(.before = 192),xyUK20_sf)
tmpsf <- st_as_sf(tmpab %>% mutate(adiff = a_star-a, bdiff = b_star-b))
toplabel <- c("After transformation","Original method","Difference")

misscol <- "aquamarine"
  #mu_limits <- c(-2.61,1)
t1 <- tmpsf %>% dplyr::select(a_star,a,adiff) %>% pivot_longer(cols=c(a_star,a,adiff),names_to = "parameter", values_to = "value" ) %>% mutate(parameter=factor(parameter,levels=c("a_star","a","adiff")) ) %>% dplyr::filter(parameter=="adiff") %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="-brewer.rd_bu",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\alpha$"))) + tm_facets("parameter") +
  tm_layout(legend.position=c("right","top"),legend.height = 12, panel.labels = toplabel[3],legend.reverse = TRUE) 
t2 <- tmpsf %>% dplyr::select(b_star,b,bdiff) %>% pivot_longer(cols=c(b_star,b,bdiff),names_to = "parameter", values_to = "value" ) %>% mutate(parameter=factor(parameter,levels=c("b_star","b","bdiff")) ) %>% dplyr::filter(parameter=="bdiff") %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="-brewer.rd_bu",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\beta$"))) + tm_facets("parameter") +
  tm_layout(legend.position=c("right","top"),legend.height = 12, panel.labels = toplabel[3],legend.reverse = TRUE) 
t <- tmap_arrange(t1,t2,ncol=2)
tmap_save(t,filename=paste0("../Documents/ab_compare_PIT_map.png"),width=8,height=8)
# plot also original values
t1 <- tmpsf %>% dplyr::select(a_star,a,adiff) %>% pivot_longer(cols=c(a_star,a,adiff),names_to = "parameter", values_to = "value" ) %>% mutate(parameter=factor(parameter,levels=c("a_star","a","adiff")) ) %>% dplyr::filter(parameter %in% c("a_star","a")) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="-brewer.rd_bu",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\alpha$"))) + tm_facets("parameter") +
  tm_layout(legend.position=c("right","top"),legend.height = 12, panel.labels = toplabel[1:2],legend.reverse = TRUE) 
t2 <- tmpsf %>% dplyr::select(b_star,b,bdiff) %>% pivot_longer(cols=c(b_star,b,bdiff),names_to = "parameter", values_to = "value" ) %>% mutate(parameter=factor(parameter,levels=c("b_star","b","bdiff")) ) %>% dplyr::filter(parameter %in% c("b_star","b")) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="-brewer.rd_bu",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\beta$"))) + tm_facets("parameter") +
  tm_layout(legend.position=c("right","top"),legend.height = 12, panel.labels = toplabel[1:2],legend.reverse = TRUE) 
t <- tmap_arrange(t1,t2,ncol=2)
tmap_save(t,filename=paste0("../Documents/ab_compare_PIT_map_values.png"),width=12,height=6)

# plot observed residuals
Z <- observed_residuals(df=data_mod_Lap,given=192,v = q,a=a,b=b)
Z_star <- observed_residuals(df=data_mod_Lap_star,given=192,v = q,a=a_star,b=b_star)

# calculate mean
Zmean <- apply(Z,MARGIN=c(2),FUN=mean)
Zmean_star <- apply(Z_star,MARGIN=c(2),FUN=mean)
# calculate variance
Zvar <- apply(Z,MARGIN=c(2),FUN=var)
Zvar_star <- apply(Z_star,MARGIN=c(2),FUN=var)
tmpmean <- rbind(data.frame("Zmean"=Zmean,"method"="original_estimate"),data.frame("Zmean"=Zmean_star,"method"="new_estimate"))
tmpvar <- rbind(data.frame("Zvar"=Zvar,"method"="original_estimate"),data.frame("Zvar"=Zvar_star,"method"="new_estimate"))
p1 <- ggplot(tmpmean) + geom_density(aes(x=Zmean,fill=method),alpha=0.5) + scale_fill_manual(values=c("original_estimate"="black","new_estimate"="#C11432"))
p2 <- ggplot(tmpvar) + geom_density(aes(x=Zvar,fill=method),alpha=0.5) + scale_fill_manual(values=c("original_estimate"="black","new_estimate"="#C11432"))
p <- grid.arrange(p1,p2,ncol=2)
ggsave(p,filename="../Documents/Z_mean_var_compare_PIT.png",width=10,height=5)

# plot selected sites
p260 <- ggplot() + geom_density(Z,mapping=aes(x=Z260),fill="black",alpha=0.5) + geom_density(Z_star,mapping=aes(x=Z260),fill="#C11432",alpha=0.5)
p321 <- ggplot() + geom_density(Z,mapping=aes(x=Z321),fill="black",alpha=0.5) + geom_density(Z_star,mapping=aes(x=Z321),fill="#C11432",alpha=0.5)
p <- grid.arrange(p260,p321,ncol=2)
ggsave(p,filename="../Documents/Z_selected_sites_compare_PIT.png",width=10,height=5)


library(tmap) # spatial map plots
library(sf) # for handling spatial sf objects
library(viridis)
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
#Leeds <- c(-1.5410242288355958,53.80098118214994)
df_sites <- data.frame(Birmingham,Glasgow,London,Inverness,Lancaster,Newcastle,Cromer,Hull,Lowestoft,Truro,Dolgellau,Bournemouth)
#spatial_par_est saves parameter estimates as est_all_sf sf object in ../Documents folder
q <- 0.9 # quantile threshold
# load all three parameter estimates sf objects
load(paste0("data_processed/N9000_sequential2_AGG_all12sites",q*100,".RData"))
est_all <- as.data.frame(est_all_sf)

# calculate the observed residuals (set one of 12 sites as conditioning site)
# cond_site_name <- "Newcastle"
# cond_site_coord <- df_sites %>% dplyr::select(all_of(cond_site_name)) %>% pull()
# cond_site <- find_site_index(cond_site_coord,grid_uk = grid)
# Z <- observed_residuals(df = data_mod_Lap, given = cond_site, v = q,a= discard(as.numeric(est_all$a),is.na),b = discard(as.numeric(est_all$b),is.na))
# 
# # estimate parameters iteratively --------------------------------------------
# Nite <- 20
# tmp1 <- res_margin_par_est_ite(df=Z,show_ite = TRUE,N=Nite)
# 
# # examine output
# 
# # plot delta estimates
# tmp_delta <- rbind(data.frame(delta=as.numeric(unlist(tmp1[[4]][1,])),iteration=1:(Nite+1),parameter = "delta_lower"),
#                    data.frame(delta=as.numeric(unlist(tmp1[[5]][1,])),iteration=1:(Nite+1),parameter = "delta_upper"))
# pd <- ggplot(tmp_delta) + geom_point(aes(x=iteration,y=delta,col=parameter))
# ggsave(pd,filename=paste0("../Documents/iterative_deltas_res_margin/deltas_",cond_site_name,".png"))
# 
# # plot across iterations for a selected site
# random_site <- 100
# tmp_all <- rbind(data.frame(delta=as.numeric(unlist(tmp1[[1]][random_site,])),iteration=1:(Nite+1),parameter = "mu_agg"),
#                    data.frame(delta=as.numeric(unlist(tmp1[[2]][random_site,])),iteration=1:(Nite+1),parameter = "sigma_lower"),
#                    data.frame(delta=as.numeric(unlist(tmp1[[3]][random_site,])),iteration=1:(Nite+1),parameter = "sigma_upper"),
#                    data.frame(delta=as.numeric(unlist(tmp1[[4]][random_site,])),iteration=1:(Nite+1),parameter = "delta_lower"),
#                    data.frame(delta=as.numeric(unlist(tmp1[[5]][random_site,])),iteration=1:(Nite+1),parameter = "delta_upper")
#                    
#                    )
# # plot the final iteration
# p1 <- ggplot(tmp_all) + geom_point(aes(x=iteration,y=delta,col=parameter))
# ggsave(p1,filename=paste0("../Documents/iterative_deltas_res_margin/random_site_par_",cond_site_name,".png"))
# # map final iteration or several
# est_site <- est_all_sf %>% filter(cond_site==cond_site_name)
# est_ite <- tmp[[6]] %>% add_row(.before=find_site_index(cond_site_coord,grid_uk = xyUK20_sf))
# names(est_ite) <- paste0(names(est_ite),"_ite")
# tmpsf <- cbind(est_site,est_ite)
# t <- tmpsf %>% dplyr::select(mu_agg,mu_agg_ite) %>% pivot_longer(cols=c(mu_agg,mu_agg_ite),names_to = "parameter", values_to = "value" ) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="-rd_bu")) + tm_facets("parameter") +
#   tm_layout(legend.position=c("right","top"),legend.height = 12) 
# 
# tmap_save(t,filename=paste0("../Documents/iterative_deltas_res_margin/mu_agg_",cond_site_name,".png"),width=8,height=6)
# 
# t <- tmpsf %>% dplyr::select(sigl,sigl_ite) %>% pivot_longer(cols=c(sigl,sigl_ite),names_to = "parameter", values_to = "value" ) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="brewer.blues")) + tm_facets("parameter") +
#   tm_layout(legend.position=c("right","top"),legend.height = 12) 
# tmap_save(t,filename=paste0("../Documents/iterative_deltas_res_margin/sigma_lower_",cond_site_name,".png"),width=8,height=6)
# 
# t <- tmpsf %>% dplyr::select(sigu,sigu_ite) %>% pivot_longer(cols=c(sigu,sigu_ite),names_to = "parameter", values_to = "value" ) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="brewer.blues")) + tm_facets("parameter") +
#   tm_layout(legend.position=c("right","top"),legend.height = 12) 
# tmap_save(t,filename=paste0("../Documents/iterative_deltas_res_margin/sigma_upper_",cond_site_name,".png"),width=8,height=6)

iter_delta_site_wrapper <- function(i) {
  iter_delta_site(j=i)
}
sapply(10,iter_delta_site,simplify=FALSE)
result <- sapply(1:ncol(df_sites),iter_delta_site,Nite=50,simplify = FALSE)
library(parallel)
result <- parallel::mclapply(1:ncol(df_sites),iter_delta_site,Nite=50,mc.cores=3)

deltal <- deltau <-  c()
for (i in 1:12) {
  deltal[i] <- result[[i]]$deltal_ite[1]
  deltau[i] <- result[[i]]$deltau_ite[1]
}

deltaldf <- cbind(data.frame(value=deltal),as.data.frame(t(df_sites))) %>% mutate(parameter="deltal")
deltaudf <- cbind(data.frame(value=deltau),as.data.frame(t(df_sites))) %>% mutate(parameter="deltau")
deltasf <- st_as_sf(rbind(deltaldf,deltaudf),coords =c(2:3))
# plot deltas on a spatial scale
st_crs(deltasf) <- st_crs(est_all_sf)
t <- tm_shape(est_all_sf) + tm_dots(size=0.2,fill_alpha=0.3) +  tm_shape(deltasf) + tm_dots(fill="value",size=1) + tm_facets(by="parameter") +
tm_layout(legend.position=c("right","top"),legend.height = 12)
tmap_save(t,filename=paste0("../Documents/all_deltas.png"),width=8,height=6)

# explore points above and below a line -----------------------
is_above <- function(x,y,x1,x2,y1,y2) { y>(y2-y1)/(x2-x1)*(x-x1)+y1
}
plot_is_above <- function(x1,x2,y1,y2,cond_var) {
  sites <- c("Birmingham", "Glasgow", "London")
  df <- data.frame(x=d[[cond_var]],y=x[[cond_var]],z=as.character(is_above(d[[cond_var]],x[[cond_var]],x1=x1,x2=x2,y1=y1,y2=y2)))
  p <- ggplot(df) + 
    geom_segment(x=x1,y=y1,xend=x2,yend=y2) +
    geom_point(aes(x=x,y=y,colour=factor(z))) + 
    ylab(TeX("$\\alpha$")) + xlab("Distance") + scale_color_manual(values = c("black", "#C11432")) + ggtitle(sites[cond_var])
  return(p)
}
sigu_above_below <- function(cond_site_name = "Birmingham",sites=df_sites,result.=result,x1,x2,y1,y2, grid = xyUK20_sf) {
  # filter conditioning site
  cond_site_coord <- sites %>% dplyr::select(all_of(cond_site_name)) %>% pull()
  cond_site <- find_site_index(cond_site_coord,grid_uk = grid)
 tmpsf <- result.[[ which(names(sites)==cond_site_name) ]]
  
  sigud <- data.frame(sigu=tmpsf$sigu_ite,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% mutate(is.above=is_above(x=dist,y=sigu,x1=x1,y1=y1,x2=x2,y2=y2))
 p <- ggplot(sigud) + 
   geom_segment(x=x1,y=y1,xend=x2,yend=y2) +
   geom_point(aes(x=dist,y=sigu,col=factor(is.above)),size=0.5) + 
   ylab(TeX("$\\sigma_u$")) + xlab("Distance") + scale_color_manual(values = c("black", "#C11432")) + ggtitle(cond_site_name) + guides(col="none")
   ggsave(p,filename=paste0("../Documents/iterative_deltas_res_margin/abovebelow_sigu_distance_",cond_site_name,".png"),width=4,height=4) 
  # plot also spatially
   sigusf <- cbind(tmpsf,sigud %>% select(dist,is.above))
   t <- tm_shape(sigusf) + tm_dots(fill="is.above",size=0.6,fill.scale = tm_scale_categorical(values=c("TRUE" = "#C11432", "FALSE" = "black")))  + tm_title(cond_site_name) + tm_layout(legend.position=c("right","top"),legend.height = 12) 
   tmap_save(t,filename=paste0("../Documents/iterative_deltas_res_margin/abovebelow_sigu_map_",cond_site_name,".png"),width=3,height=6) 
   return(sigud$is.above)
}

x1 <- rep(0,12)
x2 <- c(rep(600000,8),800000,900000,600000,600000)
y1 <- c(0.5,1,0.5,0.8,1,1,1.2,1,1.35,1.1,0.9,1.5)
y2 <- c(3,2.2,2.5,1.5,2.2,1.5,1.1,1.05,1.25,1.8,2.5,2)
# look at all sites
sapply(1:ncol(df_sites),FUN = function(i) {sigu_above_below(cond_site_name = names(df_sites)[i], x1 = x1[i], x2 = x2[i], y1 = y1[i], y2 = y2[i])})


# repeat for beta ------------------------------------------------
beta_above_below <- function(cond_site_name = "Birmingham",sites=df_sites,est_all_sf.=est_all_sf,x1,x2,y1,y2,grid = xyUK20_sf) {
  # filter conditioning site
  cond_site_coord <- sites %>% dplyr::select(all_of(cond_site_name)) %>% pull()
  cond_site <- find_site_index(cond_site_coord,grid_uk = grid)
  tmpsf <- est_all_sf. %>% filter(cond_site==cond_site_name)
  
  betad <- data.frame(b=tmpsf$b,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% mutate(is.above=is_above(x=dist,y=b,x1=x1,y1=y1,x2=x2,y2=y2))
  p <- ggplot(betad) + 
    geom_segment(x=x1,y=y1,xend=x2,yend=y2) +
    geom_point(aes(x=dist,y=b,col=factor(is.above)),size=0.5) + 
    ylab(TeX("$\\beta$")) + xlab("Distance") + scale_color_manual(values = c("black", "#C11432")) + ggtitle(cond_site_name) + guides(col="none")
  ggsave(p,filename=paste0("../Documents/iterative_deltas_res_margin/abovebelow_beta_distance_",cond_site_name,".png"),width=4,height=4) 
  # plot also spatially
  betasf <- cbind(tmpsf,betad %>% select(dist,is.above))
  t <- tm_shape(betasf) + tm_dots(fill="is.above",size=0.6,fill.scale = tm_scale_categorical(values=c("TRUE" = "#C11432", "FALSE" = "black")))  + tm_title(cond_site_name) + tm_layout(legend.position=c("right","top"),legend.height = 12) 
    tmap_save(t,filename=paste0("../Documents/iterative_deltas_res_margin/abovebelow_beta_map_",cond_site_name,".png"),width=3,height=6) 

    t <- tm_shape(betasf) + tm_dots(fill="b",size=0.6)  + tm_title(cond_site_name) + tm_layout(legend.position=c("right","top"),legend.height = 12) 
    tmap_save(t,filename=paste0("../Documents/iterative_deltas_res_margin/beta_map_",cond_site_name,".png"),width=3,height=6) 
    
    }

# look at all sites
x1 <- rep(0,ncol(df_sites))
x2 <- c(rep(60000,8),800000,900000,600000,600000)
y1 <- c(0.2,0.2,0.2,0.2,0.3,0.3,0.2,0.25,0.3,0.3,0.2,0.25)
y2 <- c(0.5,0.5,0.5,0.5,0.4,0.4,0.5,0.3,0.4,0.4,0.5,0.5)
sapply(1:ncol(df_sites),FUN = function(i) {beta_above_below(cond_site_name = names(df_sites)[i], x1 = x1[i], x2 = x2[i], y1 = y1[i], y2 = y2[i])})

# save the estimates to pass into iterative sigma_u
save(result, file="data_processed/iterative_delta_estimates.RData")

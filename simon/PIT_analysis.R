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
load("data_processed/N9000_sequential2_AGG_diagonal_sites_Birmingham_Cromer90.RData",verbose=TRUE)

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

# estimate all parameters for all sites
# find indeces of start and end sites
source("spatial_parameter_estimation.R") # for spatial_par_est function
Birmingham <- c(-1.9032,52.4806)
Cromer <- c(1.28486,53.05349)
site_start <- find_site_index(site=Birmingham,grid_uk = xyUK20_sf)
site_end <- find_site_index(site=Cromer,grid_uk = xyUK20_sf)

sites_index_diagonal <- c(192,193,194,195,196,197,218,219,220,221,242,263) # first is Birmingham and last is Cromer
site_name_diagonal <- c("Birmingham", paste0("diagonal",1:(length(sites_index_diagonal)-2)),"Cromer")
# plot these points on a map
uk_diag <- xyUK20_sf %>% mutate(siteID=as.numeric(1:nrow(xyUK20_sf)))
uk_diag <- uk_diag %>% mutate(sites_diagonal=factor(case_match(siteID,c(site_start) ~ "Birmingham",c(site_end)~"Cromer",sites_index_diagonal[2:(length(site_name_diagonal)-1)]~"diagonal_sites")))
tmap_mode("plot")

# estimate parameters along
spatial_par_est(data_Lap = data_mod_Lap_star,cond_sites = sites_index_diagonal,cond_site_names = site_name_diagonal,dayshift = 0,v=q,Ndays_season = 90,title = paste0("diagonal_sites_Birmingham_Cromer_star_",q*100))
# load estimates to calculate observed residuals
load("data_processed/N9000_sequential2_AGG_diagonal_sites_Birmingham_Cromer_star_90.RData",verbose=TRUE)
est_all_sf_star <- est_all_sf
for (i in 1:12) {
  print(site_name_diagonal[i])
  print(summary(est_all_sf_star %>% filter(cond_site==site_name_diagonal[i]) %>% select(a,b,mu_agg,sigl,sigu,deltal,deltau)-est_all_sf %>% filter(cond_site==site_name_diagonal[i]) %>% select(a,b,mu_agg,sigl,sigu,deltal,deltau)))
}

est_all_diag1 <- est_all_sf_star %>% mutate(cond_site = factor(as.character(cond_site),levels=as.character(site_name_diagonal)))
tm <- map_param(tmp_est=est_all_diag1, method = "AGG", facet_var = c("cond_site"), title_map = "Sites from Birmingham to Cromer")
condmodel_params <- c("a","b","mu","sig","muagg","sigl","sigu","sigdiff","deltal","deltau","deltadiff")

save_map_i <- function(i,tm_list=tm,doc_folder = "Birmingham_Cromer_diagonal_star",w = 8,h = 8,par_vect = condmodel_params) {
  tmap_save(tm_list[[i]],filename = paste0("../Documents/",doc_folder,"/allsites_",par_vect[i],"_original.png"),width = w,height = h) 
}
sapply(1:length(condmodel_params),FUN = save_map_i)

# estimate 6 phis
AGG_par_est_ite <- function(data_mod_Lap,site,v=0.9,Nite=10,sites = sites_index_diagonal,cond_site_names = site_name_diagonal,q=0.9,grid20km=xyUK20_sf,est_all_sf,deltal=NULL,deltau=NULL,phi0l=NULL,folder_name=NULL) {
  # filter conditioning site name from a vector of names
  if(is.null(cond_site_names)) {
    cond_site_name <- names(sites)[site]
    cond_site_names <- names(sites)
  } else {  cond_site_name <- cond_site_names[site] }
  # filter conditioning site index from a vector of indeces
  if (is.numeric(sites)) {cond_site <- sites[site]} else{
    cond_site_coord <- sites %>% dplyr::select(all_of(cond_site_name)) %>% pull()
    cond_site <- find_site_index(cond_site_coord,grid_uk = grid20km)    }
  print(cond_site)
  # calculate distance from the conditioning site
  dist_tmp <- as.numeric(unlist(st_distance(grid20km[cond_site,],grid20km)))
  # remove zero distance (conditioning site)
  dist_tmp <- dist_tmp[dist_tmp>1]
  # normalise distance using a common constant
  distnorm <- dist_tmp/1000000
  print(summary(distnorm))
  parest_site <- est_all_sf %>% dplyr::filter(cond_site==cond_site_name) %>% dplyr::select(sigl,sigu,deltal,deltau) %>% na.omit()
  names(parest_site) <- c("sigl","sigu","deltal","deltau")
  # calculate observed residuals
  aest <- discard(est_all_sf %>% dplyr::filter(cond_site==cond_site_name) %>% pull(a),is.na)
  print("Check alpha")
  print(summary(aest))
  print(length(aest))
  best <- discard(est_all_sf %>% dplyr::filter(cond_site==cond_site_name) %>% pull(b),is.na)
  print("Check beta")
  print(summary(best))
  print(length(best))
  Z <- observed_residuals(df=data_mod_Lap,given=cond_site,v = q,a=aest,b=best)
  print(names(Z)[180:220])
  try7 <- par_est_ite(z=Z,v=v,given=cond_site,cond_site_dist=distnorm, parest_site = parest_site,Nite=Nite,show_ite=TRUE,deltal=deltal,deltau= deltau,phi0l=phi0l) 
  
  # separate parameter estimation and analysis
  p <- ggplot(try7$mu_agg %>% pivot_longer(everything(),names_to = "iteration",values_to = "par") %>% mutate(iteration=factor(iteration,levels=paste0("X",1:Nite)))) + geom_boxplot(aes(x=iteration,y=par)) +ylab(TeX("$\\mu_{AGG}$"))
  ggsave(p,file=paste0("../Documents/",folder_name,"/boxplot_mu_agg_",cond_site_name,".png"),height=5,width=10)
  p <- ggplot(try7$sigl %>% pivot_longer(everything(),names_to = "iteration",values_to = "par")%>% mutate(iteration=factor(iteration,levels=paste0("X",1:Nite)))) + geom_boxplot(aes(x=iteration,y=par)) +ylab(TeX("$\\sigma_l$"))
  ggsave(p,file=paste0("../Documents/",folder_name,"/boxplot_sigmal_",cond_site_name,".png"),height=5,width=10)
  p <- ggplot(try7$sigu %>% pivot_longer(everything(),names_to = "iteration",values_to = "par")%>% mutate(iteration=factor(iteration,levels=paste0("X",1:Nite)))) + geom_boxplot(aes(x=iteration,y=par)) +ylab(TeX("$\\sigma_u$"))
  ggsave(p,file=paste0("../Documents/",folder_name,"/boxplot_sigmau_",cond_site_name,".png"),height=5,width=10)
  
  tmp_phis <- data.frame("param_value"=c(try7$phi0l,try7$phi1l,try7$phi2l,try7$phi0u,try7$phi1u,try7$phi2u),"iteration"=paste0("X",0:Nite), "phi" = rep(c("phi0l","phi1l","phi2l","phi0u","phi1u","phi2u"),each=Nite+1))
  p <- ggplot(tmp_phis %>% mutate(iteration=factor(iteration,levels=paste0("X",0:Nite)))) + geom_point(aes(x=iteration,y=param_value,col=phi)) +ylab(TeX("$\\phi$"))
  ggsave(p,file=paste0("../Documents/",folder_name,"/iteration_phis_",cond_site_name,".png"),height=5,width=10)
  
  
  # look spatially to check
  # explore also spatial parameters
  est_ite <- try7$par_sum %>% add_row(.before=cond_site)
  names(est_ite) <- paste0(names(est_ite),"_ite")
  tmpsf <- cbind(est_all_sf %>% filter(cond_site==cond_site_name),est_ite)
  # plot parameter estimates against distance
  # calculate distance from a conditioning site with st_distance()
  
  mud <- data.frame(mu=tmpsf$mu_agg_ite,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=mu,x=dist))
  ggsave(mud,filename=paste0("../Documents/",folder_name,"/muagg_distance_",cond_site_name,".png")) 
  
  sigld <- data.frame(sigl=tmpsf$sigl_ite,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=sigl,x=dist))
  ggsave(sigld,filename=paste0("../Documents/",folder_name,"/sigl_distance_",cond_site_name,".png")) 
  
  sigud <- data.frame(sigu=tmpsf$sigu_ite,dist=as.numeric(unlist(st_distance(tmpsf[cond_site,],tmpsf)))) %>% ggplot() + geom_point(aes(y=sigu,x=dist))
  ggsave(sigud,filename=paste0("../Documents/",folder_name,"/sigu_distance_",cond_site_name,".png")) 
  
  tmpsf <- tmpsf %>% mutate(mudiff = mu_agg_ite - mu_agg, sigldiff = sigl_ite - sigl, sigudiff = sigu_ite - sigu)
  toplabel <- c("New iterative approach","Original method","Difference")
  
  misscol <- "aquamarine"
    #mu_limits <- c(-2.61,1)
  t <- tmpsf %>% dplyr::select(mu_agg_ite,mu_agg,mudiff) %>% pivot_longer(cols=c(mu_agg_ite,mu_agg,mudiff),names_to = "parameter", values_to = "value" ) %>% mutate(parameter=factor(parameter,levels=c("mu_agg_ite","mu_agg","mudiff"))) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="-brewer.rd_bu",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\mu_{AGG}$"))) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 12, panel.labels = toplabel,legend.reverse = TRUE) 
  
  tmap_save(t,filename=paste0("../Documents/",folder_name,"/mu_agg_",cond_site_name,".png"),width=8,height=6)
  
  #sigma_limits <- c(-1.32,2.3)
  t <- tmpsf %>% dplyr::select(sigl_ite,sigl,sigldiff) %>% pivot_longer(cols=c(sigl_ite,sigl,sigldiff),names_to = "parameter", values_to = "value" ) %>% mutate(parameter=factor(parameter,levels=c("sigl_ite","sigl","sigldiff"))) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="-brewer.rd_bu",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\sigma_l$"))) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 9, panel.labels = toplabel,legend.reverse = TRUE) 
  tmap_save(t,filename=paste0("../Documents/",folder_name,"/sigma_lower_",cond_site_name,".png"),width=8,height=6)
  
  t <- tmpsf %>% dplyr::select(sigu_ite,sigu,sigudiff) %>% pivot_longer(cols=c(sigu_ite,sigu,sigudiff),names_to = "parameter", values_to = "value" ) %>% mutate(parameter=factor(parameter,levels=c("sigu_ite","sigu","sigudiff"))) %>% tm_shape() + tm_dots(fill="value",size=0.5,fill.scale =tm_scale_continuous(values="-brewer.rd_bu",value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$\\sigma_u$"))) + tm_facets("parameter") +
    tm_layout(legend.position=c("right","top"),legend.height = 9, panel.labels = toplabel,legend.reverse = TRUE) 
  tmap_save(t,filename=paste0("../Documents/",folder_name,"/sigma_upper_",cond_site_name,".png"),width=8,height=6)
  return(try7)
}
# load result_new for deltal and deltau
load("data_processed/iterative_phi0l_phi0u_estimates_Birmingham_Cromer_diagonal.RData",verbose=TRUE)
# requires loading functions from iterative_sigmas_Jan_2026
deltal <- result_new[[1]]$deltal[1]
deltau <- result_new[[1]]$deltau[1]
folder_name <- "Birmingham_Cromer_diagonal/new_iterative_sigmas_mu_star"
result_new_star <- sapply(1:12,FUN = function(site_order){AGG_par_est_ite(site=site_order,data_mod_Lap = data_mod_Lap_star,sites = sites_index_diagonal,cond_site_names = site_name_diagonal,est_all_sf = est_all_sf_star,Nite=10,deltal=deltal,deltau=deltau,folder_name = folder_name,phi0l=NULL)},simplify = FALSE)

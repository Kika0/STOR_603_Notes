library(tmap) # spatial map plots
library(sf) # for handling spatial sf objects
library(viridis)
library(tidyverse)
library(mvtnorm) # for multivariate normal distribution functions
library(units) # manage units for distance
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

sites_index_diagonal <- c(192,193,194,195,196,197,218,219,220,221,242,263) # first is Birmingham and last is Cromer
site_name_diagonal <- c("Birmingham", paste0("diagonal",1:(length(sites_index_diagonal)-2)),"Cromer")

# residuals: new iterative approach with 4 phi parameters ---------------------
load(file="data_processed/iterative_sigmas_estimates_Birmingham_Cromer_diagonal.RData",verbose = TRUE)
pe <- as.data.frame(result[[1]] %>% select(a,b,mu_agg_ite_sig,sigl_ite_sig,sigu_ite_sig,deltal_ite,deltau_ite)) %>% drop_na()
names(pe) <- c("a","b","mu","sigl","sigu","deltal","deltau")
# calculate observed residuals
Z <- observed_residuals(df=data_mod_Lap,given=sites_index_diagonal[1],v = q,a=pe$a,b=pe$b)

# transform to standard Normal
AGG_Normal_PIT <- function(z,theta) {
  return(qnorm(pAGG(x=z,theta=theta)))
}

ZU <- sapply(1:ncol(Z),FUN=function(k) {pAGG(x = Z[,k],theta=c(pe$mu[k],pe$sigl[k],pe$sigu[k],pe$deltal[1],pe$deltau[1]))})
ZN <- sapply(1:ncol(Z),FUN=function(k) {AGG_Normal_PIT(z = Z[,k],theta=c(pe$mu[k],pe$sigl[k],pe$sigu[k],pe$deltal[1],pe$deltau[1]))})
# add empty row for Birmingham and transform back to spatial
est_iteN <- as.data.frame(t(ZN)) %>% add_row(.before=sites_index_diagonal[1])
est_iteU <- as.data.frame(t(ZU)) %>% add_row(.before=sites_index_diagonal[1])

# find 10 largest days
hot_temps <- sort(as.numeric(unlist(data_mod_Lap[,sites_index_diagonal[1]][data_mod_Lap[,sites_index_diagonal[1]]>quantile(data_mod_Lap[,sites_index_diagonal[1]],q)])),decreasing=TRUE)[1:10]
# find the corresponding rows
hot_temps_index <- sort(as.numeric(unlist(data_mod_Lap[,sites_index_diagonal[1]][data_mod_Lap[,sites_index_diagonal[1]]>quantile(data_mod_Lap[,sites_index_diagonal[1]],q)])),decreasing=TRUE,index.return=TRUE)$ix[1:10]
# check data
for (i in hot_temps_index) {
print(summary(as.numeric(unlist(Z[i,]))) ) }
Zemp <- as.data.frame((Z %>% apply(c(2),FUN=row_number))/(nrow(Z)+1)) 
for (i in hot_temps_index) {
  print(summary(as.numeric(unlist(Zemp[i,]))) ) 
  print(summary(as.numeric(unlist(ZU[i,]))) ) }

tmpN <- st_as_sf(cbind(est_iteN %>% dplyr::select(all_of(hot_temps_index)),result[[1]] %>% dplyr::select(geometry)) %>% pivot_longer(cols=starts_with("V")))
tmpU <- st_as_sf(cbind(est_iteU %>% dplyr::select(all_of(hot_temps_index)),result[[1]] %>% dplyr::select(geometry)) %>% pivot_longer(cols=starts_with("V")))
misscol <- "aquamarine"
tN <- tm_shape(tmpN %>% mutate("name"=factor(name, levels=unique(tmpN$name)))) + tm_dots(fill="value",size=0.5,fill.scale = tm_scale_continuous(values="-brewer.rd_bu",limits=c(-8,8),value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$Z^N$"))) + tm_facets(by="name",ncol=5)
tU <- tm_shape(tmpU %>% mutate("name"=factor(name, levels=unique(tmpU$name)))) + tm_dots(fill="value",size=0.5,fill.scale = tm_scale_continuous(values="-brewer.rd_bu",limits=c(0,1),value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$Z^U$"))) + tm_facets(by="name",ncol=5)
tmap_save(tm=tN, filename=paste0("../Documents/residual_dependence_Birmingham_normal.png"),width=10,height=8)
tmap_save(tm=tU, filename=paste0("../Documents/residual_dependence_Birmingham_uniform.png"),width=10,height=8)

# repeat with residuals: new iterative approach with 6 phi parameters -------------------------------
load(file="data_processed/iterative_phi0l_phi0u_estimates_Birmingham_Cromer_diagonal.RData",verbose = TRUE)
aest <- discard(est_all_sf %>% filter(cond_site==site_name_diagonal[1]) %>% pull(a),is.na)
best <- discard(est_all_sf %>% filter(cond_site==site_name_diagonal[1]) %>% pull(b),is.na)
pe <- as.data.frame(result_new[[1]][[12]] %>% dplyr::select(mu_agg,sigl,sigu,deltal,deltau))
names(pe) <- c("mu","sigl","sigu","deltal","deltau")
# calculate observed residuals
Z <- observed_residuals(df=data_mod_Lap,given=sites_index_diagonal[1],v = q,a=aest,b=best)
ZU <- sapply(1:ncol(Z),FUN=function(k) {pAGG(x = Z[,k],theta=c(pe$mu[k],pe$sigl[k],pe$sigu[k],pe$deltal[1],pe$deltau[1]))})
ZN <- sapply(1:ncol(Z),FUN=function(k) {AGG_Normal_PIT(z = Z[,k],theta=c(pe$mu[k],pe$sigl[k],pe$sigu[k],pe$deltal[1],pe$deltau[1]))})
# add empty row for Birmingham and transform back to spatial
est_iteU <- as.data.frame(t(ZU)) %>% add_row(.before=sites_index_diagonal[1])
est_iteN <- as.data.frame(t(ZN)) %>% add_row(.before=sites_index_diagonal[1])
# check data
for (i in hot_temps_index) {
  print(summary(as.numeric(unlist(Z[i,]))) ) }
Zemp <- as.data.frame((Z %>% apply(c(2),FUN=row_number))/(nrow(Z)+1)) 
for (i in hot_temps_index) {
  print(summary(as.numeric(unlist(Zemp[i,]))) ) 
  print(summary(as.numeric(unlist(ZU[i,]))) ) }

tmpN <- st_as_sf(cbind(est_iteN %>% dplyr::select(all_of(hot_temps_index)),result[[1]] %>% dplyr::select(geometry)) %>% pivot_longer(cols=starts_with("V")))
tmpU <- st_as_sf(cbind(est_iteU %>% dplyr::select(all_of(hot_temps_index)),result[[1]] %>% dplyr::select(geometry)) %>% pivot_longer(cols=starts_with("V")))
misscol <- "aquamarine"
tN <- tm_shape(tmpN %>% mutate("name"=factor(name, levels=unique(tmpN$name)))) + tm_dots(fill="value",size=0.5,fill.scale = tm_scale_continuous(values="-brewer.rd_bu",limits=c(-8,8),value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$Z^N$"))) + tm_facets(by="name",ncol=5)
tU <- tm_shape(tmpU %>% mutate("name"=factor(name, levels=unique(tmpU$name)))) + tm_dots(fill="value",size=0.5,fill.scale = tm_scale_continuous(values="-brewer.rd_bu",limits=c(0,1),value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$Z^U$"))) + tm_facets(by="name",ncol=5)
tmap_save(tm=tN, filename=paste0("../Documents/new_residual_dependence_Birmingham_normal.png"),width=10,height=8)
tmap_save(tm=tU, filename=paste0("../Documents/new_residual_dependence_Birmingham_uniform.png"),width=10,height=8)

# repeat with residuals: separate parameters -------------------------------
aest <- discard(est_all_sf %>% filter(cond_site==site_name_diagonal[1]) %>% pull(a),is.na)
best <- discard(est_all_sf %>% filter(cond_site==site_name_diagonal[1]) %>% pull(b),is.na)
pe <- as.data.frame(result[[1]] %>% select(mu_agg,sigl,sigu,deltal,deltau)) %>% drop_na()
names(pe) <- c("mu","sigl","sigu","deltal","deltau")
# calculate observed residuals
Z <- observed_residuals(df=data_mod_Lap,given=sites_index_diagonal[1],v = q,a=aest,b=best)
ZU <- sapply(1:ncol(Z),FUN=function(k) {pAGG(x = Z[,k],theta=c(pe$mu[k],pe$sigl[k],pe$sigu[k],pe$deltal[1],pe$deltau[1]))})
ZN <- sapply(1:ncol(Z),FUN=function(k) {AGG_Normal_PIT(z = Z[,k],theta=c(pe$mu[k],pe$sigl[k],pe$sigu[k],pe$deltal[1],pe$deltau[1]))})
# add empty row for Birmingham and transform back to spatial
est_iteU <- as.data.frame(t(ZU)) %>% add_row(.before=sites_index_diagonal[1])
est_iteN <- as.data.frame(t(ZN)) %>% add_row(.before=sites_index_diagonal[1])
# check data
for (i in hot_temps_index) {
  print(summary(as.numeric(unlist(Z[i,]))) ) }
Zemp <- as.data.frame((Z %>% apply(c(2),FUN=row_number))/(nrow(Z)+1)) 
for (i in hot_temps_index) {
  print(summary(as.numeric(unlist(Zemp[i,]))) ) 
  print(summary(as.numeric(unlist(ZU[i,]))) ) }

tmpN <- st_as_sf(cbind(est_iteN %>% dplyr::select(all_of(hot_temps_index)),result[[1]] %>% dplyr::select(geometry)) %>% pivot_longer(cols=starts_with("V")))
tmpU <- st_as_sf(cbind(est_iteU %>% dplyr::select(all_of(hot_temps_index)),result[[1]] %>% dplyr::select(geometry)) %>% pivot_longer(cols=starts_with("V")))
misscol <- "aquamarine"
  tN <- tm_shape(tmpN %>% mutate("name"=factor(name, levels=unique(tmpN$name)))) + tm_dots(fill="value",size=0.5,fill.scale = tm_scale_continuous(values="-brewer.rd_bu",limits=c(-8,8),value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$Z^N$"))) + tm_facets(by="name",ncol=5)
  tU <- tm_shape(tmpU %>% mutate("name"=factor(name, levels=unique(tmpU$name)))) + tm_dots(fill="value",size=0.5,fill.scale = tm_scale_continuous(values="-brewer.rd_bu",limits=c(0,1),value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$Z^U$"))) + tm_facets(by="name",ncol=5)
  tmap_save(tm=tN, filename=paste0("../Documents/separate_residual_dependence_Birmingham_normal.png"),width=10,height=8)
  tmap_save(tm=tU, filename=paste0("../Documents/separate_residual_dependence_Birmingham_uniform.png"),width=10,height=8)

# also include empirical residuals ---------------------------------------------
ZN <- sapply(1:ncol(Z),FUN=function(k) {qnorm(Zemp[,k])})
est_iteU <- as.data.frame(t(Zemp)) %>% add_row(.before=sites_index_diagonal[1])
est_iteN <- as.data.frame(t(ZN)) %>% add_row(.before=sites_index_diagonal[1])
tmpN <- st_as_sf(cbind(est_iteN %>% dplyr::select(all_of(hot_temps_index)),result[[1]] %>% dplyr::select(geometry)) %>% pivot_longer(cols=starts_with("V")))
tmpU <- st_as_sf(cbind(est_iteU %>% dplyr::select(all_of(hot_temps_index)),result[[1]] %>% dplyr::select(geometry)) %>% pivot_longer(cols=starts_with("V")))
misscol <- "aquamarine"
tN <- tm_shape(tmpN %>% mutate("name"=factor(name, levels=unique(tmpN$name)))) + tm_dots(fill="value",size=0.5,fill.scale = tm_scale_continuous(values="-brewer.rd_bu",limits=c(-8,8),value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$Z^N$"))) + tm_facets(by="name",ncol=5)
tU <- tm_shape(tmpU %>% mutate("name"=factor(name, levels=unique(tmpU$name)))) + tm_dots(fill="value",size=0.5,fill.scale = tm_scale_continuous(values="-brewer.rd_bu",limits=c(0,1),value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = TeX("$Z^U$"))) + tm_facets(by="name",ncol=5)
tmap_save(tm=tN, filename=paste0("../Documents/empirical_residual_dependence_Birmingham_normal.png"),width=10,height=8)
tmap_save(tm=tU, filename=paste0("../Documents/empirical_residual_dependence_Birmingham_uniform.png"),width=10,height=8)
  
  
# calculate correlation matrix -------------------------------------------------
gaus_cor <- function(i,j,h,cond_index=192) {
  (mat_cor(h[i,j]) - mat_cor(h[cond_index,i])* mat_cor(h[cond_index,j]) )/ ((1-(mat_cor(h[cond_index,i])^2))^(1/2) * (1-mat_cor(h[cond_index,j])^2)^(1/2))
}
library(fields)
mat_cor <- function(x) {
  fields::Matern(x,range=200,smoothness=1)
}

gaus_cov <- function(i,j,h,cond_index=192) {
  (mat_cor(h[i,j]) - mat_cor(h[cond_index,i])* mat_cor(h[cond_index,j]) )
}
# calculate distance matrix
h <- (st_distance(xyUK20_sf,xyUK20_sf) %>% drop_units() )/1000

Zcov <- matrix(ncol=ncol(Z)+1,nrow=ncol(Z)+1)
for (i in 1:(ncol(Z)+1)) {
  for (j in 1:(ncol(Z)+1)) {
    Zcov[i,j] <- gaus_cov(i=i,j=j,h=h)
  }
}
summary(Zcov[,190:199])
# replace NA with 0
cond_index <- 192
Zcov[,cond_index] <- Zcov[cond_index,] <- 0
Zcov[cond_index,cond_index] <- 1
# random10 <- as.data.frame(t(  rmvnorm(n=100000,Sigma = Zcov)  ))
# top10_index <- sort(abs(as.numeric(unlist(random10[192,]))), index.return=TRUE,decreasing = FALSE)$ix[1:10]
# random10 <- random10[,top10_index]
random10 <- as.data.frame(t(  rmvnorm(n=10,Sigma = Zcov)  ))
names(random10) <- paste0("random",1:10)
random10[cond_index,] <- NA
tmpsf <- st_as_sf(cbind(random10,xyUK20_sf)) %>% pivot_longer(cols=contains("random"))

t <- tm_shape(tmpsf %>% mutate("name"=factor(name, levels=unique(tmpsf$name)))) + tm_dots(fill="value",size=0.5,fill.scale = tm_scale_continuous(values="-brewer.rd_bu",limits=c(-4,4),value.na=misscol,label.na = "Conditioning\n site"),fill.legend = tm_legend(title = "")) + tm_facets(by="name",ncol=5)
tmap_save(tm=t, filename=paste0("../Documents/random10_residual_dependence_Birmingham_normal.png"),width=10,height=8)

NLL_range <- function(x) {
  dmvnorm()
}

library(tidyverse)
library(tmap)
library(tmaptools)
library(sf)
library(lubridate)
library(gridExtra)
library(LaplacesDemon)
library(latex2exp)
library(evd)
theme_set(theme_bw())
theme_replace(
  panel.spacing = unit(2, "lines"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA) )
folder_name <- "../Documents/spatial_model_final_steps/"

cond_site_name <- "London"
# load observed data
#source("spatial_parameter_estimation.R") # for spatial_par_est function
load("data_processed/temperature_data.RData",verbose = TRUE)
load("data_processed/spatial_helper.RData", verbose = TRUE)
source("simon/P2q_function_helpers.R")
load("data_processed/P2qselected_helpers.RData", verbose = TRUE)
load(paste0("data_processed/N9000_sequential2_AGG_all12sites",q*100,".RData"),verbose = TRUE) # original parameter estimates
load("data_processed/iterative_phi0l_phi0u_estimates_London.RData",verbose=TRUE) # residual margin parameters
load("data_processed/residual_dependence_pars.RData", verbose = TRUE) # residual dependence parameters
deltal <- result_new[[12]]$deltal[1]
deltau <- result_new[[12]]$deltau[1]
source("simon/final_model_helpers.R")
# Model 3: parameter estimation ------------------------------------------------
par_est_model_3 <- function(cond_index,v=0.9,data_Lap=data_mod_Lap,grid20km=xyUK20_sf,deltal,deltau,Nite_2_3=10,Nite_phi=10) {
  # calculate distance from the conditioning site
  dist_tmp <- as.numeric(unlist(st_distance(grid20km[cond_index,],grid20km[-cond_index,])))
  distnorm <- dist_tmp/1000000  # normalise distance using a common constant
  # subset conditioning dataset
  data_Lapv <- data_Lap %>% filter(data_Lap[,cond_index]>quantile(data_Lap[,cond_index],v))
  # initiate objects
  x2 <- x3 <- list() # a list of output at steps 2 and 3
  x2_df <- data.frame("mu_agg"=numeric(),"sigl"=numeric(),"sigu"=numeric(),"phi0u"=numeric(),"phi1u"=numeric(),"phi2u"=numeric(),"phi0l"=numeric(),"phi1l"=numeric(),"phi2l"=numeric(),"iteration"=numeric())
  x3_df <- data.frame("a"=numeric(),"b"=numeric(),"mu"=numeric(),"iteration"=numeric())
  x_time <- data.frame("step2"=numeric(),"step3"=numeric(),"iteration"=numeric())
  # 1. estimate alpha and beta
  x1 <- par_est(df=data_Lapv,v = v,given = cond_index,margin = "Normal",method = "sequential2",keef_constraints = c(1,2))
  # iterate between steps 2. and 3.
  for (Nite_i in 1:Nite_2_3) {
  if (Nite_i>1) {
    a <- tmp$a[!is.na(tmp$a)]
    b <- tmp$b[!is.na(tmp$b)]
    pe_res <- pe_phi_ite %>% dplyr::select(sigl,sigu) %>% na.omit() 
  } else {
    pe_res <- data.frame("sigl" = rep(1,length(a)),"sigu" = rep(1,length(a)))
    a <- x1$a[!is.na(x1$a)]
    b <- x1$b[!is.na(x1$b)]
  }
  # update residuals
  Z <- observed_residuals(df=data_Lap,given=cond_index,v = v,a=a,b=b)
  # 2. estimate phis
  x_time_dummy <- Sys.time()
  x2 <- par_est_ite(z=Z,given=cond_index,cond_site_dist = distnorm, parest_site = pe_res,Nite = Nite_phi,show_ite=TRUE,deltal = deltal,deltau = deltau) 
  x_time2 <- Sys.time()-x_time_dummy
  x2_df <- rbind(x2_df,x2[[12]] %>% mutate("iteration"=Nite_i))
  # 3. reestimate a,b,mu
  pe_phi_ite <- x2[[12]] %>% dplyr::select(mu_agg,sigl,sigu,deltal,deltau)
  x_time_dummy <- Sys.time()
  x3 <- sapply(1:ncol(data_Lap),FUN=NLL_AGG_wrapper,data_Lapv=data_Lapv,cond_index=cond_index,pe_res = pe_phi_ite%>% add_row(.before = cond_index))
  x_time3 <- Sys.time()-x_time_dummy
  tmp <- as.data.frame(do.call(rbind,x3)) %>% mutate("iteration"=Nite_i)
  names(tmp)[1:3] <- c("a","b","mu")
  x3_df <- rbind(x3_df,tmp)
  x_time <- rbind(x_time,data.frame("step2"=x_time2,"step3"=x_time3,"iteration"=Nite_i))
  }
  return(list(x1,x2_df,x3_df,x_time))
}
q <- 0.9 # set quantile threshold
s <- Sys.time()
y <- par_est_model_3(cond_index=London_index,v=q,data_Lap = data_mod_Lap,deltal = deltal,deltau = deltau)
Sys.time()-s

# plot diagnostics
folder_name <- "../Documents/final_model_3/"
ggplot(y[[3]]) + geom_point(aes(x=a,y=b)) + facet_wrap(~iteration)
ggplot(y[[3]] %>% mutate(iteration=factor(iteration))) + geom_boxplot(aes(x=iteration,y=a)) 
ggplot(y[[3]] %>% mutate(iteration=factor(iteration))) + geom_boxplot(aes(x=iteration,y=b)) 

# 1. map alpha and beta new and original estimates ----------------------------
a_new <- y[[3]] %>% dplyr::filter(iteration==Nite_2_3) %>% dplyr::select(a) %>% pull(a)
b_new <- y[[3]] %>% dplyr::filter(iteration==Nite_2_3) %>% dplyr::select(b)  %>% pull(b)
title_map <- ""
misscol <- "aquamarine"
legend_text_size <- 0.7
point_size <- 0.5
legend_title_size <- 0.9
limsa <- c(0,1)
limsb <- c(0,0.65)
nrow_facet <- 1
estsf <- cbind(est_all_sf %>% filter(cond_site %in% cond_site_name), data.frame("a_new"=a_new,"b_new" = b_new))
p1 <- tm_shape(estsf) + tm_dots(fill="a",fill.scale = tm_scale_continuous(limits=limsa,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\alpha$"))) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text=TeX("$\\tilde{\\alpha}$")) 
p2 <- tm_shape(estsf) + tm_dots(fill="b",fill.scale = tm_scale_continuous(limits=limsb,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\beta$"))) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text=TeX("$\\tilde{\\beta}$")) 
p3 <- tm_shape(estsf) + tm_dots(fill="a_new",fill.scale = tm_scale_continuous(limits=limsa,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\alpha$"))) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text=TeX("$\\hat{\\alpha}$")) 
p4 <- tm_shape(estsf) + tm_dots(fill="b_new",fill.scale = tm_scale_continuous(limits=limsb,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\beta$"))) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text=TeX("$\\hat{\\beta}$")) 
tmap_save(tmap_arrange(p1,p2,p3,p4,ncol=4),filename=paste0(folder_name,"alpha_beta_fixed_res",cond_site_name,".png"),height=6,width=11)

tmp1 <- rbind(data.frame(a=aest,b=best,"method"="aoriginal"),data.frame(a=na.omit(tmp$a),b=na.omit(tmp$b),"method"="new")) %>% mutate("iteration"=rep(1:length(aest),2))
# map alpha and beta original and new estimates
plot_ab <- function(tmp) { ggplot(tmp) + 
    geom_line(aes(x=a,y=b,group=iteration),linewidth=0.1) +
    geom_point(aes(x=a,y=b,col=method),alpha=0.7,size=1) +
    xlab(TeX("${\\alpha}$")) +
    ylab(TeX("${\\beta}$")) + 
    scale_color_manual(values = c("aoriginal" = "#009ADA", "new" = "#C11432"),labels = c("Model 1","Model 2")) + coord_fixed() + theme(axis.text.y = element_text(angle = 90, vjust = 0.5)) + labs(col="")
}
p <- plot_ab(tmp=tmp1)
ggsave(p,filename=paste0(folder_name,"plot_ab_new_original.png"),width=7,height=4) 
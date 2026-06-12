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
# load observed data
#source("spatial_parameter_estimation.R") # for spatial_par_est function
load("data_processed/temperature_data.RData",verbose = TRUE)
load("data_processed/spatial_helper.RData", verbose = TRUE)
source("simon/P2q_function_helpers.R")
load("data_processed/P2qselected_helpers.RData", verbose = TRUE)
load(paste0("data_processed/N9000_sequential2_AGG_all12sites",q*100,".RData"),verbose = TRUE) # original parameter estimates
load("data_processed/iterative_phi0l_phi0u_estimates_London.RData",verbose=TRUE) # residual margin parameters for delta constants
load("data_processed/residual_dependence_pars.RData", verbose = TRUE) # residual dependence parameters
deltal <- result_new[[12]]$deltal[1]
deltau <- result_new[[12]]$deltau[1]
source("simon/final_model_helpers.R")
site_i <- 7
cond_index <- df_sites[3,site_i]
cond_site_name <- names(df_sites)[site_i]

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
    a <- x1$a[!is.na(x1$a)]
    b <- x1$b[!is.na(x1$b)]
    pe_res <- data.frame("sigl" = rep(1,length(a)),"sigu" = rep(1,length(a)))
  }
  # update residuals
  Z <- observed_residuals(df=data_Lap,given=cond_index,v = v,a=a,b=b)
  print(v) # check value of quantile threshold
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
y <- par_est_model_3(cond_index=cond_index,v=q,data_Lap = data_mod_Lap,deltal = deltal,deltau = deltau)
Sys.time()-s

# plot diagnostics
folder_name <- "../Documents/final_model_3/"
ggplot(y[[3]]) + geom_point(aes(x=a,y=b)) + facet_wrap(~iteration)
ggplot(y[[3]] %>% mutate(iteration=factor(iteration))) + geom_boxplot(aes(x=iteration,y=a)) 
ggplot(y[[3]] %>% mutate(iteration=factor(iteration))) + geom_boxplot(aes(x=iteration,y=b)) 

# 1. map alpha and beta new and original estimates ----------------------------
a_orig <- y[[1]]$a
b_orig <- y[[1]]$b
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

# 2. scatterplot comparing alpha and beta original and new estimates ----------
tmp1 <- rbind(data.frame(a=na.omit(a_orig),b=na.omit(b_orig),"method"="aoriginal"),data.frame(a=na.omit(tmp$a),b=na.omit(tmp$b),"method"="new")) %>% mutate("iteration"=rep(1:length(a_orig),2))
plot_ab <- function(tmp) { ggplot(tmp) + 
    geom_line(aes(x=a,y=b,group=iteration),linewidth=0.1) +
    geom_point(aes(x=a,y=b,col=method),alpha=0.7,size=1) +
    xlab(TeX("${\\alpha}$")) +
    ylab(TeX("${\\beta}$")) + 
    scale_color_manual(values = c("aoriginal" = "#009ADA", "new" = "#C11432"),labels = c("Model 1","Model 3")) + coord_fixed() + theme(axis.text.y = element_text(angle = 90, vjust = 0.5)) + labs(col="")
}
p <- plot_ab(tmp=tmp1)
ggsave(p,filename=paste0(folder_name,"plot_ab_new_original",cond_site_name,".png"),width=7,height=4) 

# 3. NLL comparing original and new estimates
pe <- y[[2]] %>% dplyr::filter(iteration==Nite_2_3) %>% dplyr::select(sigl,sigu,deltal,deltau)
abmu <- y[[3]] %>% dplyr::filter(iteration==Nite_2_3) %>% dplyr::select(a,b,mu) %>% na.omit()
# calculate log-likelihood
NLL_AGG_wrapper <- function(data_Lap=data_mod_Lap,i,pe_res,abmu,cond_index,v=0.9,a=na.omit(a_new),b=na.omit(b_new),mu=na.omit(mu_new)) {
  pe_i <- as.numeric(unlist(pe_res[i,]))
  abmu_i <- as.numeric(unlist(abmu[i,]))
  data_Lapv <- data_Lap %>% filter(data_Lap[,cond_index]>quantile(data_Lap[,cond_index],v))
  x1 <- as.numeric(unlist(data_Lapv[,cond_index]))
  ires <- c(1:length(data_Lap))[-cond_index][i]
  x2 <- as.numeric(unlist(data_Lapv[,ires]))
  y <- NLL_AGG_onestep(x=data.frame(x1,x2),theta=c(),a_hat=abmu_i[1],b_hat=abmu_i[2],mu_hat=abmu_i[3],sigl_hat = pe_i[1], sigu_hat = pe_i[2], deltal_hat = pe_i[3], deltau_hat = pe_i[4])
  return(y) 
}
nll3 <- sapply(1:(ncol(data_mod_Lap)-1),FUN=NLL_AGG_wrapper,data_Lap=data_mod_Lap,cond_index=cond_index,pe_res=pe,abmu=abmu)

# calculate also Model 1 for comparison
pe_orig <- y[[1]] %>% dplyr::select(a,b,mu,sig) %>% na.omit()
NLL_AGG_wrapper <- function(data_Lap=data_mod_Lap,i,cond_index,v=0.9,pe_orig) {
  data_Lapv <- data_Lap %>% filter(data_Lap[,cond_index]>quantile(data_Lap[,cond_index],v))
  pe_i <- as.numeric(unlist(pe_orig[i,]))
  x1 <- as.numeric(unlist(data_Lapv[,cond_index]))
  ires <- c(1:length(data_Lap))[-cond_index][i]
  x2 <- as.numeric(unlist(data_Lapv[,ires]))
  y <- NLL_AGG_onestep(x=data.frame(x1,x2),theta=c(),a_hat=pe_i[1],b_hat=pe_i[2],mu_hat=pe_i[3],sigl_hat = sqrt(2)* pe_i[4], sigu_hat = sqrt(2)* pe_i[4], deltal_hat = 2, deltau_hat = 2)
  return(y) 
}
nll1 <- sapply(1:(ncol(data_mod_Lap)-1),FUN=NLL_AGG_wrapper,data_Lap=data_mod_Lap,cond_index=cond_index,pe_orig=pe_orig)

# calculate AIC
a1 <- 2*nll1+2*4
a3 <- 2*nll3+2*3
tmp <- data.frame(nll1,nll3,a1,a3)

# 3. boxplots of NLL and AIC --------------------------------------------------
tmp1 <- tmp %>% dplyr::select(contains("nll")) %>% pivot_longer(cols=everything())
p <- ggplot(tmp1) + geom_boxplot(aes(y=value,x=name)) + labs(y="Negative log-likelihood",x="")
ggsave(p,filename=paste0(folder_name,"nll_compare_",cond_site_name,".png"),width=9,height=5)
# boxplots of AIC
tmp1 <- tmp %>% dplyr::select(contains("a")) %>% pivot_longer(cols=everything())
p <- ggplot(tmp1) + geom_boxplot(aes(y=value,x=name)) + labs(y="AIC",x="")
ggsave(p,filename=paste0(folder_name,"AIC_compare_",cond_site_name,".png"),width=9,height=5)

# 4. map NLL and AIC ------------------------------------------------
title_map <- ""
misscol <- "aquamarine"
legend_text_size <- 0.7
point_size <- 0.5
legend_title_size <- 0.9
limsnll <- c(min(nll1,nll3),max(nll1,nll3))
limsa <- c(min(a1,a3),max(a1,a3))
nrow_facet <- 1
estsf <- cbind(est_all_sf %>% filter(cond_site %in% cond_site_name) %>% dplyr::select(c()), tmp %>% add_row(.before=cond_index))
p1 <- tm_shape(estsf) + tm_dots(fill="nll1",fill.scale = tm_scale_continuous(limits=limsnll,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="NLL")) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Model 1") 
p3 <- tm_shape(estsf) + tm_dots(fill="nll3",fill.scale = tm_scale_continuous(limits=limsnll,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="NLL")) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Model 3") 
tmap_save(tmap_arrange(p1,p3,ncol=3),filename=paste0(folder_name,"nll_compare_map_",cond_site_name,".png"),height=6,width=9)

p1 <- tm_shape(estsf) + tm_dots(fill="a1",fill.scale = tm_scale_continuous(limits=limsa,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="AIC")) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Model 1") 
p3 <- tm_shape(estsf) + tm_dots(fill="a3",fill.scale = tm_scale_continuous(limits=limsa,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="AIC")) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Model 3") 
tmap_save(tmap_arrange(p1,p3,ncol=3),filename=paste0(folder_name,"AIC_compare_map_",cond_site_name,".png"),height=6,width=9)

# 5. map and boxplot of AIC difference ----------------------------------------
diff13=a1-a3
tmp2 <- tmp %>% mutate(diff13)
limsad <- c(min(diff13,max(diff13)))
estsf <- cbind(est_all_sf %>% filter(cond_site %in% cond_site_name) %>% dplyr::select(c()), tmp2 %>% add_row(.before=cond_index))
p1 <- tm_shape(estsf) + tm_dots(fill="diff13",fill.scale = tm_scale_continuous(limits=limsad,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="AIC difference")) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Model 1 - Model 3") 
tmap_save(tmap_arrange(p1,ncol=1),filename=paste0(folder_name,"AIC_difference_map_",cond_site_name,".png"),height=6,width=3)
# boxplot of AIC difference
tmp3 <- tmp2 %>% dplyr::select(contains("diff")) %>% pivot_longer(cols=everything())
p <- ggplot(tmp3 %>% mutate(name=factor(name,levels=c("diff13")))) + geom_boxplot(aes(y=value,x=name)) + labs(y="AIC difference",x="") + scale_x_discrete(labels=c("Model 1 - Model 3"))
ggsave(p,filename=paste0(folder_name,"AIC_compare_diff_boxplot_",cond_site_name,".png"),width=9,height=5)

# 6. AIC difference positive and negative -------------------------------------
tmp4 <- tmp2 %>% mutate("is_big"=(diff13>0))
lims13 <- c(0,max(tmp4$diff13))
estsf <- cbind(est_all_sf %>% filter(cond_site %in% cond_site_name) %>% dplyr::select(c()), tmp4 %>% add_row(.before=cond_index))
p1 <- tm_shape(estsf) + tm_dots(fill="black",size=0.1) +  tm_shape(estsf %>% dplyr::filter(is_big %in% c(TRUE,NA))) + tm_dots(fill="diff13",fill.scale = tm_scale_continuous(values="viridis",value.na=misscol,limits = lims13,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="AIC")) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Model 1 - Model 3") 
p2 <- tm_shape(estsf) + tm_dots(fill="is_big",fill.scale = tm_scale_categorical(values=c("black","#C11432"),value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="")) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Is Model 3 better?") 
tmap_save(tmap_arrange(p2,p1,ncol=2),filename=paste0(folder_name,"AIC_difference_1_3_",cond_site_name,".png"),height=6,width=6)

tmp4 <- tmp2 %>% mutate("is_big"=(diff13<0))
lims31 <- c(-max(tmp4$diff13),0)
estsf <- cbind(est_all_sf %>% filter(cond_site %in% cond_site_name) %>% dplyr::select(c()), tmp4 %>% add_row(.before=cond_index))
p1 <- tm_shape(estsf) + tm_dots(fill="black",size=0.1) +  tm_shape(estsf %>% dplyr::filter(is_big %in% c(TRUE,NA))) + tm_dots(fill="diff13",fill.scale = tm_scale_continuous(values="-viridis",value.na=misscol,limits=lims31,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="AIC",reverse=TRUE)) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Model 1 - Model 3") 
p2 <- tm_shape(estsf) + tm_dots(fill="is_big",fill.scale = tm_scale_categorical(values=c("black","#C11432"),value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="")) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Is Model 1 better?") 
tmap_save(tmap_arrange(p2,p1,ncol=2),filename=paste0(folder_name,"AIC_difference_1_3_neg_",cond_site_name,".png"),height=6,width=6)

# 7. map AIC difference above certain values ----------------------------------
howbetter1 <- 20
howbetter2 <- 100
tmp4 <- tmp2 %>% mutate("is_big1"=(diff13>howbetter1),"is_big2"=(diff13>howbetter2))
estsf <- cbind(est_all_sf %>% filter(cond_site %in% cond_site_name) %>% dplyr::select(c()), tmp4 %>% add_row(.before=cond_index))
p1 <- tm_shape(estsf) + tm_dots(fill="black",size=0.1) +  tm_shape(estsf %>% dplyr::filter(is_big1 %in% c(TRUE,NA))) + tm_dots(fill="diff13",fill.scale = tm_scale_continuous(values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="AIC difference")) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="AIC difference > 20") 
p2 <- tm_shape(estsf) + tm_dots(fill="black",size=0.1) +  tm_shape(estsf %>% dplyr::filter(is_big2 %in% c(TRUE,NA))) + tm_dots(fill="diff13",fill.scale = tm_scale_continuous(values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="AIC difference")) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="AIC difference > 100") 
tmap_save(tmap_arrange(p1,p2,ncol=2),filename=paste0(folder_name,"AIC_difference_1_3_",howbetter1,"_",howbetter2,"_AICdiff_",cond_site_name,".png"),height=6,width=6)

library(tidyverse)
library(latex2exp)
library(gridExtra)
library(sf)
library(tmap)

theme_set(theme_bw())
theme_replace(
  panel.spacing = unit(2, "lines"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA) )
load("data_processed/spatial_helper.RData", verbose = TRUE)
folder_name <- "../Documents/final_model_diagnostics/"
# load data to illustrate plot
q <- 0.9 # set quantile threshold
load(paste0("data_processed/N9000_sequential2_AGG_all12sites",q*100,".RData"),verbose = TRUE) # original parameter estimates
# load model 3 estimates
load("data_processed/final_model_3_parameter_estimates.RData",verbose=TRUE)

# 1. diagnostic plots from final iterative model ------------------------------
plot_diagnostics_final_model <- function(y_mod3,xyUK20_sf,site_i) {
folder_name <- "../Documents/final_model_3/"
cond_index <- y_mod3[[1]]$given[1]
cond_site_name <- names(df_sites)[site_i]
# 1. map alpha and beta new and original estimates ----------------------------
a_orig <- append(y_mod3[[1]]$a,NA,after=cond_index-1)
b_orig <- append(y_mod3[[1]]$b,NA,after=cond_index-1)
a_new <- y_mod3[[3]] %>% dplyr::filter(iteration==max(iteration)) %>% dplyr::select(a) %>% pull(a)
b_new <- y_mod3[[3]] %>% dplyr::filter(iteration==max(iteration)) %>% dplyr::select(b)  %>% pull(b)
title_map <- ""
misscol <- "aquamarine"
  legend_text_size <- 0.7
  point_size <- 0.5
  legend_title_size <- 0.9
  limsa <- c(0,1)
  limsb <- c(0,0.75)
  nrow_facet <- 1
  estsf <- cbind(xyUK20_sf, data.frame("a"=a_orig,"b"=b_orig,"a_new"=a_new,"b_new" = b_new))
  p1 <- tm_shape(estsf) + tm_dots(fill="a",fill.scale = tm_scale_continuous(limits=limsa,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\alpha$"))) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text=TeX("$\\tilde{\\alpha}$")) 
  p2 <- tm_shape(estsf) + tm_dots(fill="b",fill.scale = tm_scale_continuous(limits=limsb,values="Blues",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\beta$"))) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text=TeX("$\\tilde{\\beta}$")) 
  p3 <- tm_shape(estsf) + tm_dots(fill="a_new",fill.scale = tm_scale_continuous(limits=limsa,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\alpha$"))) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text=TeX("$\\hat{\\alpha}$")) 
  p4 <- tm_shape(estsf) + tm_dots(fill="b_new",fill.scale = tm_scale_continuous(limits=limsb,values="Blues",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\beta$"))) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text=TeX("$\\hat{\\beta}$")) 
  tmap_save(tmap_arrange(p1,p2,p3,p4,ncol=4),filename=paste0(folder_name,"alpha_beta_fixed_res_",cond_site_name,".png"),height=6,width=11)
  
  # 2. scatterplot comparing alpha and beta original and new estimates ----------
  tmp1 <- rbind(data.frame(a=na.omit(a_orig),b=na.omit(b_orig),"method"="aoriginal"),data.frame(a=na.omit(a_new),b=na.omit(b_new),"method"="new")) %>% mutate("iteration"=rep(1:length(na.omit(a_orig)),2))
  plot_ab <- function(tmp) { ggplot(tmp) + 
      geom_line(aes(x=a,y=b,group=iteration),linewidth=0.1) +
      geom_point(aes(x=a,y=b,col=method),alpha=0.7,size=1) +
      xlab(TeX("${\\alpha}$")) +
      ylab(TeX("${\\beta}$")) + 
      scale_color_manual(values = c("aoriginal" = "#009ADA", "new" = "#C11432"),labels = c("Model 1","Model 3")) + coord_fixed() + theme(axis.text.y = element_text(angle = 90, vjust = 0.5)) + labs(col="")
  }
  p <- plot_ab(tmp=tmp1)
  ggsave(p,filename=paste0(folder_name,"plot_ab_new_original_",cond_site_name,".png"),width=7,height=4) 
  
  # 3. NLL comparing original and new estimates
  pe <- y_mod3[[2]] %>% dplyr::filter(iteration==max(iteration)) %>% dplyr::select(sigl,sigu,deltal,deltau)
  abmu <- y_mod3[[3]] %>% dplyr::filter(iteration==max(iteration)) %>% dplyr::select(a,b,mu) %>% na.omit()
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
  pe_orig <- y_mod3[[1]] %>% dplyr::select(a,b,mu,sig) %>% na.omit()
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
    estsf <- cbind(xyUK20_sf %>% dplyr::select(c()), tmp %>% add_row(.before=cond_index))
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
    estsf <- cbind(xyUK20_sf %>% dplyr::select(c()), tmp2 %>% add_row(.before=cond_index))
    p1 <- tm_shape(estsf) + tm_dots(fill="diff13",fill.scale = tm_scale_continuous(limits=limsad,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="AIC difference")) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Model 1 - Model 3") 
    tmap_save(tmap_arrange(p1,ncol=1),filename=paste0(folder_name,"AIC_difference_map_",cond_site_name,".png"),height=6,width=3)
    # boxplot of AIC difference
    tmp3 <- tmp2 %>% dplyr::select(contains("diff")) %>% pivot_longer(cols=everything())
    p <- ggplot(tmp3 %>% mutate(name=factor(name,levels=c("diff13")))) + geom_boxplot(aes(y=value,x=name)) + labs(y="AIC difference",x="") + scale_x_discrete(labels=c("Model 1 - Model 3"))
    ggsave(p,filename=paste0(folder_name,"AIC_compare_diff_boxplot_",cond_site_name,".png"),width=9,height=5)
    
    # 6. AIC difference positive and negative -------------------------------------
    tmp4 <- tmp2 %>% mutate("is_big"=(diff13>0))
    lims13 <- c(0,max(tmp4$diff13))
    estsf <- cbind(xyUK20_sf %>% dplyr::select(c()), tmp4 %>% add_row(.before=cond_index))
    p1 <- tm_shape(estsf) + tm_dots(fill="black",size=0.1) +  tm_shape(estsf %>% dplyr::filter(is_big %in% c(TRUE,NA))) + tm_dots(fill="diff13",fill.scale = tm_scale_continuous(values="viridis",value.na=misscol,limits = lims13,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="AIC")) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Model 1 - Model 3") 
    p2 <- tm_shape(estsf) + tm_dots(fill="is_big",fill.scale = tm_scale_categorical(values=c("FALSE"="black","TRUE"="#C11432"),value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="")) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Is Model 3 better?") 
    tmap_save(tmap_arrange(p2,p1,ncol=2),filename=paste0(folder_name,"AIC_difference_1_3_",cond_site_name,".png"),height=6,width=6)
    
    tmp4 <- tmp2 %>% mutate("is_big"=(diff13<0))
    lims31 <- c(-max(tmp4$diff13),0)
    estsf <- cbind(xyUK20_sf %>% dplyr::select(c()), tmp4 %>% add_row(.before=cond_index))
    p1 <- tm_shape(estsf) + tm_dots(fill="black",size=0.1) +  tm_shape(estsf %>% dplyr::filter(is_big %in% c(TRUE,NA))) + tm_dots(fill="diff13",fill.scale = tm_scale_continuous(values="-viridis",value.na=misscol,limits=lims31,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="AIC",reverse=TRUE)) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Model 1 - Model 3") 
    p2 <- tm_shape(estsf) + tm_dots(fill="is_big",fill.scale = tm_scale_categorical(values=c("FALSE"="black","TRUE"="#C11432"),value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="")) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Is Model 1 better?") 
    tmap_save(tmap_arrange(p2,p1,ncol=2),filename=paste0(folder_name,"AIC_difference_1_3_neg_",cond_site_name,".png"),height=6,width=6)
    
    # 7. map AIC difference above certain values ----------------------------------
    howbetter1 <- 20
    howbetter2 <- 100
    tmp4 <- tmp2 %>% mutate("is_big1"=(diff13>howbetter1),"is_big2"=(diff13>howbetter2))
    estsf <- cbind(xyUK20_sf %>% dplyr::select(c()), tmp4 %>% add_row(.before=cond_index))
    p1 <- tm_shape(estsf) + tm_dots(fill="black",size=0.1) +  tm_shape(estsf %>% dplyr::filter(is_big1 %in% c(TRUE,NA))) + tm_dots(fill="diff13",fill.scale = tm_scale_continuous(values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="AIC difference")) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="AIC difference > 20") 
    p2 <- tm_shape(estsf) + tm_dots(fill="black",size=0.1) +  tm_shape(estsf %>% dplyr::filter(is_big2 %in% c(TRUE,NA))) + tm_dots(fill="diff13",fill.scale = tm_scale_continuous(values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="AIC difference")) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="AIC difference > 100") 
    tmap_save(tmap_arrange(p1,p2,ncol=2),filename=paste0(folder_name,"AIC_difference_1_3_",howbetter1,"_",howbetter2,"_AICdiff_",cond_site_name,".png"),height=6,width=6)
}
sapply(1:13,FUN=function(i){plot_diagnostics_final_model(y_mod3 = par_est_model_3[[i]],xyUK20_sf = xyUK20_sf,site_i=i)})

# 2. plot of beta against latitude difference from the conditioning site ------
plot_beta_latitude <- function(b,given,res,gridUK=xyUK20_sf,cond_site,east_noeast=NULL,folder_name,plot_name="beta_latitude_difference",comb_sites=FALSE,line_mean=FALSE) {
  stopifnot(length(b)==length(given))
lat_diff <-   sapply(1:length(b),FUN=function(i) {
    as.numeric(gridUK$lat)[res[i]] - as.numeric(gridUK$lat)[given[i]]
  })
if (is.null(east_noeast)) {
  tmp <- data.frame("b"=b,"cond_site"=cond_site,"lat_diff"=lat_diff)  
} else {
tmp <- data.frame("b"=b,"cond_site"=cond_site,"lat_diff"=lat_diff,"east_noeast"=east_noeast)
}
c13 <- c(
  "#009ADA", "#C11432", # red
           "green4",
           "#6A3D9A", # purple
           "#FF7F00", # orange
           "black", "gold1",
           "#FB9A99", # lt pink
           "gray70", 
           "darkturquoise", "green1", 
           "darkorange4","#F6A7B8"
)
if (comb_sites==FALSE) {
p <- ggplot(tmp) + geom_point(aes(x=lat_diff,y=b,col=cond_site),size=0.3) + facet_wrap(cond_site) + scale_color_manual(values=c13) + labs(x="Latitude difference",y=TeX("$\\beta$"),col="Conditioning site") + theme(axis.title.y = element_text(angle = 0,vjust=0.5)) 
ggsave(p,filename=paste0(folder_name,plot_name,".png"),width=10,height=7)
}
if (comb_sites==TRUE) {
  p <- ggplot(tmp) + geom_point(aes(x=lat_diff,y=b,col=cond_site),size=0.3) + scale_color_manual(values=c13) + labs(x="Latitude difference",y=TeX("$\\beta$"),col="Conditioning site") + theme(axis.title.y = element_text(angle = 0,vjust=0.5)) 
  if (line_mean==TRUE) {
    p <- p + geom_smooth(aes(y=b,x=lat_diff,col=cond_site),method="loess")
  }
  if (line_mean=="all") {
    p <- ggplot(tmp) + geom_point(aes(x=lat_diff,y=b),size=0.3) + labs(x="Latitude difference",y=TeX("$\\beta$"),col="Conditioning site") + 
      theme(axis.title.y = element_text(angle = 0,vjust=0.5)) + geom_smooth(aes(y=b,x=lat_diff,col=east_noeast),method="loess") +
      geom_smooth(aes(y=b,x=lat_diff,col="black"),method="loess") +
      scale_color_manual(values=c("#009ADA","#C11432","black"),labels=c("East coast", "Not east coast","Combined"))
  }
  ggsave(p,filename=paste0(folder_name,plot_name,".png"),width=5,height=4)
}
}
# first draft apply to the original estimates
#plot_beta_latitude(b=est_all_sf$b,given=est_all_sf$given,res=est_all_sf$res,cond_site = est_all_sf$cond_site,folder_name = folder_name)

tmp <- data.frame("b_new"=numeric(),"b_old"=numeric(),"given"=numeric(),"res"=numeric(),"cond_site"=character())
for (i in 1:length(par_est_model_3)) {
  tmp <- rbind(tmp,data.frame("b_new"=par_est_model_3[[i]][[3]] %>% filter(iteration==10) %>% pull(b) %>% na.omit(),
                              "b_old"=par_est_model_3[[i]][[1]] %>% pull(b),
                              "given"=par_est_model_3[[i]][[1]] %>% pull(given),
                              "res"=par_est_model_3[[i]][[1]] %>% pull(res),
                              "cond_site" = names(df_sites)[i]))
}


tmp <- tmp %>% mutate(cond_site=factor(cond_site))
coastal_point <- function(grid) {
  sapply(1:nrow(grid),FUN = function(i) {sum(as.vector(st_distance(grid[i,],grid))<20500)<5 & grid$lat[i]+2*grid$lon[i]>48.5})
}
cp <- coastal_point(grid = xyUK20_sf) 
cp[df_sites[3,5]] <- FALSE
east_coast <-   sapply(1:nrow(tmp),FUN=function(i) {
  cp[tmp$res[i]]
})
tmp <- tmp %>% mutate("east_coast"=east_coast)
tmp <- tmp %>% dplyr::filter(east_coast=="TRUE")
plot_beta_latitude(b=tmp$b_new,given=tmp$given,res=tmp$res,cond_site = tmp$cond_site,folder_name = folder_name,plot_name = "beta_model_3")
plot_beta_latitude(b=tmp$b_old,given=tmp$given,res=tmp$res,cond_site = tmp$cond_site,folder_name = folder_name,plot_name = "beta_model_1")

plot_beta_latitude(b=tmp$b_new,given=tmp$given,res=tmp$res,cond_site = tmp$cond_site,folder_name = folder_name,plot_name = "beta_model_3_comb",comb_sites=TRUE)

east_coast_sites <- c("Newcastle","Lowestoft","Cromer","Hull","Inverness")
tmp_ec <- tmp %>% filter(cond_site %in% east_coast_sites)
tmp_nec <- tmp %>% filter(!cond_site %in% east_coast_sites)
plot_beta_latitude(b=tmp_ec$b_new,given=tmp_ec$given,res=tmp_ec$res,cond_site = tmp_ec$cond_site,folder_name = folder_name,plot_name = "beta_model_3_east_coast",comb_sites=TRUE)
plot_beta_latitude(b=tmp_nec$b_new,given=tmp_nec$given,res=tmp_nec$res,cond_site = tmp_nec$cond_site,folder_name = folder_name,plot_name = "beta_model_3_not_east_coast",comb_sites=TRUE)

# add regression lines
plot_beta_latitude(b=tmp_ec$b_new,given=tmp_ec$given,res=tmp_ec$res,cond_site = tmp_ec$cond_site,folder_name = folder_name,plot_name = "beta_model_3_east_coast_mean_smooth",comb_sites=TRUE,line_mean = TRUE)
plot_beta_latitude(b=tmp_nec$b_new,given=tmp_nec$given,res=tmp_nec$res,cond_site = tmp_nec$cond_site,folder_name = folder_name,plot_name = "beta_model_3_not_east_coast_mean_smooth",comb_sites=TRUE,line_mean=TRUE)

# combine into one plot
tmp <- tmp %>% mutate("east_noeast" = factor(ifelse(cond_site %in% east_coast_sites,"east_coast","not_east_coast")))
plot_beta_latitude(b=tmp$b_new,given=tmp$given,res=tmp$res,cond_site = tmp$cond_site,east_noeast=tmp$east_noeast,folder_name = folder_name,plot_name = "beta_model_3_all_mean_smooth",comb_sites=TRUE,line_mean="all")

# fit a beta distribution to beta values on the east coast for Birmingham
NLL_beta_beta <- function(theta,a=NULL,b=NULL,beta=NULL,data_Lap,d_latitude,mu,phi,deltal,deltau,dij,alpha,cond_index,east_coast,res) {
print(theta)
data_Lapv <- data_Lap %>% filter(data_Lap[,cond_index]>quantile(data_Lap[,cond_index],0.9))
#data_Lapv1<- data_Lapv[,c(east_coast,tmp$cond_index[1])]

c <- theta[1]
#a <- theta[2]
#b <- theta[3]
gamma1 <- theta[2]
gamma2 <- theta[3]
# if (gamma1<0 | gamma2<0 | c<0 | min((d_latitude-a)/(b-a))<0 | max((d_latitude-a)/(b-a))>1) {return(10e10)}
if (gamma1<0 | gamma2<0 | c<0 ) {return(10e10)}
else {
  if (is.null(beta)) {
beta <- c*stats::dbeta(x=(d_latitude-a)/(b-a),shape1=gamma1,shape2=gamma2)
  }
sigu <- phi[1] + phi[2]*(1-exp(-(phi[3]*dij)))
sigl <- phi[4] + phi[5]*(1-exp(-(phi[6]*dij)))
y_res <- sapply(1:length(res),FUN=function(j) {NLL_AGG_onestep(theta=c(),x=data.frame("Y1"=data_Lapv[,cond_index],"Y2"=data_Lapv[,res[j]]),a_hat=alpha[j],b_hat=beta[j],mu_hat = mu[j],sigl_hat = sigl[j],sigu_hat = sigu[j],deltal_hat = deltal[j],deltau_hat = deltau[j]) } )
print(y_res)
y <- sum(y_res)
return(y_res)
}
}

est_beta <- function(beta=NULL,data_Lap,d_latitude,mu,phi,deltal,deltau,dij,alpha,cond_index,east_coast,res) {

  data_Lapv <- data_Lap %>% filter(data_Lap[,cond_index]>quantile(data_Lap[,cond_index],0.9))
  #data_Lapv1<- data_Lapv[,c(east_coast,tmp$cond_index[1])]
    sigu <- phi[1] + phi[2]*(1-exp(-(phi[3]*dij)))
    sigl <- phi[4] + phi[5]*(1-exp(-(phi[6]*dij)))
    beta <- sapply(1:length(sigl),FUN=function(j){
      x <- optim(par=c(0.2),fn=NLL_AGG_onestep,x=data.frame("Y1"=data_Lapv[,cond_index],"Y2"=data_Lapv[,res[j]]),a_hat=alpha[j],b_hat=beta[j],mu_hat = mu[j],sigl_hat = sigl[j],sigu_hat = sigu[j],deltal_hat = deltal[j],deltau_hat = deltau[j])
      return(x)})
    print(beta)
    
    return(beta)
}

beta_model_prepare_dataset <- function(data,gridUK=xyUK20_sf,site_i,df_sites) {
  tmp <- data.frame("alpha"=numeric(),"beta"=numeric(),"given"=numeric(),"res"=numeric(),"cond_site"=character(),
                    "mu"=numeric(),"sigl"=numeric(),"sigu"=numeric(),"deltal"=numeric(),"deltau"=numeric()
                    )
    tmp <- rbind(tmp,data.frame("alpha"=as.numeric(data[[site_i]][[3]] %>% dplyr::filter(iteration==max(iteration)) %>% pull(a) %>% na.omit()),
                                "beta"=as.numeric(data[[site_i]][[3]] %>% dplyr::filter(iteration==max(iteration)) %>% pull(b) %>% na.omit()),
                                "given"=data[[site_i]][[1]] %>% pull(given),
                                "res"=data[[site_i]][[1]] %>% pull(res),
                                "cond_site" = names(df_sites)[site_i],
                                "mu"=data[[site_i]][[2]] %>% dplyr::filter(iteration==max(iteration)) %>% pull(mu_agg),   
                                "sigl"=data[[site_i]][[2]] %>% dplyr::filter(iteration==max(iteration)) %>% pull(sigl),   
                                "sigu"=data[[site_i]][[2]] %>% dplyr::filter(iteration==max(iteration)) %>% pull(sigu),   
                                "deltal"=data[[site_i]][[2]] %>% dplyr::filter(iteration==max(iteration)) %>% pull(deltal),   
                                "deltau"=data[[site_i]][[2]] %>% dplyr::filter(iteration==max(iteration)) %>% pull(deltau)  ))

  # calculate additional variables
  cond_index <- tmp$given[1]
  # distance from the conditioning site
  dist_tmp <- as.numeric(unlist(st_distance(gridUK[cond_index,],gridUK[-cond_index,])))
  distnorm <- dist_tmp/1000000  # normalise distance using a common constant
  # latitude difference from the conditioning site
  lat_diff <-   sapply(1:nrow(tmp),FUN=function(k) {
    as.numeric(gridUK$lat)[tmp$res[k]] - as.numeric(gridUK$lat)[tmp$given[k]]
  })
  tmp <- tmp %>% mutate("dij"=distnorm,"d_latitude"=lat_diff)
  return(tmp)
}

# test new function
site_i <- 1
y1 <- beta_model_prepare_dataset(data=par_est_model_3,site_i=site_i,df_sites=df_sites)
# reestimate all betas
xb_all <- est_beta(data_Lap=data_mod_Lap,phi = phis,east_coast=east_coast,res=y1$res,mu=y1$mu,deltal=y1$deltal,deltau=y1$deltau,alpha=y1$alpha,dij=y1$dij,d_latitude = y1$d_latitude,cond_index=y1$given[1])
bjall <- likball <- c()
for (i in 1:nrow(y1)) {
  bjall <- append(bjall,xb_all[,i]$par)
  likball <- append(likball,xb_all[,i]$value)
}
title_map <- ""
misscol <- "aquamarine"
  legend_text_size <- 0.7
  point_size <- 0.5
  legend_title_size <- 0.9
  limsa <- c(0,1)
  limsb <- c(0,0.75)
  nrow_facet <- 1
  
estsf <- cbind(xyUK20_sf, data.frame("a"=y1$alpha,"b"=y1$beta,"b_new" = bjall) %>% add_row(.before = y1$given[1]))
p1 <- tm_shape(estsf) + tm_dots(fill="a",fill.scale = tm_scale_continuous(limits=limsa,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\alpha$"))) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text=TeX("$\\tilde{\\alpha}$")) 
p2 <- tm_shape(estsf) + tm_dots(fill="b",fill.scale = tm_scale_continuous(limits=limsb,values="Blues",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\beta$"))) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text=TeX("$\\tilde{\\beta}$")) 
p4 <- tm_shape(estsf) + tm_dots(fill="b_new",fill.scale = tm_scale_continuous(limits=limsb,values="Blues",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\beta$"))) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text=TeX("$\\hat{\\beta}$")) 
tmap_save(tmap_arrange(p1,p2,p4,ncol=4),filename=paste0(folder_name,"new_alpha_beta_fixed_res_","Birmingham",".png"),height=6,width=11)


phis <- as.numeric(((par_est_model_3[[site_i]][[2]] %>% dplyr::select(phi0u,phi1u,phi2u,phi0l,phi1l,phi2l))[1,]))
# subset for only east coast
east_coast <-   sapply(1:nrow(y),FUN=function(j_ec) {
  cp[y$res[j_ec]]
})
y <- y1[east_coast,]
y <- y %>% dplyr::filter(beta>0.001)
#x <- optim(par = c(0.5,-4,4,2,2),fn=NLL_beta_beta,data_Lap=data_mod_Lap,phi = phis,east_coast=east_coast,res=y$res,mu=y$mu,deltal=y$deltal,deltau=y$deltau,alpha=y$alpha,dij=y$dij,d_latitude = y$d_latitude,cond_index=y$given[1],control = list(maxit=2000))
a <- -2
b <- 4
x <- optim(par = c(0.5,2,2),fn=NLL_beta_beta,data_Lap=data_mod_Lap,a=a,b=b,phi = phis,east_coast=east_coast,res=y$res,mu=y$mu,deltal=y$deltal,deltau=y$deltau,alpha=y$alpha,dij=y$dij,d_latitude = y$d_latitude,cond_index=y$given[1],control = list(maxit=2000))
# c <- x$par[1]
# a <- x$par[2]
# b <- x$par[3]
# gamma1 <- x$par[4]
# gamma2 <- x$par[5]
c <- x$par[1]
gamma1 <- x$par[2]
gamma2 <- x$par[3]
# compare for different values of c
c_seq <- seq(0.01,0.5,length.out=50)
#c_likb <- sapply(1:length(c_seq),FUN = function(i) {NLL_beta_beta(theta = c(c_seq[i],gamma1,gamma2),beta=y$beta,data_Lap=data_mod_Lap,a=a,b=b,phi = phis,east_coast=east_coast,res=y$res,mu=y$mu,deltal=y$deltal,deltau=y$deltau,alpha=y$alpha,dij=y$dij,d_latitude = y$d_latitude,cond_index=y$given[1])})
#c_likb <- sapply(1:length(c_seq),FUN = function(i) {NLL_beta_beta(theta = c(c_seq[i],gamma1,gamma2),beta=y$beta,data_Lap=data_mod_Lap,a=a,b=b,phi = phis,east_coast=east_coast,res=y$res,mu=y$mu,deltal=y$deltal,deltau=y$deltau,alpha=y$alpha,dij=y$dij,d_latitude = y$d_latitude,cond_index=y$given[1])})
ybeta <- NLL_beta_beta(theta = c(c_seq[i],gamma1,gamma2),beta=y$beta,data_Lap=data_mod_Lap,a=a,b=b,phi = phis,east_coast=east_coast,res=y$res,mu=y$mu,deltal=y$deltal,deltau=y$deltau,alpha=y$alpha,dij=y$dij,d_latitude = y$d_latitude,cond_index=y$given[1])
yc <- NLL_beta_beta(theta = c(c,gamma1,gamma2),beta=NULL,data_Lap=data_mod_Lap,a=a,b=b,phi = phis,east_coast=east_coast,res=y$res,mu=y$mu,deltal=y$deltal,deltau=y$deltau,alpha=y$alpha,dij=y$dij,d_latitude = y$d_latitude,cond_index=y$given[1])
p <- ggplot(data.frame(ybeta,yc)) + geom_point(aes(x=ybeta,y=yc)) + geom_abline(slope=1,intercept=0,linetype="dashed") + labs(x="Model 3 beta NLLs",y="Parametric beta function NLLs")
plot_name <- "different_beta_function_compare_separate_beta"
ggsave(p,filename=paste0(folder_name,plot_name,".png"),width=5,height=4)

p <- ggplot(data.frame(c_seq,c_lik) %>% dplyr::
              filter(c_lik<10e10)) + geom_point(aes(x=c_seq,y=c_lik))
plot_name <- "different_beta_function_c_likelihood_b_4_bigb_y"
ggsave(p,filename=paste0(folder_name,plot_name,".png"),width=5,height=4)


beta_new <- c*stats::dbeta(x=(latseq-a)/(b-a),shape1=gamma1,shape2=gamma2) 
# reestimate beta separately at each site
xb <- est_beta(data_Lap=data_mod_Lap,phi = phis,east_coast=east_coast,res=y$res,mu=y$mu,deltal=y$deltal,deltau=y$deltau,alpha=y$alpha,dij=y$dij,d_latitude = y$d_latitude,cond_index=y$given[1])
bj <- likb <- c()
for (i in 1:nrow(y)) {
  bj <- append(bj,xb[,i]$par)
  likb <- append(likb,xb[,i]$value)
}

# compare Model 3, Model 3 separate and parametric form
latseq <- seq(-2,7,length.out=500)
beta_par <- c*stats::dbeta(x=(latseq-a)/(b-a),shape1=gamma1,shape2=gamma2) 
p <- ggplot() + geom_point(data=y %>% mutate("beta"=bj) ,aes(x=d_latitude,y=beta),fill="#009ADA") + geom_line(data=data.frame("latseq"=latseq,"beta"=beta_par),aes(x=latseq,y=beta),color="#C11432") + labs(x="Latitude difference",y=TeX("$\\beta$")) + 
  scale_color_manual(values=c10) +  theme(axis.title.y = element_text(angle = 0,vjust=0.5))
#plot_name <- "different_beta_function_initial_values_SANN_fix_ab_gamma2"
plot_name <- "different_beta_function_mod3_revisited"
ggsave(p,filename=paste0(folder_name,plot_name,".png"),width=5,height=4)


# plot Model 3 and Model 3 estimates
p <- ggplot(data.frame("y_mod3"=y$beta,"y_mod3_sep"=bj)) + geom_point(aes(x=y_mod3,y=y_mod3_sep)) + geom_abline(slope=1,intercept=0,linetype="dashed") + labs(x="Model 3 beta NLLs",y="Reestimated beta NLLs")
plot_name <- "model_3_compare_separate_beta"
ggsave(p,filename=paste0(folder_name,plot_name,".png"),width=5,height=4)

p <- ggplot(data.frame(ybeta,yc,likb)) + geom_point(aes(x=ybeta,y=yc)) + geom_point(aes(x=likb,y=yc),color="#C11432") + geom_abline(slope=1,intercept=0,linetype="dashed") + labs(x="Model 3 beta NLLs",y="Parametric beta function NLLs")
plot_name <- "different_beta_function_compare_separate_beta"
ggsave(p,filename=paste0(folder_name,plot_name,".png"),width=5,height=4)

# if (gamma1<0 | gamma2<0 | c<0 | min((d_latitude-a)/(b-a))<0 | max((d_latitude-a)/(b-a))>1) {return(10e10)}
 # beta_new <- c*stats::dbeta(x=(y$d_latitude-a)/(b-a),shape1=gamma1,shape2=gamma2)
beta_new <- c*stats::dbeta(x=(y$d_latitude-a)/(b-a),shape1=gamma1,shape2=gamma2)
p <- seq(0,1,length.out=10000)
plot(x=p,y=y_result)
# plot fitted density and observed
ggplot(y %>% mutate("b_new"=beta_new)) + geom_point(aes(x=d_latitude,y=beta),fill="#009ADA") + geom_line(aes(x=d_latitude,y=b_new))

# check whether initial values matter
latseq <- seq(-2,7,length.out=500)
init_beta <- data.frame(matrix(data=c(0.05,-4,4,2,2,0.1,-4,4,2,2,0.5,-4,4,2,2,0.5,-2,2,2,2,0.5,-2,4,2,2,0.5,-2,6,2,2,0.5,0,1,2,2,0.1,-4,4,2,1,0.1,-4,4,1,2,0.1,-4,4,1,1),byrow=TRUE,ncol=5))
names(init_beta) <- c("c","a","b","gamma1","gamma2")
init_beta <- init_beta[,c(1,4,5)]
xpar <- init_beta
xall <- list()
tmpi <- data.frame("latseq"=numeric(),"beta"=numeric(),"beta_parameters"=character())
for (i in 1:nrow(init_beta)) {
  # x <- optim(par = as.numeric(init_beta[i,]),fn=NLL_beta_beta,data_Lap=data_mod_Lap,phi = phis,east_coast=east_coast,res=y$res,mu=y$mu,deltal=y$deltal,deltau=y$deltau,alpha=y$alpha,dij=y$dij,d_latitude = y$d_latitude,cond_index=y$given[1],control = list(maxit=2000),method="SANN")
  # c <- x$par[1]
  # a <- x$par[2]
  # b <- x$par[3]
  # gamma1 <- x$par[4]
  # gamma2 <- x$par[5]
  # a <- -2
  # b <- 4
 x <- optim(par =as.numeric(init_beta[i,]),fn=NLL_beta_beta,data_Lap=data_mod_Lap,a=a,b=b,phi = phis,east_coast=east_coast,res=y$res,mu=y$mu,deltal=y$deltal,deltau=y$deltau,alpha=y$alpha,dij=y$dij,d_latitude = y$d_latitude,cond_index=y$given[1],control = list(maxit=2000))
  
#  x <- optim(par =as.numeric(init_beta[i,]),fn=NLL_beta_beta,data_Lap=data_mod_Lap,a=a,b=b,phi = phis,east_coast=east_coast,res=y$res,mu=y$mu,deltal=y$deltal,deltau=y$deltau,alpha=y$alpha,dij=y$dij,d_latitude = y$d_latitude,cond_index=y$given[1],control = list(maxit=2000),method="SANN")
  
  # c <- x$par[1]
  # a <- x$par[2]
  # b <- x$par[3]
  # gamma1 <- x$par[4]
  # gamma2 <- x$par[5]
  c <- x$par[1]
  gamma1 <- x$par[2]
  gamma2 <- x$par[3]
  
  xpar[i,] <- x$par
  xall[[i]] <- x
  # if (gamma1<0 | gamma2<0 | c<0 | min((d_latitude-a)/(b-a))<0 | max((d_latitude-a)/(b-a))>1) {return(10e10)}
  #beta_new <- c*stats::dbeta(x=(latseq-a)/(b-a),shape1=gamma1,shape2=gamma2) 
  beta_new <- c*stats::dbeta(x=(latseq-a)/(b-a),shape1=gamma1,shape2=gamma2) 
  tmpi <- rbind(tmpi,data.frame("latseq"=latseq,"beta"=beta_new,"beta_parameters"=paste("(",paste(as.character(as.numeric(init_beta[i,])), sep="' '", collapse=", "),")")))
}

# plot
c10 <- c(
  "#009ADA", "#C11432", # red
           "green4",
           "#6A3D9A", # purple
           "#FF7F00", # orange
            "gold1",
           "black", # lt pink
           "gray70", 
           "darkorange4","#F6A7B8"
)
p <- ggplot() + geom_point(data=y ,aes(x=d_latitude,y=beta),fill="#009ADA") + geom_line(data=tmpi %>% mutate(beta_parameters=factor(beta_parameters)),aes(x=latseq,y=beta,col=beta_parameters)) + labs(x="Latitude difference",y=TeX("$\\beta$"),col="Initial parameter values") + 
  scale_color_manual(values=c10) +  theme(axis.title.y = element_text(angle = 0,vjust=0.5))
#plot_name <- "different_beta_function_initial_values_SANN_fix_ab_gamma2"
plot_name <- "different_beta_function_initial_values_Nelder_Mead_fix_ab_gamma2_bigb4"
ggsave(p,filename=paste0(folder_name,plot_name,".png"),width=5,height=4)

# plot on a map
cond_index <- y$given[1]
t1 <- xyUK20_sf %>% mutate("b_model3"=NA,"b_parametric"=NA)
t1$b_model3[-cond_index][east_coast] <- y$beta
t1$b_parametric[-cond_index][east_coast] <- beta_new
t2 <- tm_shape(t1) + tm_dots("b_model3")
t3 <- tm_shape(t1) + tm_dots("b_parametric")
tmap_arrange(t2,t3,ncol=2)

# combine estimation for sites not on the east coast ---------------------------------


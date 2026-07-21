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

folder_name <- "../Documents/final_model_diagnostics/"
# load data to illustrate plot
q <- 0.9 # set quantile threshold
load(paste0("data_processed/N9000_sequential2_AGG_all12sites",q*100,".RData"),verbose = TRUE) # original parameter estimates

#

# 2. plot of beta against latitude difference from the conditioning site ------
plot_beta_latitude <- function(b,given,res,gridUK=xyUK20_sf,cond_site,east_noeast=NULL,folder_name,plot_name="beta_latitude_difference",comb_sites=FALSE,line_mean=FALSE) {
  stopifnot(length(b)==length(given))
lat_diff <-   sapply(1:length(b),FUN=function(i) {
    as.numeric(gridUK$lat)[res[i]] - as.numeric(gridUK$lat)[given[i]]
  })
tmp <- data.frame("b"=b,"cond_site"=cond_site,"lat_diff"=lat_diff,"east_noeast"=east_noeast)

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
      scale_color_manual(values=c("#009ADA","#C11432","black"),labels=c("East coast", "Not East coast","Combined"))
  }
  ggsave(p,filename=paste0(folder_name,plot_name,".png"),width=5,height=4)
}
}

plot_beta_latitude(b=est_all_sf$b,given=est_all_sf$given,res=est_all_sf$res,cond_site = est_all_sf$cond_site,folder_name = folder_name)

# try with new data
load("data_processed/final_model_3_parameter_estimates.RData",verbose=TRUE)
for (i in 1:13) {
  print(par_est_model_3[[i]][[4]])
}

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


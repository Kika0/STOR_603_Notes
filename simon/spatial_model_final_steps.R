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
load("data_processed/iterative_phi0l_phi0u_estimates_London.RData",verbose=TRUE) # residual margin parameters
load("data_processed/residual_dependence_pars.RData", verbose = TRUE) # residual dependence parameters

# 1. reestimate alpha and beta ------------------------------------------------
to_opt <- function(x1,x2,theta,i,cond_index=London_index,pe_i) {
  # a <- theta[1]
  # b <- theta[2]
  # if (a<0 | a>1 | b<0 | b>1) {return(10^6)}
  y <-  NLL_AGG_onestep(x=data.frame(x1,x2),theta=theta,mu_hat=pe_i[1],sigl_hat = pe_i[2], sigu_hat = pe_i[3], deltal_hat = pe_i[4], deltau_hat = pe_i[5])
  return(y)
}
NLL_AGG_wrapper <- function(data_Lap=data_mod_Lap,i,pe_res=pe %>% add_row(.before = London_index),cond_index,v=0.9) {
  pe_i <- as.numeric(unlist(pe_res[i,]))
  data_Lapv <- data_Lap %>% filter(data_Lap[,cond_index]>quantile(data_Lap[,cond_index],v))
  x1 <- as.numeric(unlist(data_Lapv[,London_index]))
  x2 <- as.numeric(unlist(data_Lapv[,i]))
  if (is.na(pe_i[1])) {return(NA)}
  y <- optim(par=c(0.8,0.3),fn= to_opt,x1=x1,x2=x2,i=i,cond_index=cond_index,pe_i=pe_i)
  return(y$par) 
}
x <- sapply(1:ncol(data_mod_Lap),FUN=NLL_AGG_wrapper,data_Lap=data_mod_Lap,cond_index=London_index)
tmp <- as.data.frame(do.call(rbind,x))
names(tmp) <- c("a","b")
a_new <- tmp$a
b_new <- tmp$b
summary(tmp$a)
summary(tmp$b)
# get original parameter values for comparison
aest <- est_all_sf %>% filter(cond_site=="London") %>% pull(a) %>% na.omit()
best <- est_all_sf %>% filter(cond_site=="London") %>% pull(b) %>% na.omit()
muest <- est_all_sf %>% filter(cond_site=="London") %>% pull(mu) %>% na.omit()
sigest <- est_all_sf %>% filter(cond_site=="London") %>% pull(sig) %>% na.omit()

summary(aest)
summary(best)
#to_opt(x1=x1,x2=x2,i=1,theta=c(0.8,0.3),pe_i=pe_i)
tmp1 <- rbind(data.frame(a=aest,b=best,"method"="original"),data.frame(a=na.omit(tmp$a),b=na.omit(tmp$b),"method"="new")) %>% mutate("iteration"=rep(1:length(aest),2))
# map alpha and beta original and new estimates
plot_ab <- function(tmp) { ggplot(tmp) + 
    geom_line(aes(x=a,y=b,group=iteration),linewidth=0.1) +
    geom_point(aes(x=a,y=b,col=method),alpha=0.7,size=1) +
    xlab(TeX("${\\alpha}$")) +
    ylab(TeX("${\\beta}$")) + 
    scale_color_manual(values = c("original" = "#009ADA", "new" = "#C11432")) + coord_fixed() + theme(axis.text.y = element_text(angle = 90, vjust = 0.5))
}
p <- plot_ab(tmp=tmp1)
ggsave(p,filename=paste0(folder_name,"plot_ab_new_original.png"),width=9,height=5)

# map alpha and beta
title_map <- ""
misscol <- "aquamarine"
legend_text_size <- 0.7
point_size <- 0.5
legend_title_size <- 0.9
limsa <- c(0,1)
limsb <- c(0,0.65)
nrow_facet <- 1
estsf <- cbind(est_all_sf %>% filter(cond_site %in% "London"), data.frame("a_new"=a_new,"b_new" = b_new))
p1 <- tm_shape(estsf) + tm_dots(fill="a",fill.scale = tm_scale_continuous(limits=limsa,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\alpha$"))) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text=TeX("$\\tilde{\\alpha}$")) 
p2 <- tm_shape(estsf) + tm_dots(fill="b",fill.scale = tm_scale_continuous(limits=limsb,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\beta$"))) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text=TeX("$\\tilde{\\beta}$")) 
p3 <- tm_shape(estsf) + tm_dots(fill="a_new",fill.scale = tm_scale_continuous(limits=limsa,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\alpha$"))) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text=TeX("$\\hat{\\alpha}$")) 
p4 <- tm_shape(estsf) + tm_dots(fill="b_new",fill.scale = tm_scale_continuous(limits=limsb,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\beta$"))) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text=TeX("$\\hat{\\beta}$")) 
tmap_save(tmap_arrange(p1,p2,p3,p4,ncol=4),filename=paste0(folder_name,"alpha_beta_fixed_res.png"),height=6,width=11)

# 4. Model comparison ---------------------------------------------------------
# model 2b estimates
to_opt <- function(x1,x2,theta,i,cond_index=London_index,pe_i) {
  # a <- theta[1]
  # b <- theta[2]
  # if (a<0 | a>1 | b<0 | b>1) {return(10^6)}
  y <-  NLL_AGG_onestep(x=data.frame(x1,x2),theta=theta,b_hat=0,mu_hat=pe_i[1],sigl_hat = pe_i[2], sigu_hat = pe_i[3], deltal_hat = pe_i[4], deltau_hat = pe_i[5])
  return(y)
}
NLL_AGG_wrapper <- function(data_Lap=data_mod_Lap,i,pe_res=pe %>% add_row(.before = London_index),cond_index,v=0.9) {
  pe_i <- as.numeric(unlist(pe_res[i,]))
  data_Lapv <- data_Lap %>% filter(data_Lap[,cond_index]>quantile(data_Lap[,cond_index],v))
  x1 <- as.numeric(unlist(data_Lapv[,London_index]))
  x2 <- as.numeric(unlist(data_Lapv[,i]))
  if (is.na(pe_i[1])) {return(NA)}
  y <- optim(par=c(0.8),fn= to_opt,x1=x1,x2=x2,i=i,cond_index=cond_index,pe_i=pe_i)
  return(y$par) 
}
x <- sapply(1:ncol(data_mod_Lap),FUN=NLL_AGG_wrapper,data_Lap=data_mod_Lap,cond_index=London_index)
a_new2 <- na.omit(x)
b_new2 <- 0
#to_opt(x1=x1,x2=x2,i=1,theta=c(0.8,0.3),pe_i=pe_i)
tmp1 <- rbind(data.frame(a=a_new2,b=b_new2,"method"="model_2b"),data.frame(a=na.omit(tmp$a),b=na.omit(tmp$b),"method"="model_2")) %>% mutate("iteration"=rep(1:length(aest),2))
# map alpha and beta original and new estimates
plot_ab <- function(tmp) { ggplot(tmp) + 
    geom_line(aes(x=a,y=b,group=iteration),linewidth=0.1) +
    geom_point(aes(x=a,y=b,col=method),alpha=0.7,size=1) +
    xlab(TeX("${\\alpha}$")) +
    ylab(TeX("${\\beta}$")) + 
    scale_color_manual(values = c("model_2b" = "#009ADA", "model_2" = "#C11432")) + coord_fixed() + theme(axis.text.y = element_text(angle = 90, vjust = 0.5))
}
p <- plot_ab(tmp=tmp1)
ggsave(p,filename=paste0(folder_name,"plot_ab_new_original_model2b.png"),width=9,height=5)

# calculate log-likelihood of Model 1
NLL_AGG_wrapper <- function(data_Lap=data_mod_Lap,i,cond_index,v=0.9,a=aest,b=best,mu=muest,sig=sigest) {
  data_Lapv <- data_Lap %>% filter(data_Lap[,cond_index]>quantile(data_Lap[,cond_index],v))
  x1 <- as.numeric(unlist(data_Lapv[,London_index]))
  ires <- c(1:length(data_Lap))[-London_index][i]
  x2 <- as.numeric(unlist(data_Lapv[,ires]))
  y <- NLL_AGG_onestep(x=data.frame(x1,x2),theta=c(),a_hat=a[i],b_hat=b[i],mu_hat=pe_i[1],sigl_hat = sqrt(2)* sig[i], sigu_hat = sqrt(2)* sig[i], deltal_hat = 2, deltau_hat = 2)
  return(y) 
}
nll1 <- sapply(1:(ncol(data_mod_Lap)-1),FUN=NLL_AGG_wrapper,data_Lap=data_mod_Lap,cond_index=London_index)

NLL_AGG_wrapper <- function(data_Lap=data_mod_Lap,i,cond_index,v=0.9,a=aest,b=best,mu=muest,sig=sigest) {
  data_Lapv <- data_Lap %>% filter(data_Lap[,cond_index]>quantile(data_Lap[,cond_index],v))
  x1 <- as.numeric(unlist(data_Lapv[,London_index]))
  ires <- c(1:length(data_Lap))[-London_index][i]
  x2 <- as.numeric(unlist(data_Lapv[,ires]))
y <- Y_likelihood(theta = c(a[i],b[i],mu[i],sig[i]),df = data_Lapv, given = cond_index,sim=ires)
    return(y) 
}
nll1 <- -sapply(1:(ncol(data_mod_Lap)-1),FUN=NLL_AGG_wrapper,data_Lap=data_mod_Lap,cond_index=London_index)

NLL_AGG_wrapper <- function(data_Lap=data_mod_Lap,i,pe_res=pe,cond_index,v=0.9,a=na.omit(a_new),b=na.omit(b_new)) {
  pe_i <- as.numeric(unlist(pe_res[i,]))
  data_Lapv <- data_Lap %>% filter(data_Lap[,cond_index]>quantile(data_Lap[,cond_index],v))
  x1 <- as.numeric(unlist(data_Lapv[,London_index]))
  ires <- c(1:length(data_Lap))[-London_index][i]
  x2 <- as.numeric(unlist(data_Lapv[,ires]))
  y <- NLL_AGG_onestep(x=data.frame(x1,x2),theta=c(),a_hat=a[i],b_hat=b[i],mu_hat=pe_i[1],sigl_hat = pe_i[2], sigu_hat = pe_i[3], deltal_hat = pe_i[4], deltau_hat = pe_i[5])
  return(y) 
}
nll2 <- sapply(1:(ncol(data_mod_Lap)-1),FUN=NLL_AGG_wrapper,data_Lap=data_mod_Lap,cond_index=London_index)

NLL_AGG_wrapper <- function(data_Lap=data_mod_Lap,i,pe_res=pe,cond_index,v=0.9,a=a_new2) {
  pe_i <- as.numeric(unlist(pe_res[i,]))
  data_Lapv <- data_Lap %>% filter(data_Lap[,cond_index]>quantile(data_Lap[,cond_index],v))
  x1 <- as.numeric(unlist(data_Lapv[,London_index]))
  ires <- c(1:length(data_Lap))[-London_index][i]
  x2 <- as.numeric(unlist(data_Lapv[,ires]))
  y <- NLL_AGG_onestep(x=data.frame(x1,x2),theta=c(),a_hat=a[i],b_hat=0,mu_hat=pe_i[1],sigl_hat = pe_i[2], sigu_hat = pe_i[3], deltal_hat = pe_i[4], deltau_hat = pe_i[5])
  return(y) 
}
nll2b <- sapply(1:(ncol(data_mod_Lap)-1),FUN=NLL_AGG_wrapper,data_Lap=data_mod_Lap,cond_index=London_index)

# calculate AIC
a1 <- 2*nll1+2*4
a2 <- 2*nll2+2*3
a2b <- 2*nll2b+2*2
tmp <- data.frame(nll1,nll2,nll2b,a1,a2,a2b)

# boxplots of negative log-likelihood
tmp1 <- tmp %>% dplyr::select(contains("nll")) %>% pivot_longer(cols=everything())
p <- ggplot(tmp1) + geom_boxplot(aes(y=value,x=name)) + labs(y="Negative log-likelihood",x="")
ggsave(p,filename=paste0(folder_name,"nll_compare.png"),width=9,height=5)
# boxplots of AIC
tmp1 <- tmp %>% dplyr::select(contains("a")) %>% pivot_longer(cols=everything())
p <- ggplot(tmp1) + geom_boxplot(aes(y=value,x=name)) + labs(y="AIC",x="")
ggsave(p,filename=paste0(folder_name,"AIC_compare.png"),width=9,height=5)

# finally, plot on a map
title_map <- ""
misscol <- "aquamarine"
legend_text_size <- 0.7
point_size <- 0.5
legend_title_size <- 0.9
limsnll <- c(min(nll1,nll2,nll2b),max(nll1,nll2,nll2b))
limsa <- c(min(a1,a2,a2b),max(a1,a2,a2b))
nrow_facet <- 1
estsf <- cbind(est_all_sf %>% filter(cond_site %in% "London") %>% dplyr::select(c()), tmp %>% add_row(.before=London_index))
p1 <- tm_shape(estsf) + tm_dots(fill="nll1",fill.scale = tm_scale_continuous(limits=limsnll,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="NLL")) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Model 1") 
p2 <- tm_shape(estsf) + tm_dots(fill="nll2",fill.scale = tm_scale_continuous(limits=limsnll,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="NLL")) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Model 2") 
p3 <- tm_shape(estsf) + tm_dots(fill="nll2b",fill.scale = tm_scale_continuous(limits=limsnll,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="NLL")) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Model 2b") 
tmap_save(tmap_arrange(p1,p2,p3,ncol=3),filename=paste0(folder_name,"nll_compare_map.png"),height=6,width=9)

p1 <- tm_shape(estsf) + tm_dots(fill="a1",fill.scale = tm_scale_continuous(limits=limsa,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="AIC")) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Model 1") 
p2 <- tm_shape(estsf) + tm_dots(fill="a2",fill.scale = tm_scale_continuous(limits=limsa,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="AIC")) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Model 2") 
p3 <- tm_shape(estsf) + tm_dots(fill="a2b",fill.scale = tm_scale_continuous(limits=limsa,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="AIC")) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Model 2b") 
tmap_save(tmap_arrange(p1,p2,p3,ncol=3),filename=paste0(folder_name,"AIC_compare_map.png"),height=6,width=9)

# plot also AIC difference
diff12=a1-a2
diff2b2=a2b-a2
diff12b=a1-a2b
tmp2 <- tmp %>% mutate(diff12=a1-a2,diff2b2=a2b-a2,diff12b=a1-a2b)
limsad <- c(min(diff12,diff2b2,diff12b),max(diff12,diff2b2,diff12b))
estsf <- cbind(est_all_sf %>% filter(cond_site %in% "London") %>% dplyr::select(c()), tmp2 %>% add_row(.before=London_index))
p1 <- tm_shape(estsf) + tm_dots(fill="diff12",fill.scale = tm_scale_continuous(limits=limsad,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="AIC")) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Model 1 - Model 2") 
p2 <- tm_shape(estsf) + tm_dots(fill="diff2b2",fill.scale = tm_scale_continuous(limits=limsad,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="AIC")) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Model 2b - Model 2") 
p3 <- tm_shape(estsf) + tm_dots(fill="diff12b",fill.scale = tm_scale_continuous(limits=limsad,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="AIC")) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Model 1 - Model 2b") 
tmap_save(tmap_arrange(p1,p2,p3,ncol=3),filename=paste0(folder_name,"AIC_difference_map.png"),height=6,width=9)

# plot also as boxplots
tmp3 <- tmp2 %>% dplyr::select(contains("diff")) %>% pivot_longer(cols=everything())
p <- ggplot(tmp3 %>% mutate(name=factor(name,levels=c("diff12","diff2b2","diff12b")))) + geom_boxplot(aes(y=value,x=name)) + labs(y="AIC difference",x="") + scale_x_discrete(labels=c("Model 1 - Model 2", "Model 2b - Model 2", "Model 1 - Model 2b"))
ggsave(p,filename=paste0(folder_name,"AIC_compare_diff_boxplot.png"),width=9,height=5)

# 5. Model comparison with mu estimated ---------------------------------------
# need to get model 2 and 2b parameter estimates
to_opt <- function(x1,x2,theta,i,cond_index=London_index,pe_i) {
  # a <- theta[1]
  # b <- theta[2]
  # if (a<0 | a>1 | b<0 | b>1) {return(10^6)}
  y <-  NLL_AGG_onestep(x=data.frame(x1,x2),theta=theta,sigl_hat = pe_i[2], sigu_hat = pe_i[3], deltal_hat = pe_i[4], deltau_hat = pe_i[5])
  return(y)
}
NLL_AGG_wrapper <- function(data_Lap=data_mod_Lap,i,pe_res=pe %>% add_row(.before = London_index),cond_index,v=0.9) {
  pe_i <- as.numeric(unlist(pe_res[i,]))
  data_Lapv <- data_Lap %>% filter(data_Lap[,cond_index]>quantile(data_Lap[,cond_index],v))
  x1 <- as.numeric(unlist(data_Lapv[,London_index]))
  x2 <- as.numeric(unlist(data_Lapv[,i]))
  if (is.na(pe_i[1])) {return(NA)}
  y <- optim(par=c(0.8,0.3,pe_i[1]),fn= to_opt,x1=x1,x2=x2,i=i,cond_index=cond_index,pe_i=pe_i)
  return(y$par) 
}
x <- sapply(1:ncol(data_mod_Lap),FUN=NLL_AGG_wrapper,data_Lap=data_mod_Lap,cond_index=London_index)
tmp <- as.data.frame(do.call(rbind,x))
names(tmp) <- c("a","b","mu")
a_new <- tmp$a
b_new <- tmp$b
mu_new <- tmp$mu
tmp1 <- rbind(data.frame(a=aest,b=best,"method"="original"),data.frame(a=na.omit(tmp$a),b=na.omit(tmp$b),"method"="new")) %>% mutate("iteration"=rep(1:length(aest),2))
# map alpha and beta original and new estimates
plot_ab <- function(tmp) { ggplot(tmp) + 
    geom_line(aes(x=a,y=b,group=iteration),linewidth=0.1) +
    geom_point(aes(x=a,y=b,col=method),alpha=0.7,size=1) +
    xlab(TeX("${\\alpha}$")) +
    ylab(TeX("${\\beta}$")) + 
    scale_color_manual(values = c("original" = "#009ADA", "new" = "#C11432")) + coord_fixed() + theme(axis.text.y = element_text(angle = 90, vjust = 0.5))
}
p <- plot_ab(tmp=tmp1)
ggsave(p,filename=paste0(folder_name,"plot_ab_new_original_mu.png"),width=9,height=5)

# map alpha and beta
title_map <- ""
misscol <- "aquamarine"
legend_text_size <- 0.7
point_size <- 0.5
legend_title_size <- 0.9
limsa <- c(0,1)
limsb <- c(0,0.65)
nrow_facet <- 1
estsf <- cbind(est_all_sf %>% filter(cond_site %in% "London"), data.frame("a_new"=a_new,"b_new" = b_new))
p1 <- tm_shape(estsf) + tm_dots(fill="a",fill.scale = tm_scale_continuous(limits=limsa,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\alpha$"))) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text=TeX("$\\tilde{\\alpha}$")) 
p2 <- tm_shape(estsf) + tm_dots(fill="b",fill.scale = tm_scale_continuous(limits=limsb,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\beta$"))) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text=TeX("$\\tilde{\\beta}$")) 
p3 <- tm_shape(estsf) + tm_dots(fill="a_new",fill.scale = tm_scale_continuous(limits=limsa,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\alpha$"))) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text=TeX("$\\hat{\\alpha}$")) 
p4 <- tm_shape(estsf) + tm_dots(fill="b_new",fill.scale = tm_scale_continuous(limits=limsb,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title=TeX("$\\beta$"))) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text=TeX("$\\hat{\\beta}$")) 
tmap_save(tmap_arrange(p1,p2,p3,p4,ncol=4),filename=paste0(folder_name,"alpha_beta_fixed_res_mu.png"),height=6,width=11)

# model 2b estimates
to_opt <- function(x1,x2,theta,i,cond_index=London_index,pe_i) {
  # a <- theta[1]
  # b <- theta[2]
  # if (a<0 | a>1 | b<0 | b>1) {return(10^6)}
  y <-  NLL_AGG_onestep(x=data.frame(x1,x2),theta=theta,b_hat=0,sigl_hat = pe_i[2], sigu_hat = pe_i[3], deltal_hat = pe_i[4], deltau_hat = pe_i[5])
  return(y)
}
NLL_AGG_wrapper <- function(data_Lap=data_mod_Lap,i,pe_res=pe %>% add_row(.before = London_index),cond_index,v=0.9) {
  pe_i <- as.numeric(unlist(pe_res[i,]))
  data_Lapv <- data_Lap %>% filter(data_Lap[,cond_index]>quantile(data_Lap[,cond_index],v))
  x1 <- as.numeric(unlist(data_Lapv[,London_index]))
  x2 <- as.numeric(unlist(data_Lapv[,i]))
  if (is.na(pe_i[1])) {return(NA)}
  y <- optim(par=c(0.8,pe_i[1]),fn= to_opt,x1=x1,x2=x2,i=i,cond_index=cond_index,pe_i=pe_i)
  return(y$par) 
}
x <- sapply(1:ncol(data_mod_Lap),FUN=NLL_AGG_wrapper,data_Lap=data_mod_Lap,cond_index=London_index)
tmp <- as.data.frame(do.call(rbind,x))
names(tmp) <- c("a","mu")
a_new2 <- na.omit(tmp$a)
b_new2 <- 0
mu_new2 <- na.omit(tmp$mu)
#to_opt(x1=x1,x2=x2,i=1,theta=c(0.8,0.3),pe_i=pe_i)
tmp1 <- rbind(data.frame(a=a_new2,b=b_new2,"method"="model_2b"),data.frame(a=na.omit(a_new),b=na.omit(b_new),"method"="model_2")) %>% mutate("iteration"=rep(1:length(aest),2))
# map alpha and beta original and new estimates
plot_ab <- function(tmp) { ggplot(tmp) + 
    geom_line(aes(x=a,y=b,group=iteration),linewidth=0.1) +
    geom_point(aes(x=a,y=b,col=method),alpha=0.7,size=1) +
    xlab(TeX("${\\alpha}$")) +
    ylab(TeX("${\\beta}$")) + 
    scale_color_manual(values = c("model_2b" = "#009ADA", "model_2" = "#C11432")) + coord_fixed() + theme(axis.text.y = element_text(angle = 90, vjust = 0.5))
}
p <- plot_ab(tmp=tmp1)
ggsave(p,filename=paste0(folder_name,"plot_ab_new_original_model2b_mu.png"),width=9,height=5)

# calculate log-likelihood

NLL_AGG_wrapper <- function(data_Lap=data_mod_Lap,i,pe_res=pe,cond_index,v=0.9,a=na.omit(a_new),b=na.omit(b_new),mu=na.omit(mu_new)) {
  pe_i <- as.numeric(unlist(pe_res[i,]))
  data_Lapv <- data_Lap %>% filter(data_Lap[,cond_index]>quantile(data_Lap[,cond_index],v))
  x1 <- as.numeric(unlist(data_Lapv[,London_index]))
  ires <- c(1:length(data_Lap))[-London_index][i]
  x2 <- as.numeric(unlist(data_Lapv[,ires]))
  y <- NLL_AGG_onestep(x=data.frame(x1,x2),theta=c(),a_hat=a[i],b_hat=b[i],mu_hat=mu[i],sigl_hat = pe_i[2], sigu_hat = pe_i[3], deltal_hat = pe_i[4], deltau_hat = pe_i[5])
  return(y) 
}
nll2 <- sapply(1:(ncol(data_mod_Lap)-1),FUN=NLL_AGG_wrapper,data_Lap=data_mod_Lap,cond_index=London_index)

NLL_AGG_wrapper <- function(data_Lap=data_mod_Lap,i,pe_res=pe,cond_index,v=0.9,a=a_new2,mu=na.omit(mu_new2)) {
  pe_i <- as.numeric(unlist(pe_res[i,]))
  data_Lapv <- data_Lap %>% filter(data_Lap[,cond_index]>quantile(data_Lap[,cond_index],v))
  x1 <- as.numeric(unlist(data_Lapv[,London_index]))
  ires <- c(1:length(data_Lap))[-London_index][i]
  x2 <- as.numeric(unlist(data_Lapv[,ires]))
  y <- NLL_AGG_onestep(x=data.frame(x1,x2),theta=c(),a_hat=a[i],b_hat=0,mu_hat=mu[i],sigl_hat = pe_i[2], sigu_hat = pe_i[3], deltal_hat = pe_i[4], deltau_hat = pe_i[5])
  return(y) 
}
nll2b <- sapply(1:(ncol(data_mod_Lap)-1),FUN=NLL_AGG_wrapper,data_Lap=data_mod_Lap,cond_index=London_index)

# calculate AIC
a1 <- 2*nll1+2*4
a2 <- 2*nll2+2*3
a2b <- 2*nll2b+2*2

tmp <- data.frame(nll1,nll2,nll2b,a1,a2,a2b)

# boxplots of negative log-likelihood
tmp1 <- tmp %>% dplyr::select(contains("nll")) %>% pivot_longer(cols=everything())
p <- ggplot(tmp1) + geom_boxplot(aes(y=value,x=name)) + labs(y="Negative log-likelihood",x="")
ggsave(p,filename=paste0(folder_name,"nll_compare_mu.png"),width=9,height=5)
# boxplots of AIC
tmp1 <- tmp %>% dplyr::select(contains("a")) %>% pivot_longer(cols=everything())
p <- ggplot(tmp1) + geom_boxplot(aes(y=value,x=name)) + labs(y="AIC",x="")
ggsave(p,filename=paste0(folder_name,"AIC_compare_mu.png"),width=9,height=5)

# finally, plot on a map
title_map <- ""
misscol <- "aquamarine"
legend_text_size <- 0.7
point_size <- 0.5
legend_title_size <- 0.9
limsnll <- c(min(nll1,nll2,nll2b),max(nll1,nll2,nll2b))
limsa <- c(min(a1,a2,a2b),max(a1,a2,a2b))
nrow_facet <- 1
estsf <- cbind(est_all_sf %>% filter(cond_site %in% "London") %>% dplyr::select(c()), tmp %>% add_row(.before=London_index))
p1 <- tm_shape(estsf) + tm_dots(fill="nll1",fill.scale = tm_scale_continuous(limits=limsnll,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="NLL")) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Model 1") 
p2 <- tm_shape(estsf) + tm_dots(fill="nll2",fill.scale = tm_scale_continuous(limits=limsnll,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="NLL")) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Model 2") 
p3 <- tm_shape(estsf) + tm_dots(fill="nll2b",fill.scale = tm_scale_continuous(limits=limsnll,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="NLL")) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Model 2b") 
tmap_save(tmap_arrange(p1,p2,p3,ncol=3),filename=paste0(folder_name,"nll_compare_map_mu.png"),height=6,width=9)

p1 <- tm_shape(estsf) + tm_dots(fill="a1",fill.scale = tm_scale_continuous(limits=limsa,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="AIC")) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Model 1") 
p2 <- tm_shape(estsf) + tm_dots(fill="a2",fill.scale = tm_scale_continuous(limits=limsa,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="AIC")) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Model 2") 
p3 <- tm_shape(estsf) + tm_dots(fill="a2b",fill.scale = tm_scale_continuous(limits=limsa,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="AIC")) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Model 2b") 
tmap_save(tmap_arrange(p1,p2,p3,ncol=3),filename=paste0(folder_name,"AIC_compare_map_mu.png"),height=6,width=9)

# plot also AIC difference
diff12=a1-a2
diff2b2=a2b-a2
diff12b=a1-a2b
tmp2 <- tmp %>% mutate(diff12=a1-a2,diff2b2=a2b-a2,diff12b=a1-a2b)
limsad <- c(min(diff12,diff2b2,diff12b),max(diff12,diff2b2,diff12b))
estsf <- cbind(est_all_sf %>% filter(cond_site %in% "London") %>% dplyr::select(c()), tmp2 %>% add_row(.before=London_index))
p1 <- tm_shape(estsf) + tm_dots(fill="diff12",fill.scale = tm_scale_continuous(limits=limsad,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="AIC")) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Model 1 - Model 2") 
p2 <- tm_shape(estsf) + tm_dots(fill="diff2b2",fill.scale = tm_scale_continuous(limits=limsad,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="AIC")) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Model 2b - Model 2") 
p3 <- tm_shape(estsf) + tm_dots(fill="diff12b",fill.scale = tm_scale_continuous(limits=limsad,values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="AIC")) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Model 1 - Model 2b") 
tmap_save(tmap_arrange(p1,p2,p3,ncol=3),filename=paste0(folder_name,"AIC_difference_map_mu.png"),height=6,width=9)

# plot middle above 20
tmp4 <- tmp2 %>% mutate("is_big"=(diff2b2>20))
estsf <- cbind(est_all_sf %>% filter(cond_site %in% "London") %>% dplyr::select(c()), tmp4 %>% add_row(.before=London_index))
p1 <- tm_shape(estsf %>% dplyr::filter(is_big %in% c(TRUE,NA))) + tm_dots(fill="diff2b2",fill.scale = tm_scale_continuous(values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="AIC")) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Model 2b - Model 2") 
p2 <- tm_shape(estsf) + tm_dots(fill="is_big",fill.scale = tm_scale_categorical(values=c("black","#C11432"),value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="mean(y_sim)")) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Is bigger than 20") 
tmap_save(tmap_arrange(p1,p2,ncol=2),filename=paste0(folder_name,"AIC_difference_2b_2.png"),height=6,width=6)

tmp4 <- tmp2 %>% mutate("is_big"=(diff2b2>1))
estsf <- cbind(est_all_sf %>% filter(cond_site %in% "London") %>% dplyr::select(c()), tmp4 %>% add_row(.before=London_index))
p1 <- tm_shape(estsf %>% dplyr::filter(is_big %in% c(TRUE,NA))) + tm_dots(fill="diff2b2",fill.scale = tm_scale_continuous(values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="AIC")) +  tm_layout(legend.position=c("right","top"),legend.height = 12,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Model 2b - Model 2") 
p2 <- tm_shape(estsf) + tm_dots(fill="is_big",fill.scale = tm_scale_categorical(values=c("black","#C11432"),value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="mean(y_sim)")) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Is bigger than 1") 
tmap_save(tmap_arrange(p1,p2,ncol=2),filename=paste0(folder_name,"AIC_difference_2b_2_1.png"),height=6,width=6)

# plot also as boxplots
tmp3 <- tmp2 %>% dplyr::select(contains("diff")) %>% pivot_longer(cols=everything())
p <- ggplot(tmp3 %>% mutate(name=factor(name,levels=c("diff12","diff2b2","diff12b")))) + geom_boxplot(aes(y=value,x=name)) + labs(y="AIC difference",x="") + scale_x_discrete(labels=c("Model 1 - Model 2", "Model 2b - Model 2", "Model 1 - Model 2b"))
ggsave(p,filename=paste0(folder_name,"AIC_compare_diff_boxplot_mu.png"),width=9,height=5)

# 2. recreate simulated fields -----------------------------------------------
# get index for London
x <- july3_obs[London_index]
# transform to Laplace
xL <- unif_laplace_pit(x)
# get residual margin parameters
pe <- as.data.frame(result_new[[12]] %>% dplyr::select(mu_agg,sigl,sigu,deltal,deltau))
names(pe) <- c("mu","sigl","sigu","deltal","deltau")
# get fields
# reconstruct the fields
y_sim <- apply(random10N,MARGIN=c(2),FUN=function(xk){xL*aest+xL^best*xk})
# add row for the conditioning site
y_sim <- as.data.frame(y_sim)
#y_sim <- y_sim %>% add_row(.before=London_index)
y_sim[London_index,] <- xL
names(y_sim) <- paste0("random",1:10)
# plot on standard Normal scale
tmpsf <- st_as_sf(cbind(y_sim,xyUK20_sf)) %>% pivot_longer(cols=contains("random"))
t <- tm_shape(tmpsf %>% mutate("name"=factor(name, levels=unique(tmpsf$name)))) + tm_dots(fill="value",size=0.5,fill.scale = tm_scale_continuous(values="-brewer.rd_bu",limits=c(-10,20)),fill.legend = tm_legend(title = "")) + tm_facets(by="name",ncol=5)
tmap_save(tm=t, filename=paste0(folder_name,"random10_y_sim_London.png"),width=10,height=8)

# transform all 10 onto 2 different scales
#unif_orig_P2q(u=july3_obs,P2q=P2q_sites2,gpdpar = gpdpar_sites2)
x_sim1 <- apply(y_sim,MARGIN=c(2),FUN=function(xk){unif_orig_P2q(u=plaplace(xk),P2q=P2q_sites2,gpdpar = gpdpar_sites2)}) %>% as.data.frame()
names(x_sim1) <- paste0("random",1:10)
# plot on standard Normal scale
tmpsf <- st_as_sf(cbind(x_sim1,xyUK20_sf)) %>% pivot_longer(cols=contains("random"))
t <- tm_shape(tmpsf %>% mutate("name"=factor(name, levels=unique(tmpsf$name)))) + tm_dots(fill="value",size=0.5,fill.scale = tm_scale_intervals(values="-brewer.rd_bu",breaks=seq(10,50,5)),fill.legend = tm_legend(title = "")) + tm_facets(by="name",ncol=5)
tmap_save(tm=t, filename=paste0(folder_name,"random10_x_sim2026.png"),width=10,height=8)
t <- tm_shape(tmpsf %>% mutate("name"=factor(name, levels=unique(tmpsf$name)))) + tm_dots(fill="value",size=0.5,fill.scale = tm_scale_continuous(values="-brewer.rd_bu",limits=c(10,50)),fill.legend = tm_legend(title = "")) + tm_facets(by="name",ncol=5)
tmap_save(tm=t, filename=paste0(folder_name,"random10_xcont_sim2026_London.png"),width=10,height=8)

x_sim2 <- apply(y_sim,MARGIN=c(2),FUN=function(xk){unif_orig_P2q(u=plaplace(xk),P2q=P2q_sites3,gpdpar = gpdpar_sites3)}) %>% as.data.frame()
names(x_sim2) <- paste0("random",1:10)
# plot on standard Normal scale
tmpsf <- st_as_sf(cbind(x_sim2,xyUK20_sf)) %>% pivot_longer(cols=contains("random"))
t <- tm_shape(tmpsf %>% mutate("name"=factor(name, levels=unique(tmpsf$name)))) + tm_dots(fill="value",size=0.5,fill.scale = tm_scale_intervals(values="-brewer.rd_bu",breaks=seq(10,50,5)),fill.legend = tm_legend(title = "")) + tm_facets(by="name",ncol=5)
tmap_save(tm=t, filename=paste0(folder_name,"random10_x_sim2076.png"),width=10,height=8)
t <- tm_shape(tmpsf %>% mutate("name"=factor(name, levels=unique(tmpsf$name)))) + tm_dots(fill="value",size=0.5,fill.scale = tm_scale_continuous(values="-brewer.rd_bu",limits=c(10,50)),fill.legend = tm_legend(title = "")) + tm_facets(by="name",ncol=5)
tmap_save(tm=t, filename=paste0(folder_name,"random10_xcont_sim2076_London.png"),width=10,height=8)

# check the mean value at each site across many fields
random10000 <- as.data.frame(  spam::rmvnorm(n=100000,Sigma = Zcov)  )
random10000N <- sapply(1:ncol(random10000),FUN=function(k) {Normal_AGG_PIT(z = random10000[,k],theta=c(pe$mu[k],pe$sigl[k],pe$sigu[k],pe$deltal[1],pe$deltau[1]))}) %>% as.data.frame()
names(random10000N) <- names(Z)
# calculate mean at each site
# reconstruct the fields
aest <- discard(est_all_sf %>% filter(cond_site=="London") %>% pull(a),is.na)
best <- discard(est_all_sf %>% filter(cond_site=="London") %>% pull(b),is.na)
y_sim <- apply(random10000N,MARGIN=c(1),FUN=function(xk){xL*aest+xL^best*xk})
# add row for the conditioning site
y_sim <- as.data.frame(y_sim)
tmp <- apply(y_sim %>% add_row(.before = London_index),MARGIN = c(1),FUN=mean)
tmpsf <- st_as_sf(cbind(tmp,xyUK20_sf))
lims <- c(1.5,6.1)
p1 <- tm_shape(tmpsf) + tm_dots(fill="tmp",fill.scale = tm_scale_continuous(values="viridis",limits=lims,value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="mean(y_sim)")) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Mean value at each site") 
# is greater than xL
tmpsf <- tmpsf %>% mutate("is_big"=(tmp>xL))
p2 <- tm_shape(tmpsf) + tm_dots(fill="is_big",fill.scale = tm_scale_categorical(values=c("black","#C11432"),value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="mean(y_sim)")) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Is bigger than cond. site") 
tmap_save(tmap_arrange(p1,p2,ncol=2),filename=paste0(folder_name,"y_sim_mean.png"),height=6,width=6)

# plot Y2,Y1 on Laplace margins
tmp <- y_sim %>% add_row(.before=London_index)

p <- data.frame("x"=as.numeric(unlist(y_sim[find_site_index(Hull,grid_uk = xyUK20_sf),]))) %>%  ggplot() + geom_density(aes(x=x))
ggsave(p, filename = paste0(folder_name,"Hull_simulated.png"))
p <- data.frame("x"=as.numeric(unlist(y_sim[find_site_index(Hull,grid_uk = xyUK20_sf),]))) %>%  ggplot() + geom_histogram(aes(x=x))
ggsave(p, filename = paste0(folder_name,"Hull_simulated_hist.png"))

# try plotting both simulated and observed values
tmpo <- data_mod_Lap %>% dplyr::select(all_of(c(find_site_index(Hull,grid=xyUK20_sf),London_index)))
tmpo_ex <- tmpo %>% dplyr::filter(Y100>quantile(Y100,0.9,na.rm=TRUE))
tmps <- data.frame("y"=as.numeric(unlist(y_sim[find_site_index(Hull,grid_uk = xyUK20_sf),])),"x"=xL)
names(tmps) <- names(tmpo)
tmpa <- rbind(tmpo %>% mutate("sim_obs"=c("observed")),tmps %>% mutate("sim_obs"=c("simulated")))

ggplot(tmpa) + geom_point(aes(x=Y100,y=Y321,col=sim_obs),size=0.2) + 
  scale_color_manual(values = c("observed" = "black", "simulated" = "#C11432")) 

# 3. repeat with expected value of the residuals ------------------------------
Ez <- apply(pe,MARGIN = c(1),FUN=AGG_mean)
# map mean on a map compared with observed residuals mean
Ez_df <- data.frame(Ez) %>% add_row(.before=London_index)
tmp <- data.frame(Zmean=apply(Z,MARGIN = c(2),FUN=mean)) %>% add_row(.before = London_index)
tmpsf <- st_as_sf(cbind(Ez_df,tmp,xyUK20_sf))
lims <- c(min(Ez, apply(Z,MARGIN = c(2),FUN=mean)),max(Ez,apply(Z,MARGIN = c(2),FUN=mean)))
p1 <- tm_shape(tmpsf) + tm_dots(fill="Zmean",fill.scale = tm_scale_continuous(values="viridis",limits=lims,value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="Value")) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Observed residual mean") 
p2 <- tm_shape(tmpsf) + tm_dots(fill="Ez",fill.scale = tm_scale_continuous(values="viridis",limits=lims,value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="Value")) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Expected value of fitted AGG") 
tmap_save(tmap_arrange(p1,p2,ncol=2),filename=paste0(folder_name,"Z_mean_obs_fitted.png"),height=6,width=6)

# reconstruct the fields
aest <- discard(est_all_sf %>% filter(cond_site=="London") %>% pull(a),is.na)
best <- discard(est_all_sf %>% filter(cond_site=="London") %>% pull(b),is.na)
y_sim <- sapply(Ez,FUN=function(xk){xL*aest+xL^best*xk})
# add row for the conditioning site
y_sim <- as.data.frame(y_sim)
tmp <- apply(y_sim %>% add_row(.before = London_index),MARGIN = c(1),FUN=mean)
tmpsf <- st_as_sf(cbind(tmp,xyUK20_sf))
lims <- c(1.5,6.1) # set common limits for comparison
p1 <- tm_shape(tmpsf) + tm_dots(fill="tmp",fill.scale = tm_scale_continuous(values="viridis",value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="mean(y_sim)")) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Mean value at each site") 
# is greater than xL
tmpsf <- tmpsf %>% mutate("is_big"=(tmp>xL))
p2 <- tm_shape(tmpsf) + tm_dots(fill="is_big",fill.scale = tm_scale_categorical(values=c("black","#C11432"),value.na=misscol,label.na = "Conditioning\n site"),size=point_size, fill.legend = tm_legend(title="mean(y_sim)")) +  tm_layout(legend.position=c("right","top"),legend.height = 10,legend.text.size = legend_text_size,legend.title.size=legend_title_size,legend.reverse=TRUE,frame=FALSE) + tm_title(text="Is bigger than cond. site") 
tmap_save(tmap_arrange(p1,p2,ncol=2),filename=paste0(folder_name,"y_sim_mean_exact.png"),height=6,width=6)





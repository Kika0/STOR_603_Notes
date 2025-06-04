# iterative estimation for AGG parameters -------

library(tmap)
library(tidyverse)
library(latex2exp)
theme_set(theme_bw())
theme_replace(
  panel.spacing = unit(2, "lines"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA) )
file.sources = list.files(pattern="*helpers.R")
sapply(file.sources,source,.GlobalEnv)
# source("spatial_parameter_estimation.R") # for spatial_par_est function
load("data_processed/temperature_data.RData")
load("data_processed/spatial_helper.RData")

res_margin_par_est_iterative <- function(obs_res,method="AGG",Nite=10,deltal_init=2,deltau_init=2) {
  d <- ncol(obs_res)
    if (method=="AGG") {
      mu <- sigl <- sigu <- deltal <- deltau <- data.frame(matrix(ncol=(Nite+1),nrow = (d)))
      if (is.numeric(deltal_init) & is.numeric(deltau_init)) {
        deltal[,1] <- deltal_init
        deltau[,1] <- deltau_init
      } else {
        deltal[,1] <- 2
        deltau[,1] <- 2
      }
      for (j in 1:ncol(obs_res)) {
        Z2 <- as.numeric(unlist(obs_res[,j]))
      # estimate initial values
      opt <- optim(fn=NLL_AGG,x=Z2,deltal_hat=as.numeric(deltal[1,1]),deltau_hat=as.numeric(deltau[1,1]),par=c(mean(Z2),sd(Z2),sd(Z2)),control=list(maxit=2000),method = "Nelder-Mead")
      mu[j,1] <- opt$par[1]
      sigl[j,1] <- opt$par[2]
      sigu[j,1] <- opt$par[3]
      for (i in 1:Nite) {
    #    init_parAGG <- c(mu[j,i],sigl[j,i],sigu[j,i])
        # update deltas
        opt <- optim(fn=NLL_AGG,x=Z2,mu_hat=as.numeric(mu[j,i]),sigl_hat=as.numeric(sigl[j,i]),sigu_hat=as.numeric(sigu[j,i]),par=c(as.numeric(deltal[j,i]),as.numeric(deltau[j,i])),control=list(maxit=2000),method = "Nelder-Mead")
        deltal[j,i+1] <- opt$par[length(opt$par)-1]
        deltau[j,i+1] <- opt$par[length(opt$par)]
        # reestimate mu,sigl and sigu
        opt <- optim(fn=NLL_AGG,x=Z2,deltal_hat=as.numeric(deltal[j,i+1]),deltau_hat=as.numeric(deltau[j,i+1]),par=c(as.numeric(mu[j,i]),as.numeric(sigl[j,i]),as.numeric(sigu[j,i])),control=list(maxit=2000),method = "Nelder-Mead")
        mu[j,i+1] <- opt$par[1]
        sigl[j,i+1] <- opt$par[2]
        sigu[j,i+1] <- opt$par[3]
      }
    }
    
  }
  return(list(mu,sigl,sigu,deltal,deltau))
}


# load all three parameter estimates sf objects ----
load("data_processed/N9000_sequential2_AGG_all12sites.RData")
est_all <- est_all_sf

# check current function for estimating parameters of residual margins ----
# focus on Birmingham conditioning site
a <- na.omit(as.data.frame(est_all %>% dplyr::filter(cond_site=="Birmingham"))$a)
b <- na.omit(as.data.frame(est_all_sf %>% dplyr::filter(cond_site=="Birmingham"))$b)
res <- observed_residuals(df=data_mod_Lap,a=a,b=b,v = 0.9)
#old_par <- res_margin_par_est(obs_res=res)

# try new function
Nite <- 50
d <- ncol(res)
time.start <- Sys.time()
ite_par <- res_margin_par_est_iterative(obs_res=res,Nite=Nite)
time.end <- Sys.time()
time.start-time.end
# print for a selected residual site
mu_randomsite <- data.frame(ite=1:(Nite+1),mu=as.numeric(ite_par[[1]][1,]))
# boxplots of how estimates change through iterations
tmpmu <- ite_par[[1]]%>% pivot_longer(cols=everything(),names_to="iteration",values_to="count") %>% mutate("iteration"=factor(rep(1:(Nite+1),d),levels=1:(Nite+1)))
tmpsigl <- ite_par[[2]]%>% pivot_longer(cols=everything(),names_to="iteration",values_to="count") %>% mutate("iteration"=factor(rep(1:(Nite+1),d),levels=1:(Nite+1)))
tmpsigu <- ite_par[[3]]%>% pivot_longer(cols=everything(),names_to="iteration",values_to="count") %>% mutate("iteration"=factor(rep(1:(Nite+1),d),levels=1:(Nite+1)))
tmpdeltal <- ite_par[[4]]%>% pivot_longer(cols=everything(),names_to="iteration",values_to="count") %>% mutate("iteration"=factor(rep(1:(Nite+1),d),levels=1:(Nite+1)))
tmpdeltau <- ite_par[[5]]%>% pivot_longer(cols=everything(),names_to="iteration",values_to="count") %>% mutate("iteration"=factor(rep(1:(Nite+1),d),levels=1:(Nite+1)))


ggplot(tmpmu) + geom_boxplot(aes(x=iteration,y=count))                     

# compare a boxplot of the previous approach, first iteration and last iteration
mu_agg <- na.omit(as.data.frame(est_all %>% dplyr::filter(cond_site=="Birmingham"))$mu_agg)
sigl <- na.omit(as.data.frame(est_all_sf %>% dplyr::filter(cond_site=="Birmingham"))$sigl)
sigu <- na.omit(as.data.frame(est_all_sf %>% dplyr::filter(cond_site=="Birmingham"))$sigu)
deltal <- na.omit(as.data.frame(est_all_sf %>% dplyr::filter(cond_site=="Birmingham"))$deltal)
deltau <- na.omit(as.data.frame(est_all_sf %>% dplyr::filter(cond_site=="Birmingham"))$deltau)

# join with the first and the last iteration
iterative_AGG_boxplot <- function(tmp=tmpmu,par_AGG=mu_agg,y_label=TeX("$\\mu_{AGG}$"),N=Nite) {
tmpdf <- rbind(data.frame(count=par_AGG,iteration="AGG"),tmp %>% filter(iteration %in% c(1,N+1)))
p <- ggplot(tmpdf) + geom_boxplot(aes(x=iteration,y=count,col=factor(iteration))) +  scale_color_manual(name="Method",labels=c("First iteration","Last iteration","AGG"),                                                                  
                                                                            values=c("black","black","#C11432")) +ylab(y_label)
ggsave(p,filename=paste0("../Documents/iterativeAGGboxplot",gsub("\\[|]","",as.character(y_label)),".png"),width=5,height=5)
}

iterative_AGG_boxplot()
iterative_AGG_boxplot(tmp=tmpsigl,par_AGG = sigl,y_label=TeX("$\\sigma_{l}$"))
iterative_AGG_boxplot(tmp=tmpsigu,par_AGG = sigu,y_label=TeX("$\\sigma_{u}$"))
iterative_AGG_boxplot(tmp=tmpdeltal,par_AGG = deltal,y_label=TeX("$\\delta_{l}$"))
iterative_AGG_boxplot(tmp=tmpdeltau,par_AGG = deltau,y_label=TeX("$\\delta_{u}$"))

# join back to spatial

iterative_AGG_map <- function(tmp=tmpmu,par_AGG=mu_agg,y_label=TeX("$\\mu_{AGG}$"),N=Nite) {
  tmpdf <- rbind(data.frame(count=par_AGG,iteration="AGG"),tmp %>% filter(iteration %in% c(1,N+1)) %>% arrange(iteration))
estBirm <- est_all %>% dplyr::filter(cond_site=="Birmingham")
 tmpsf <-  rbind(estBirm %>% mutate(iteration="AGG"),estBirm %>% mutate(iteration=as.character(1)),estBirm %>% mutate(iteration=as.character(N+1))) %>% mutate(value=res)
tmpsf$value[!is.na( tmpsf$res)] <- tmpdf$count

tm <- tm_shape(tmpsf)  + tm_dots(col="value",palette="viridis",n=10,size=0.3,colorNA="aquamarine",title=y_label,textNA = "Conditioning site") + tm_facets(by= c("iteration"))+  tm_layout(legend.outside.size=0.3,asp=0.5,legend.text.size = 1,legend.title.size=1.5)

    tmap_save(tm,filename=paste0("../Documents/iterativeAGGmap",gsub("\\[|]","",as.character(y_label)),".png"),width=10,height=6)
}

iterative_AGG_map()
iterative_AGG_map(tmp=tmpsigl,par_AGG = sigl,y_label=TeX("$\\sigma_{l}$"))
iterative_AGG_map(tmp=tmpsigu,par_AGG = sigu,y_label=TeX("$\\sigma_{u}$"))
iterative_AGG_map(tmp=tmpdeltal,par_AGG = deltal,y_label=TeX("$\\delta_{l}$"))
iterative_AGG_map(tmp=tmpdeltau,par_AGG = deltau,y_label=TeX("$\\delta_{u}$"))


#library(tmap) # spatial map plots
library(gridExtra)
library(sf) # for handling spatial sf objects
library(tidyverse)
library(latex2exp) # latex expressions for plot labels
library(texmex) # for pollution data
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
df_sites <- data.frame(Birmingham,Glasgow,London,Inverness,Lancaster,Newcastle,Cromer,Hull,Lowestoft,Truro,Dolgellau,Bournemouth)

# 1. plot of London/Birmingham and London/Inverness on original margins -------
Birm_index <- find_site_index(site=Birmingham,grid_uk = xyUK20_sf)
London_index <- find_site_index(site=London,grid_uk = xyUK20_sf)
Inverness_index <- find_site_index(site=Inverness,grid_uk = xyUK20_sf)
tmp <- data_mod[,c(Birm_index,Glasgow_index,Inverness_index)]
names(tmp) <- names(df_sites)[c(1,3,4)]
p1 <- ggplot(tmp) + geom_point(aes(x=London,y=Birmingham),size=0.5)
p2 <- ggplot(tmp) + geom_point(aes(x=London,y=Inverness),size=0.5)
p <- grid.arrange(p1,p2,ncol=2)
ggsave(p,filename= "../Documents/London_Birmingham_Inverness_original.png",width=8,height=4)

# 2. plot illustrating conditioning on a variable ------
tmp <- data_mod_Lap[,c(Birm_index,Glasgow_index,Inverness_index)]
names(tmp) <- names(df_sites)[c(1,3,4)]
q <- 0.9
vL <- qlaplace(q)
tmp$tf <- (tmp$London>vL)
vL <- quantile(tmp[,1],q)
p1 <- ggplot(tmp) + geom_point(aes(x=London,y=Birmingham),size=0.5) + xlab(TeX("$Y_1$")) + ylab(TeX("$Y_2$")) + scale_color_manual(values = c("black", "#009ADA")) + theme(legend.position="none") + coord_fixed() + xlim(-10,10) + ylim(-10,10) 
p2 <- ggplot(tmp) + geom_point(aes(x=London,y=Inverness),size=0.5) + xlab(TeX("$Y_1$")) + ylab(TeX("$Y_3$")) + scale_color_manual(values = c("black", "#009ADA")) + theme(legend.position="none") + coord_fixed() + xlim(-10,10) + ylim(-10,10)  
p1n <- ggplot(tmp) + geom_point(aes(x=London,y=Birmingham,col=tf),size=0.5) + scale_color_manual(values = c("black", "#009ADA")) + theme(legend.position="none") + coord_fixed() + xlim(-10,10) + ylim(-10,10) +geom_vline(xintercept=vL,color="#009ADA",linetype="dashed")
p2n <- ggplot(tmp) + geom_point(aes(x=London,y=Inverness,col=tf),size=0.5)  + scale_color_manual(values = c("black", "#009ADA")) + theme(legend.position="none") + coord_fixed() + xlim(-10,10) + ylim(-10,10) +geom_vline(xintercept=vL,color="#009ADA",linetype="dashed") 
p <- grid.arrange(p1,p2,ncol=2)
ggsave(p,filename="../Documents/condmodel_illustration.png",width=8,height=4)
pn <- grid.arrange(p1n,p2n,ncol=2)
ggsave(pn,filename="../Documents/condmodel_illustrationvL.png",width=8,height = 4)

# 3. illustrate pollution data for O3 conditioning on PM10 --------
winter_lap <- as.data.frame((winter %>% apply(c(2),FUN=row_number))/(nrow(winter)+1)) %>% apply(c(1,2),FUN=unif_laplace_pit) %>% as.data.frame()
vL <- quantile(winter$PM10,0.7)
tmp <- winter %>% mutate(above_thres= as.character(winter$PM10>vL))
vL1 <- quantile(summer$PM10,0.7)
tmp1 <- summer %>% mutate(above_thres= as.character(summer$PM10>vL1))
p <- grid.arrange(ggplot(tmp) + geom_point(aes(x=PM10,y=O3,col=above_thres),size=0.5,alpha=0.5) +
                    scale_color_manual(values = c("FALSE"="black","TRUE" = "#009ADA")) +
                    geom_vline(xintercept=vL,color="#009ADA",linetype="dashed") +
                    xlab(TeX("$PM_{10}$")) + ylab(TeX("$O_3$"))+
                    theme(legend.position = "none") ,
                  ggplot(tmp1) + geom_point(aes(x=PM10,y=O3,col=above_thres),size=0.5,alpha=0.5) +
                    scale_color_manual(values = c("FALSE"="black","TRUE" = "#C11432")) +
                    geom_vline(xintercept=vL1,color="#C11432",linetype="dashed") +
                    xlab(TeX("$PM_{10}$")) + ylab(TeX("$O_3$")) +
                    theme(legend.position="none") ,ncol=2)
ggsave(p,filename = "../Documents/PM10_O3_illustrate.png",height=3,width=6.5)

# 4. density plot illustration for pollution data residual margins -------
v <- 0.7
j <- 4 # change j to repeat analysis cond. on other variables

# transform to Laplace margins
summer_lap <- as.data.frame((summer %>% apply(c(2),FUN=row_number))/(nrow(summer)+1)) %>% apply(c(1,2),FUN=unif_laplace_pit) %>% as.data.frame()
pes <-  par_est(df=summer_lap,v=v,given=j,margin = "Normal", method = "sequential2")
obsress <- observed_residuals(df = summer_lap,given = j,v = v,a = pes$a,b=pes$b) 
Z2 <- as.numeric(unlist(obsress[,1]))

# calculate parameters of each method
L1 <- res_margin_par_est(obs_res = obsress,method="Normal")
L2 <- res_margin_par_est(obs_res = obsress,method="GenGaus")
L3 <- res_margin_par_est(obs_res = obsress,method="AGGdelta")
L4 <- res_margin_par_est(obs_res = obsress,method="AGGsig")
L5 <- res_margin_par_est(obs_res = obsress,method="AGG")

d <- ncol(winter)
L1 <- L1 %>% mutate(method="Normal",AIC=2*likres+2*2) %>% mutate(res=c(1:d)[-j]) %>% mutate(mu_agg=mu,sigl=sig,sigu=sig,deltal=2,deltau=2) %>% dplyr::select(c(likres,mu_agg,sigl,sigu,deltal,deltau,AIC,method,res))
L2 <- L2 %>% mutate(method="GenGaus",AIC=2*likres+2*3) %>% mutate(res=c(1:d)[-j]) %>% mutate(sigl=sig_agg,sigu=sig_agg,deltal=delta,deltau=delta) %>% dplyr::select(c(likres,mu_agg,sigl,sigu,deltal,deltau,AIC,method,res))
L3 <- L3 %>% mutate(method="AGGdelta",AIC=2*likres+2*4) %>% mutate(res=c(1:d)[-j]) %>% mutate(sigl=sig_agg,sigu=sig_agg) %>% dplyr::select(c(likres,mu_agg,sigl,sigu,deltal,deltau,AIC,method,res))
L4 <- L4 %>% mutate(method="AGGsig",AIC=2*likres+2*4) %>% mutate(res=c(1:d)[-j]) %>% mutate(deltal=delta,deltau=delta) %>% dplyr::select(c(likres,mu_agg,sigl,sigu,deltal,deltau,AIC,method,res))
L5 <- L5 %>% mutate(method="AGG",AIC=2*likres+2*5) %>% mutate(res=c(1:d)[-j]) %>% dplyr::select(c(likres,mu_agg,sigl,sigu,deltal,deltau,AIC,method,res))

tmp <- rbind(L1,L2,L3,L4,L5) %>% mutate(given=j) %>%  mutate(method=factor(method),given=factor(as.character(given)),res=factor(as.character(res)))

# explore fits in the first iteration
if (j==1) {pickres <- 2} else {pickres <- 1}
tmp1giv2 <- tmp %>% filter(given==as.character(j),res==as.character(pickres))
methods <- unique(tmp$method)
tr1 <- data.frame(y=as.numeric(),x=as.numeric(),Method=as.character())
nllv <- c()
method_par <- data.frame("Normal" = c(as.numeric(tmp1giv2[1,2:6])),"GenGaus" = c(as.numeric(tmp1giv2[2,2:6])),"AGGdelta" = c(as.numeric(tmp1giv2[3,2:6])),"AGGsig" = c(as.numeric(tmp1giv2[4,2:6])),"AGG" = c(as.numeric(tmp1giv2[5,2:6])))

Z <- as.numeric(unlist(obsress[,1]))
x <- seq(min(Z)-1,max(Z)+1,by=0.1)
for (i in 1:5) {
  methodpar <- as.numeric(method_par[,i])
  tr1 <- rbind(tr1,data.frame(y=AGG_density(x=x,theta = methodpar),x=x,Method=names(method_par[i])))
}

pl <- ggplot(data.frame(x=Z)) + geom_density(aes(x=x),linetype="dashed") + xlim(c(min(Z)-1,max(Z)+1))
pl1 <- pl + geom_line(data=tr1 %>% mutate(Method_optim=factor(Method)),aes(x=x,y=y,col=Method)) +
  xlab("Observed residuals kernel smoothed density") + ylab("Density") + ggtitle("Kernel smoothed residual density and AGG fits")+ scale_color_manual(values=c("Normal"="#070707","GenGaus"="#C11432","AGGdelta"="#009ADA","AGGsig"="#66A64F","AGG"="#FDD10A")) + ylim(c(0,0.5))
ggsave(filename = paste0("../Documents/AGG/density0_cond",j,".png"),plot=pl1,height=5,width=8)


pl1 <- pl + geom_line(data=tr1 %>% mutate(Method_optim=factor(Method)) %>% filter(Method %in% c()),aes(x=x,y=y,col=Method)) +
  xlab("Observed residuals kernel smoothed density") + ylab("Density") + ggtitle("Kernel smoothed residual density and AGG fits")+ scale_color_manual(values=c("Normal"="#070707","GenGaus"="#C11432","AGGdelta"="#009ADA","AGGsig"="#66A64F","AGG"="#FDD10A"),drop=FALSE) + ylim(c(0,0.5))
ggsave(filename = paste0("../Documents/AGG/density1_cond",j,".png"),plot=pl1,height=5,width=8)

pl1 <- pl + geom_line(data=tr1 %>% mutate(Method_optim=factor(Method)) %>% filter(Method %in% c("Normal")),aes(x=x,y=y,col=Method)) +
  xlab("Observed residuals kernel smoothed density") + ylab("Density") + ggtitle("Kernel smoothed residual density and AGG fits")+ scale_color_manual(values=c("Normal"="#070707","GenGaus"="#C11432","AGGdelta"="#009ADA","AGGsig"="#66A64F","AGG"="#FDD10A"),drop=FALSE) + ylim(c(0,0.5))
ggsave(filename = paste0("../Documents/AGG/density2_cond",j,".png"),plot=pl1,height=5,width=8)

pl1 <- pl + geom_line(data=tr1 %>% mutate(Method_optim=factor(Method))%>% filter(Method %in% c("Normal","GenGaus")),aes(x=x,y=y,col=Method)) +
  xlab("Observed residuals kernel smoothed density") + ylab("Density") + ggtitle("Kernel smoothed residual density and AGG fits")+ scale_color_manual(values=c("Normal"="#070707","GenGaus"="#C11432","AGGdelta"="#009ADA","AGGsig"="#66A64F","AGG"="#FDD10A")) + ylim(c(0,0.5))
ggsave(filename = paste0("../Documents/AGG/density3_cond",j,".png"),plot=pl1,height=5,width=8)

pl1 <- pl + geom_line(data=tr1 %>% mutate(Method_optim=factor(Method))%>% filter(Method %in% c("Normal","GenGaus","AGGdelta")),aes(x=x,y=y,col=Method)) +
  xlab("Observed residuals kernel smoothed density") + ylab("Density") + ggtitle("Kernel smoothed residual density and AGG fits")+ scale_color_manual(values=c("Normal"="#070707","GenGaus"="#C11432","AGGdelta"="#009ADA","AGGsig"="#66A64F","AGG"="#FDD10A")) + ylim(c(0,0.5))
ggsave(filename = paste0("../Documents/AGG/density4_cond",j,".png"),plot=pl1,height=5,width=8)

pl1 <- pl + geom_line(data=tr1 %>% mutate(Method_optim=factor(Method))%>% filter(Method %in% c("Normal","GenGaus","AGGdelta","AGGsig")),aes(x=x,y=y,col=Method)) +
  xlab("Observed residuals kernel smoothed density") + ylab("Density") + ggtitle("Kernel smoothed residual density and AGG fits")+ scale_color_manual(values=c("Normal"="#070707","GenGaus"="#C11432","AGGdelta"="#009ADA","AGGsig"="#66A64F","AGG"="#FDD10A")) + ylim(c(0,0.5))
ggsave(filename = paste0("../Documents/AGG/density5_cond",j,".png"),plot=pl1,height=5,width=8)

# 5a. explore bootstrap of summer residuals: boxplots of AIC --------
Nrep <- 100
for (i in 2:Nrep) {
  set.seed(12*i)
  # N <- 50000
  # v <- 0.99 # threshold for conditioning variable
  # sims <- generate_Y(N=N) %>% link_log(dep=1/2) %>%
  #   link_log(dep=1/2) %>% link_log(dep=1/2) %>% link_log(dep=1/2) %>%
  #   apply(c(1,2),FUN=frechet_laplace_pit) %>% as.data.frame()
  # pe <-  par_est(df=sims,v=v,given=j,margin = "Normal", method = "sequential2")
  # obsr <- observed_residuals(df=sims,given=j,v=v,a=pe$a,b=pe$b)
  obsr <- apply(obsress,MARGIN = 2,FUN = function(x) {sample(x,replace=TRUE)})
  L1 <- res_margin_par_est(obs_res = obsr,method="Normal")
  L2 <- res_margin_par_est(obs_res = obsr,method="GenGaus")
  L3 <- res_margin_par_est(obs_res = obsr,method="AGGdelta")
  L4 <- res_margin_par_est(obs_res = obsr,method="AGGsig")
  L5 <- res_margin_par_est(obs_res = obsr,method="AGG")
  
  d <- ncol(summer)
  L1 <- L1 %>% mutate(method="Normal",AIC=2*likres+2*2) %>% mutate(res=c(1:d)[-j]) %>% mutate(mu_agg=mu,sigl=sig,sigu=sig,deltal=2,deltau=2) %>% dplyr::select(c(likres,mu_agg,sigl,sigu,deltal,deltau,AIC,method,res))
  L2 <- L2 %>% mutate(method="GenGaus",AIC=2*likres+2*3) %>% mutate(res=c(1:d)[-j]) %>% mutate(sigl=sig_agg,sigu=sig_agg,deltal=delta,deltau=delta) %>% dplyr::select(c(likres,mu_agg,sigl,sigu,deltal,deltau,AIC,method,res))
  L3 <- L3 %>% mutate(method="AGGdelta",AIC=2*likres+2*4) %>% mutate(res=c(1:d)[-j]) %>% mutate(sigl=sig_agg,sigu=sig_agg) %>% dplyr::select(c(likres,mu_agg,sigl,sigu,deltal,deltau,AIC,method,res))
  L4 <- L4 %>% mutate(method="AGGsig",AIC=2*likres+2*4) %>% mutate(res=c(1:d)[-j]) %>% mutate(deltal=delta,deltau=delta) %>% dplyr::select(c(likres,mu_agg,sigl,sigu,deltal,deltau,AIC,method,res))
  L5 <- L5 %>% mutate(method="AGG",AIC=2*likres+2*5) %>% mutate(res=c(1:d)[-j]) %>% dplyr::select(c(likres,mu_agg,sigl,sigu,deltal,deltau,AIC,method,res))
  
  tmp1 <- rbind(L1,L2,L3,L4,L5) %>% mutate(given=j) %>%  mutate(method=factor(method),given=factor(as.character(given)),res=factor(as.character(res)))
tmp <- rbind(tmp,tmp1)
}

pollutant_names <- c("O_3","NO_2","NO","SO_2","PM_{10}")
label_residual_pollutant <- c(TeX("$O_3$"),TeX("$NO_2$"),TeX("$NO$"),TeX("$SO_2$"),TeX("$PM_{10}$"))
ylims <- c(min(tmp$AIC),max(tmp$AIC))
plottitle <- TeX(paste0("Conditioning on $",pollutant_names[j],"> F^{-1}_{",pollutant_names[j],"} (0.7)$"))
# do a series of plots for the presentation
pAIC1 <- ggplot(tmp %>% filter(given==as.character(j)) %>% filter(res %in% as.character(c(1:d)[-j][1]))) + geom_boxplot(aes(y=AIC,x=res,col=method)) +  scale_color_manual(values=c("Normal"="#070707","GenGaus"="#C11432","AGGdelta"="#009ADA","AGGsig"="#66A64F","AGG"="#FDD10A"),drop=FALSE) + xlab("Residuals") + ggtitle(plottitle) + scale_x_discrete(labels=label_residual_pollutant[-j],drop=FALSE) + ylim(ylims)
pAIC2 <- ggplot(tmp %>% filter(given==as.character(j)) %>% filter(res %in% as.character(c(1:d)[-j][1:2]))) + geom_boxplot(aes(y=AIC,x=res,col=method)) +  scale_color_manual(values=c("Normal"="#070707","GenGaus"="#C11432","AGGdelta"="#009ADA","AGGsig"="#66A64F","AGG"="#FDD10A"),drop=FALSE) + xlab("Residuals") + ggtitle(plottitle) + scale_x_discrete(labels=label_residual_pollutant[-j],drop=FALSE) + ylim(ylims)
pAIC3 <- ggplot(tmp %>% filter(given==as.character(j)) %>% filter(res%in%  as.character(c(1:d)[-j][1:3]))) + geom_boxplot(aes(y=AIC,x=res,col=method)) +  scale_color_manual(values=c("Normal"="#070707","GenGaus"="#C11432","AGGdelta"="#009ADA","AGGsig"="#66A64F","AGG"="#FDD10A"),drop=FALSE) + xlab("Residuals") + ggtitle(plottitle) + scale_x_discrete(labels=label_residual_pollutant[-j],drop=FALSE) + ylim(ylims)
pAIC4 <- ggplot(tmp %>% filter(given==as.character(j))) + geom_boxplot(aes(y=AIC,x=res,col=method)) +  scale_color_manual(values=c("Normal"="#070707","GenGaus"="#C11432","AGGdelta"="#009ADA","AGGsig"="#66A64F","AGG"="#FDD10A")) + xlab("Residuals") + ggtitle(plottitle) + scale_x_discrete(labels=label_residual_pollutant[-j],drop=FALSE) + ylim(ylims)
ggsave(filename = paste0("../Documents/AGG/AICZ2_cond",j,".png"),plot=pAIC1,height=5,width=6)
ggsave(filename = paste0("../Documents/AGG/AICZ3_cond",j,".png"),plot=pAIC2,height=5,width=6)
ggsave(filename = paste0("../Documents/AGG/AICZ4_cond",j,".png"),plot=pAIC3,height=5,width=6)
ggsave(filename = paste0("../Documents/AGG/AICZ5_cond",j,".png"),plot=pAIC4,height=5,width=6)

# 5b. likelihood ratio test ------
#ablik <- c()
GenGaus <- c()
AGGsig <- c()
AGGdelta <- c()
AGGsigdelta <- c()
res <- rep(paste0("Z",c(1:d)[-j]),length(L1))
for (i in 1:length(L1)) {
  GenGaus <- tmp %>% filter(method=="GenGaus") %>% pull(likres)
  AGGsig <- tmp %>% filter(method=="AGGsig") %>% pull(likres)  
  AGGdelta <- tmp %>% filter(method=="AGGdelta") %>% pull(likres)
  AGGsigdelta <- tmp %>% filter(method=="AGG") %>% pull(likres)
}

AGG_GG <- 2*((tmp %>% filter(method=="GenGaus",given==as.character(j)) %>% pull(likres))-(tmp %>% filter(method=="AGG",given==as.character(j)) %>% pull(likres)) )
AGGsig_GG <- 2*((tmp %>% filter(method=="GenGaus",given==as.character(j)) %>% pull(likres))-(tmp %>% filter(method=="AGGsig",given==as.character(j)) %>% pull(likres)) )
AGGdelta_GG <- 2*((tmp %>% filter(method=="GenGaus",given==as.character(j)) %>% pull(likres))-(tmp %>% filter(method=="AGGdelta",given==as.character(j)) %>% pull(likres)) )
tmpagg <- data.frame(AGG_GG=AGG_GG,res=tmp %>% filter(method=="GenGaus",given==as.character(j)) %>% pull(res),given=tmp %>% filter(method=="GenGaus",given==as.character(j)) %>% pull(given))
ylims <- c(0,ceiling(max(c(AGG_GG,AGGsig_GG,AGGdelta_GG))))
p1 <- ggplot(tmpagg,aes(x=res,y=AGG_GG))+
  geom_boxplot() + xlab("Residuals") +
  # ylab("Different scale vs shape")
  ylab(TeX("$2(l_{AGG}(\\hat{\\mu},\\hat{\\sigma}_l,\\hat{\\sigma}_u,\\hat{\\delta}_l,\\hat{\\delta}_u)-l_{GG}(\\hat{\\mu},\\hat{\\sigma},\\hat{\\delta}))$")) + ylim(ylims) +
  geom_hline(yintercept = qchisq(0.95,2),linetype="dashed",col = "#66A64F") +  geom_hline(yintercept = qchisq(0.99,2),linetype="dashed", col= "#FDD10A") + scale_x_discrete(labels=label_residual_pollutant[-j],drop=FALSE)

tmpagg <- data.frame(AGG_GG=AGGsig_GG,res=tmp %>% filter(method=="GenGaus",given==as.character(j)) %>% pull(res),given=tmp %>% filter(method=="AGGsig",given==as.character(j)) %>% pull(given))
p2 <- ggplot(tmpagg,aes(x=res,y=AGG_GG))+ 
  geom_boxplot() + xlab("Residuals") +
  # ylab("Different scale vs shape")
  ylab(TeX("$2(l_{AGG}(\\hat{\\mu},\\hat{\\sigma}_l,\\hat{\\sigma}_u,\\hat{\\delta}_c,\\hat{\\delta}_c)-l_{GG}(\\hat{\\mu},\\hat{\\sigma},\\hat{\\delta}))$"))  + ylim(ylims) +
  geom_hline(yintercept = qchisq(0.95,1),linetype="dashed",col = "#66A64F") +  geom_hline(yintercept = qchisq(0.99,1),linetype="dashed", col= "#FDD10A") + scale_x_discrete(labels=label_residual_pollutant[-j],drop=FALSE)

tmpagg <- data.frame(AGG_GG=AGGdelta_GG,res=tmp %>% filter(method=="GenGaus",given==as.character(j)) %>% pull(res),given=tmp %>% filter(method=="GenGaus",given==as.character(j)) %>% pull(given))
p3 <- ggplot(tmpagg,aes(x=res,y=AGG_GG))+
  geom_boxplot() + xlab("Residuals") +
  # ylab("Different scale vs shape")
  ylab(TeX("$2(l_{AGG}(\\hat{\\mu},\\hat{\\sigma}_c,\\hat{\\sigma}_c,\\hat{\\delta}_l,\\hat{\\delta}_u)-l_{GG}(\\hat{\\mu},\\hat{\\sigma},\\hat{\\delta}))$"))+ ylim(ylims) +
  geom_hline(yintercept = qchisq(0.95,1),linetype="dashed",col = "#66A64F") +  geom_hline(yintercept = qchisq(0.99,1),linetype="dashed", col= "#FDD10A") + scale_x_discrete(labels=label_residual_pollutant[-j],drop=FALSE)
p <- grid.arrange(p1,p2,p3,ncol=3)

ggsave(filename = paste0("../Documents/AGG/AGGmargin_likratiotest_cond",j,".png"),plot=p,height=3,width=8)

# 5c. count the minima of each method ----
tmpc <- tmp %>% mutate(ite=rep(1:Nrep,each=20)) %>% group_by(ite,res) %>% slice(which.min(AIC))
p <- ggplot(tmpc) + geom_bar(aes(x=method,fill=res)) +  scale_fill_manual(values=c("#C11432","#009ADA","#66A64F","#FDD10A"),labels=label_residual_pollutant[-j],drop=FALSE) + labs(fill=c("Residual variable"))
ggsave(filename = paste0("../Documents/AGG/countAIC_cond",j,".png"),plot=p,height=3,width=5)

# 5d. combine AGG methods -----
tmpc <- tmpc %>% mutate(method=recode(method,AGG="AGG",AGGsig="AGG",AGGdelta="AGG")) 
p <- ggplot(tmpc) + geom_bar(aes(x=method,fill=res)) +  scale_fill_manual(values=c("#C11432","#009ADA","#66A64F","#FDD10A"),labels=label_residual_pollutant[-j],drop=FALSE) + labs(fill=c("Residual variable"))
ggsave(filename = paste0("../Documents/AGG/countAIC_AGGcond",j,".png"),plot=p,height=3,width=5)

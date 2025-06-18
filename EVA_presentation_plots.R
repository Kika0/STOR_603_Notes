#library(tmap) # spatial map plots
library(gridExtra)
library(sf) # for handling spatial sf objects
library(tidyverse)
library(latex2exp) # latex expressions for plot labels
library(texmex) # for pollution data

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

# 1. plot of Birmingham and Glasgow on original and on Laplace margins
Birm_index <- find_site_index(site=Birmingham,grid_uk = xyUK20_sf)
Glasgow_index <- find_site_index(site=Glasgow,grid_uk = xyUK20_sf)

tmp <- data_mod[,c(Birm_index,Glasgow_index)]
names(tmp) <- names(df_sites)[1:2]
p1 <- ggplot(tmp) + geom_point(aes(x=Birmingham,y=Glasgow),size=0.5)
tmpL <- data_mod_Lap[,c(Birm_index,Glasgow_index)]
names(tmpL) <- names(df_sites)[1:2]
p2 <- ggplot(tmpL) + geom_point(aes(x=Birmingham,y=Glasgow),size=0.5)
p <- grid.arrange(p1,p2,ncol=2)
ggsave(p,filename= "../Documents/Birmingham_Glasgow_original_Laplace.png",width=8,height=4)

# 2. plot illustrating conditioning on a variable
London_index <- find_site_index(site=London,grid_uk = xyUK20_sf)
Inverness_index <- find_site_index(site=Inverness,grid_uk = xyUK20_sf)
tmp <- data_mod_Lap[,c(Birm_index,Glasgow_index,Inverness_index)]
names(tmp) <- names(df_sites)[c(1,3,4)]
q <- 0.9
tmp$tf <- (tmp$London>vL)
vL <- quantile(tmp[,1],q)
p1 <- ggplot(tmp) + geom_point(aes(x=London,y=Birmingham),size=0.5) + xlab(TeX("$Y_1$")) + ylab(TeX("$Y_2$")) + scale_color_manual(values = c("black", "#009ADA")) + theme(legend.position="none") + coord_fixed() + xlim(-10,10) + ylim(-10,10) 
p2 <- ggplot(tmp) + geom_point(aes(x=London,y=Inverness),size=0.5) + xlab(TeX("$Y_1$")) + ylab(TeX("$Y_3$")) + scale_color_manual(values = c("black", "#009ADA")) + theme(legend.position="none") + coord_fixed() + xlim(-10,10) + ylim(-10,10)  
p1n <- ggplot(tmp) + geom_point(aes(x=London,y=Birmingham,col=tf),size=0.5) + xlab(TeX("$Y_1$")) + ylab(TeX("$Y_2$")) + scale_color_manual(values = c("black", "#009ADA")) + theme(legend.position="none") + coord_fixed() + xlim(-10,10) + ylim(-10,10) +geom_vline(xintercept=vL,color="#009ADA",linetype="dashed")
p2n <- ggplot(tmp) + geom_point(aes(x=London,y=Inverness,col=tf),size=0.5) + xlab(TeX("$Y_1$")) + ylab(TeX("$Y_3$")) + scale_color_manual(values = c("black", "#009ADA")) + theme(legend.position="none") + coord_fixed() + xlim(-10,10) + ylim(-10,10) +geom_vline(xintercept=vL,color="#009ADA",linetype="dashed") 
p <- grid.arrange(p1,p2,ncol=2)
ggsave(p,filename="../Documents/condmodel_illustration.png",width=8,height=4)
pn <- grid.arrange(p1n,p2n,ncol=2)
ggsave(pn,filename="../Documents/condmodel_illustrationvL.png",width=8,height = 4)

# 3. update density plot
v <- 0.7
j <- 1
# transform to Laplace margins
summer_lap <- as.data.frame((summer %>% apply(c(2),FUN=row_number))/(nrow(summer)+1)) %>% apply(c(1,2),FUN=unif_laplace_pit) %>% as.data.frame()
pes <-  par_est(df=summer_lap,v=v,given=j,margin = "Normal", method = "sequential2")
obsress <- observed_residuals(df = summer_lap,given = j,v = v,a = pes$a,b=pes$b) 
x <- obsress[,1]

# calculate parameters of each method
margin <- c("GenGaus","AGGdelta","AGGsig","AGG")
tr1 <- data.frame(y=as.numeric(),x=as.numeric(),Method=as.character())
for (i in 1:5) {
  methodpar <- as.numeric(res_margin_par_est(x,method=margin[i]))[-1]
  tr1 <- rbind(tr1,data.frame(y=AGG_density(x=x,theta = methodpar),x=x,Method=margin[i]))
}

pl <- ggplot(data.frame(x=Z2)) + geom_density(aes(x=x),linetype="dashed")
pl1 <- pl + geom_line(data=tr1 %>% mutate(Method_optim=factor(Method)),aes(x=x,y=y,col=Method)) +
  xlab("Observed residuals kernel smoothed density") + ylab("Density") + ggtitle("Kernel smoothed residual density and AGG fits")+ scale_color_manual(values=c("Normal"="#070707","GenGaus"="#C11432","AGGdelta"="#009ADA","AGGsig"="#66A64F","AGG"="#FDD10A")) + ylim(c(0,0.5))
ggsave(filename = "../Documents/AGG/density0.png",plot=pl1,height=5,width=8)


pl1 <- pl + geom_line(data=tr1 %>% mutate(Method_optim=factor(Method)) %>% filter(Method %in% c()),aes(x=x,y=y,col=Method)) +
  xlab("Observed residuals kernel smoothed density") + ylab("Density") + ggtitle("Kernel smoothed residual density and AGG fits")+ scale_color_manual(values=c("Normal"="#070707","GenGaus"="#C11432","AGGdelta"="#009ADA","AGGsig"="#66A64F","AGG"="#FDD10A"),drop=FALSE) + ylim(c(0,0.5))
ggsave(filename = "../Documents/AGG/density1.png",plot=pl1,height=5,width=8)

pl1 <- pl + geom_line(data=tr1 %>% mutate(Method_optim=factor(Method)) %>% filter(Method %in% c("Normal")),aes(x=x,y=y,col=Method)) +
  xlab("Observed residuals kernel smoothed density") + ylab("Density") + ggtitle("Kernel smoothed residual density and AGG fits")+ scale_color_manual(values=c("Normal"="#070707","GenGaus"="#C11432","AGGdelta"="#009ADA","AGGsig"="#66A64F","AGG"="#FDD10A"),drop=FALSE) + ylim(c(0,0.5))
ggsave(filename = "../Documents/AGG/density2.png",plot=pl1,height=5,width=8)

pl1 <- pl + geom_line(data=tr1 %>% mutate(Method_optim=factor(Method))%>% filter(Method %in% c("Normal","GenGaus")),aes(x=x,y=y,col=Method)) +
  xlab("Observed residuals kernel smoothed density") + ylab("Density") + ggtitle("Kernel smoothed residual density and AGG fits")+ scale_color_manual(values=c("Normal"="#070707","GenGaus"="#C11432","AGGdelta"="#009ADA","AGGsig"="#66A64F","AGG"="#FDD10A")) + ylim(c(0,0.5))
ggsave(filename = "../Documents/AGG/density3.png",plot=pl1,height=5,width=8)

pl1 <- pl + geom_line(data=tr1 %>% mutate(Method_optim=factor(Method))%>% filter(Method %in% c("Normal","GenGaus","AGGdelta")),aes(x=x,y=y,col=Method)) +
  xlab("Observed residuals kernel smoothed density") + ylab("Density") + ggtitle("Kernel smoothed residual density and AGG fits")+ scale_color_manual(values=c("Normal"="#070707","GenGaus"="#C11432","AGGdelta"="#009ADA","AGGsig"="#66A64F","AGG"="#FDD10A")) + ylim(c(0,0.5))
ggsave(filename = "../Documents/AGG/density4.png",plot=pl1,height=5,width=8)

# 1. temperature large at Birmingham
data_mod[data_mod]
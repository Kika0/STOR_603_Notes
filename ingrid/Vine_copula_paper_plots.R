library(tidyverse)
library(texmex)

theme_set(theme_bw())
theme_replace(
  panel.spacing = unit(2, "lines"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA) )
file.sources = list.files(pattern="*helpers.R")
sapply(file.sources,source,.GlobalEnv)

# 4. density plot illustration for pollution data residual margins -------
v <- 0.7
j <- 3 # change j to repeat analysis cond. on other variables

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
pl1 <- pl + geom_line(data=tr1 %>% mutate(Method_optim=factor(Method)) %>% filter(Method %in% c("Normal","GenGaus","AGG")),aes(x=x,y=y,col=Method)) +
  xlab("Observed residuals kernel smoothed density") + ylab("Density") + ggtitle("Kernel smoothed residual density and AGG fits")+ scale_color_manual(values=c("Normal"="#070707","GenGaus"="#C11432","AGG"="#009ADA","AGGsig"="#66A64F","AGGdelta"="#FDD10A"),drop=FALSE) + ylim(c(0,0.5))
ggsave(filename = paste0("../Documents/Vinecopula_paper/density1_cond",j,".png"),plot=pl1,height=5,width=8)

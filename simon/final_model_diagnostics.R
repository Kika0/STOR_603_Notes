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

# 1. plot of beta against latitude difference from the conditioning site ------
plot_beta_latitude <- function(b,given,res,gridUK=xyUK20_sf,cond_site,folder_name,plot_name="beta_latitude_difference") {
  stopifnot(length(b)==length(given))
lat_diff <-   sapply(1:length(b),FUN=function(i) {
    as.numeric(gridUK$lat)[res[i]] - as.numeric(gridUK$lat)[given[i]]
  })
tmp <- data.frame("b"=b,"cond_site"=cond_site,"lat_diff"=lat_diff)
c12 <- c(
  "#009ADA", "#C11432", # red
           "green4",
           "#6A3D9A", # purple
           "#FF7F00", # orange
           "black", "gold1",
           "#FB9A99", # lt pink
           "gray70", 
           "darkturquoise", "green1", 
           "darkorange4"
)
p <- ggplot(tmp) + geom_point(aes(x=lat_diff,y=b,col=cond_site),size=0.3) + facet_wrap(cond_site) + scale_color_manual(values=c12) + labs(x="Latitude difference",y=TeX("$\\beta$"),col="Conditioning site")
ggsave(p,filename=paste0(folder_name,plot_name,".png"),width=8,height=5)
}

plot_beta_latitude(b=est_all_sf$b,given=est_all_sf$given,res=est_all_sf$res,cond_site = est_all_sf$cond_site,folder_name = folder_name)

# try with new data
load("data_processed/final_model_3_parameter_estimates.RData",verbose=TRUE)
for (i in 1:13) {
  print(par_est_model_3[[i]][[4]])
}

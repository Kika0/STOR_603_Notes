# load the data
load("../kristina/ukgd_cpm85_5k_x100y23_MSdata01.RData")
data01.std.param # not much info here
data01 <- data01 %>% mutate(year=floor(time))
dim(data01)
summary(data01)
glimpse(data01)
data01 %>% group_by(year) %>% summarize(mean_year=mean(x)) %>% ggplot(aes(x=year,y=mean_year)) + geom_point()
ggplot(data01 %>% filter(doy %in% c(1:20))) + geom_point(aes(x=time,y=x))
# to explore gpd functions for the tails
str(gpdpar)
summary(chosen.MSgpd)

# load the Gpd parameters
load("../kristina/MSGpdParam/ukgd_cpm85_5k_x100y23.MSGpdParam.2024-10-22-064747.RData")
glimpse(gpdpar) # dataframe for columns for scale, shape and threshold

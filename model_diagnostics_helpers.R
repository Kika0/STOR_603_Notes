PP_plot <- function(observed,simulated,title=NULL) {
# first do for one dimension
X1 <- observed
X2 <- simulated
# add test to check for same length
if (length(X1)!=length(X2)) {
  stop("Length of observed and simulated data must be the same")
}
X <- c(X1,X2)
x <- sort(X) # sequence to evaluate emp distribution estimates
a1 <- a2 <- c()
for (i in 1:length(x)) {
  a1[i] <- sum(X1<x[i]+10^(-8))/length(X1)
  a2[i] <- sum(X2<x[i]+10^(-8))/length(X2)
}
# add also 95% bootstrap interval
bf1 <- data.frame(x=x)
bf2 <- data.frame(x=x)
for (i in 1:1000) {
  p1 <- p2 <- c()
  # sample data
  Y1 <- sample(size=length(X1),x=X,replace = FALSE)
  Y2 <- lubridate::setdiff(X,Y1)
  # calculate empirical probabilities
  for (i in 1:length(x)) {
    p1[i] <- sum(Y1<x[i]+10^(-8))/length(Y1)
    p2[i] <- sum(Y2<x[i]+10^(-8))/length(Y2)
  }
  bf1 <- cbind(bf1,p1)
  bf2 <- cbind(bf2,p2)
}
bf1num <- as.numeric(unlist(bf1[,2:1001]))
bf2num <- as.numeric(unlist(bf2[,2:1001]))
u1 <- l1 <-  seq(0,1,length.out=length(X1)+1)
u2 <- l2 <- c()
for (i in 1:length(u1)) {
  u2[i] <- quantile(bf2num[bf1num==round(u1[i],2)],p=0.975)
  l2[i] <- quantile(bf2num[bf1num==round(u1[i],2)],p=0.025)
}

df <- data.frame(x=a1,y=a2)
dfCI <- data.frame(u1=u1,u2=u2,l1=l1,l2=l2)
  pp <- df %>% ggplot()  + 
    # geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), color = "#C11432",size=1.2) +
    # geom_point(aes(x=Empirical,y=Model)) + 
    # # geom_point(data=dftmp,mapping=aes(x=x,y=y),col="#C11432") +
    geom_line(data=dfCI,aes(x=l1,y=l2),linetype="dashed", col="#C11432") +
    geom_line(data=dfCI,aes(x=u1,y=u2),linetype="dashed", col="#C11432") +
    geom_line(aes(x=x,y=y), col="black")  +
    geom_ribbon(data=dfCI,aes(x=u1,ymin=l2,ymax=u2), fill="#C11432", alpha=0.2) +
    ggtitle(title) + 
    xlab("Empirical") + ylab("Model") + coord_fixed()
  return(pp)
}

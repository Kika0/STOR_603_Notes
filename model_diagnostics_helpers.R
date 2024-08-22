PP_plot <- function(observed,simulated) {
# first do for one dimension
X1 <- observed
X2 <- simulated
# add test to check for same length
if (length(X1)!=length(X2)) {
  stop("Length of observed and simulated data must be the same")
}
X <- c(X1,X2)
x <- sort(X) # sequence to evaluate emp distribution estimates
p1 <- p2 <- c()
p11 <- p22 <- c()
Y1 <- sample(size=length(X1),x=X,replace = FALSE)
Y2 <- lubridate::setdiff(X,Y1)
for (i in 1:length(x)) {
  p11[i] <- sum(Y1<x[i])/length(Y1)
  p22[i] <- sum(Y2<x[i])/length(Y2)
}
p1 <- append(p1,p11)
p2 <- append(p2,p22)
a1 <- a2 <- c()
for (i in 1:length(x)) {
  a1[i] <- sum(X1<x[i])/length(X1)
  a2[i] <- sum(X2<x[i])/length(X2)
}
df <- data.frame(x=a1,y=a2,x1=p11,y1=p22)
  pp <- df %>% ggplot()  + 
    # geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), color = "#C11432",size=1.2) +
    # geom_point(aes(x=Empirical,y=Model)) + 
    # # geom_point(data=dftmp,mapping=aes(x=x,y=y),col="#C11432") +
    geom_line(aes(x=x1,y=y1),linetype="dashed", col="#C11432") +
    geom_line(aes(x=x,y=y), linetype="dashed", col="black")  +
    # geom_ribbon(data=dfCI,aes(x=Empirical,ymin=Uup,ymax=Ulow), fill="#C11432", alpha=0.2) +
    ggtitle("Probability Plot") + xlab("Empirical") + ylab("Model")
  return(pp)
}

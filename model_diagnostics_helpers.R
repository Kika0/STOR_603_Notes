PP_plot <- function(observed,simulated,title=NULL,CIcol=NULL,Uup=NULL,Ulow=NULL,tol_bounds="bootstrap") {
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
N <- length(X1)
for (i in 1:length(x)) {
  a1[i] <- sum(X1<x[i]+10^(-8))/N
  a2[i] <- sum(X2<x[i]+10^(-8))/N
}
# add also 95% bootstrap interval
u1 <- l1 <-  1:N/(N+1)
# special case for custom tolerance bounds
if (is.numeric(Uup) & is.numeric(Ulow) ) {
  # add test to check for same length
  if (N!=length(Uup) & N!=length(Ulow)) {
    stop("Length of tolerance bounds must be the same as the data")
  }
  
} else if (tol_bounds=="bootstrap") {
bf1 <- data.frame(x=x)
bf2 <- data.frame(x=x)
Nboot <- 1000
for (i in 1:Nboot) {
  p1 <- p2 <- c()
  # sample data
  Y1 <- sample(size=N,x=X,replace = FALSE)
  Y2 <- lubridate::setdiff(X,Y1)
  # calculate vectors of empirical probabilities
  for (i in 1:(2*N)) {
    p1[i] <- sum(Y1<x[i]+10^(-8))/(N+1)
    p2[i] <- sum(Y2<x[i]+10^(-8))/(N+1)
  }
  bf1 <- cbind(bf1,p1)
  bf2 <- cbind(bf2,p2)
}
bf1num <- as.numeric(unlist(bf1[,2:(Nboot+1)]))
bf2num <- as.numeric(unlist(bf2[,2:(Nboot+1)]))
Uup <- Ulow <- c()
for (i in 1:N) {
  Uup[i] <- quantile(bf2num[round(bf1num,5)==round(u1[i],5)],p=0.975)
  Ulow[i] <- quantile(bf2num[round(bf1num,5)==round(u1[i],5)],p=0.025)
}
} else if (tol_bounds=="beta_dist") {
  Uup <- Ulow <- c()
  Uup <- sapply(1:N, function(i){qbeta(0.975, i, N+1-i)})
  Ulow <- sapply(1:N, function(i){qbeta(0.025, i, N+1-i)})
}
df <- data.frame(x=X1,y=X2)
if (tol_bounds!="none") {
dfCI <- data.frame(u1=u1,u2=Uup,l1=l1,l2=Ulow)
} else {
  dfCI <- data.frame(u1=Uup,u2=u1,l1=Ulow,l2=l1)
}
if (is.null(CIcol)) {
  fillcol <- "#C11432"
} else { fillcol <- CIcol}
  pp <- ggplot()  + 
    # geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), color = fillcol,size=1.2) +
    # geom_point(aes(x=Empirical,y=Model)) + 
    # # geom_point(data=dftmp,mapping=aes(x=x,y=y),col=fillcol) +
    geom_line(data=dfCI,aes(x=l1,y=l2),linetype="dashed", col=fillcol) +
    geom_line(data=dfCI,aes(x=u1,y=u2),linetype="dashed", col=fillcol) +
    geom_point(data=df,aes(x=x,y=y), col="black",size=0.2)  +
    geom_ribbon(data=dfCI,aes(x=u1,ymin=l2,ymax=u2), fill=fillcol, alpha=0.2) +
    geom_segment(data=data.frame(x1=0,x2=1,y1=0,y2=1),mapping=aes(x=x1,y=y1,xend=x2,yend=y2),linetype="dashed") +
    ggtitle(title) + 
    xlab("Model") + ylab("Empirical") + coord_fixed()
  return(pp)
}

QQ_plot <- function(observed,simulated,title=NULL,CIcol=NULL,Uup=NULL,Ulow=NULL,tol_bounds="bootstrap") {
  # first do for one dimension
  X1 <- observed
  X2 <- simulated
  # add test to check for same length
  if (length(X1)!=length(X2)) {
    stop("Length of observed and simulated data must be the same")
  } 
  X <- c(X1,X2)
  x <- sort(X) # sequence to evaluate empirical distribution estimates
  N <- length(X1)
  # add also 95% bootstrap interval
  u1 <- l1 <-  sort(X1)
  # special case for custom tolerance bounds
  if (is.numeric(Uup) & is.numeric(Ulow) ) {
    # add test to check for same length
    if (N!=length(Uup) & N!=length(Ulow)) {
      stop("Length of tolerance bounds must be the same as the data")
    }
    
  } else if (tol_bounds=="bootstrap") {
    bf1 <- data.frame(x=x)
    bf2 <- data.frame(x=x)
    Nboot <- 1000
    for (i in 1:Nboot) {
      p1 <- p2 <- c()
      # sample data
      Y1 <- sort(sample(size=N,x=X,replace = FALSE)) # y-axis
      Y2 <- sort(lubridate::setdiff(X,Y1)) # x-axis
      bf1 <- cbind(bf1,Y1) # y-axis
      bf2 <- cbind(bf2,Y2) # x-axis
    }
    bf1num <- as.numeric(unlist(bf1[,2:(Nboot+1)])) # y-axis
    bf2num <- as.numeric(unlist(bf2[,2:(Nboot+1)])) # x-axis
    Uup <- Ulow <- c()
    for (i in 1:N) {
      # tolerance bounds are on x-axis
      Uup[i] <- quantile(bf2num[round(bf1num,5)==round(u1[i],5)],p=0.975)
      Ulow[i] <- quantile(bf2num[round(bf1num,5)==round(u1[i],5)],p=0.025)
    }
  } 
  df <- data.frame(x=X2,y=X1)
  if (tol_bounds!="none") {
    dfCI <- data.frame(u1=u1,u2=Uup,l1=l1,l2=Ulow)
  } else {
    dfCI <- data.frame(u1=Uup,u2=u1,l1=Ulow,l2=l1)
  }
  if (is.null(CIcol)) {
    fillcol <- "#C11432"
  } else { fillcol <- CIcol}
  Qmin <- min(X1,X2)-0.2
  Qmax <- max(X1,X2)+0.2
  
  qq <- ggplot()  + 
    # geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), color = fillcol,size=1.2) +
    # geom_point(aes(x=Empirical,y=Model)) + 
    # # geom_point(data=dftmp,mapping=aes(x=x,y=y),col=fillcol) +
    geom_line(data=dfCI,aes(x=l1,y=l2),linetype="dashed", col=fillcol) +
    geom_line(data=dfCI,aes(x=u1,y=u2),linetype="dashed", col=fillcol) +
    geom_point(data=df,aes(x=x,y=y), col="black",size=0.2)  +
    geom_ribbon(data=dfCI,aes(x=u1,ymin=l2,ymax=u2), fill=fillcol, alpha=0.2) +
    geom_segment(data=data.frame(x1=Qmin,x2=Qmax,y1=Qmin,y2=Qmax),mapping=aes(x=x1,y=y1,xend=x2,yend=y2),linetype="dashed") +
  ggtitle(title) + 
    xlab("Model") + ylab("Empirical") + coord_equal()
  return(qq)
}


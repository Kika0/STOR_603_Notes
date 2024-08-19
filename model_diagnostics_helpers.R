PP_plot <- function(observed,simulated,CImethod = "bootstrap",CIl=NULL,CIu=NULL) {
  N <- length(unlist(observed))
  Empirical <- sort(rank(as.numeric(unlist(observed))))/(N+1)
  Model <- sort(rank(as.numeric(unlist(simulated))))/(N+1)
  df <- data.frame(Emprirical=Empirical,Model=Model)
  if (CImethod == "bootstap") {
    emp <- Empirical
    mod <- sort(rank(as.numeric(sample(unlist(observed),size=N))))/(N+1)
    dftmp <- data.frame(x=emp,y=mod)
  }
  # calculate upper and lower bounds for 95% CI to add to the plots
  Ulow<-c(0,sapply(Empirical*(N+1), function(i){qbeta(0.025, i, N+1-i)}),1)
  Uup<-c(0,sapply(Empirical*(N+1), function(i){qbeta(0.975, i, N+1-i)}),1)
  dfCI <- data.frame(Empirical=c(0,Empirical,1),Ulow=Ulow,Uup=Uup)
  pp <- df %>% ggplot()  + 
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), color = "#C11432",size=1.2) +
    geom_point(aes(x=Empirical,y=Model)) + 
    # geom_point(data=dftmp,mapping=aes(x=x,y=y),col="#C11432") +
    geom_line(data=dfCI,aes(x=Empirical,y=Uup),linetype="dashed", col="#C11432") +
    geom_line(data=dfCI,aes(x=Empirical,y=Ulow), linetype="dashed", col="#C11432")  +
    geom_ribbon(data=dfCI,aes(x=Empirical,ymin=Uup,ymax=Ulow), fill="#C11432", alpha=0.2) +
    ggtitle("Probability Plot") + xlab("Empirical") + ylab("Model")
  return(pp)
}

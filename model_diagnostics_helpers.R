PP_plot <- function(observed,simulated,CImethod = "bootstrap") {
  N <- length(unlist(observed))
  Empirical <- rank(as.numeric(unlist(observed)))/(N+1)
  Model <- rank(as.numeric(unlist(simulated)))/(N+1)
  df <- data.frame(Emprirical=Empirical,Model=Model)
  pp <- df %>% ggplot()  + 
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), color = "#C11432",size=1.2) +
    geom_point(aes(x=Empirical,y=Model)) + 
    # geom_line(aes(x=Empirical,y=Uup),linetype="dashed", color="#C11432") +
    # geom_line(aes(x=Empirical,y=Ulow), linetype="dashed", color="#C11432")  + 
    # geom_ribbon(aes(x=Empirical,ymin=Uup,ymax=Ulow), fill="#C11432", alpha=0.2) +
    ggtitle("Probability Plot") + xlab("Empirical") + ylab("Model")
  return(pp)
}

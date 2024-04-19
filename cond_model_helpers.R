#' Sample from trivariate logistic EVD
#' 
#' Generate from trivariate logistic EVD distribution with X-Y-Z dependence
#' 
#' @param N A number of generated points.
#' @param dep A vector of dependence parameters for each link.
#' @return A dataframe with 3 columns and N rows.
#' @export
#' @examples
#' generate_dep_X_Y_Y_Z(N=5000,dep=c(1/2,1/2))
generate_dep_X_Y_Y_Z <- function(N,dep=c(1/2,1/2)) {
  
  dat <- tibble(x1=numeric(),x2=numeric(),x3=numeric(),R=numeric(),u=numeric(),q=numeric(),w1=numeric(),w2=numeric(),w3=numeric())
  # set.seed(12)
  x_y <- evd::rbvevd(N,dep=dep[1],model="log")
  x <- exp(x_y[,1])
  Y <- exp(x_y[,2])
  a <- dep[2]
  # generate z
  to_opt <- function(z) {
    (  (  y^(-(1/a)+1)*(y^(-1/a)+z^(-1/a))^(a-1)*exp(-(y^(-1/a)+z^(-1/a))^a)*exp(1/y)  )-Unif)^2
  }
  z <- c()
  for (i in 1:nrow(x_y)){
    Unif <- runif(1) # generate U
    y <- Y[i]
    z[i] <- optim(par=1,fn=to_opt,lower=0,upper=10^6,method="Brent")$par
  }
  sims <- data.frame(X_1=x,X_2=Y,X_3=z)
  return(sims)
}

generate_dep_log_norm <- function(N,dep=c(1/2,1/2)) {
  dat <- tibble(x1=numeric(),x2=numeric(),x3=numeric(),R=numeric(),u=numeric(),q=numeric(),w1=numeric(),w2=numeric(),w3=numeric())
  # set.seed(12)
  x_y <- evd::rbvevd(N,dep=dep[1],model="log")
  x <- exp(x_y[,1])
  Y <- exp(x_y[,2])
  r <- dep[2]
  # generate z
  U <- runif(N)
  y_N <- qnorm(exp(-1/Y))
  z_N <- qnorm(1-U)*(1-r^2)^(1/2) + r*y_N
  z <- (-log(pnorm(z_N)))^(-1)
  sims <- data.frame(X_1=x,X_2=Y,X_3=z)
  return(sims)
}

generate_dep_norm_norm <- function(N,dep=c(1/2,1/2)) {
  dat <- tibble(x1=numeric(),x2=numeric(),x3=numeric(),R=numeric(),u=numeric(),q=numeric(),w1=numeric(),w2=numeric(),w3=numeric())
  # set.seed(12)
  x_y <- mvrnorm(n=N,mu=c(0,0),Sigma = matrix(c(1,dep[1],dep[1],1),ncol=2))
  x <- (-log(pnorm(x_y[,1])))^(-1)
  Y <- (-log(pnorm(x_y[,2])))^(-1)
  r <- dep[2]
  # generate z
  U <- runif(N)
  y_N <- qnorm(exp(-1/Y))
  z_N <- qnorm(1-U)*(1-r^2)^(1/2) + r*y_N
  z <- (-log(pnorm(z_N)))^(-1)
  sims <- data.frame(X_1=x,X_2=Y,X_3=z)
  return(sims)
}

# generate Y from standard Fréchet
generate_Y <- function(N) {
U <- runif(N)
Y <- (-log(U))^(-1)
sims <- data.frame(Y1=Y)
  return(sims)
}

link_log <- function(sims,dep=1/2) {
  Y <- sims %>% dplyr::select("Y1") %>% pull()
  N <- length(Y)
  a <- dep
  # generate z
  to_opt <- function(x) {
    (  (  y^(-(1/a)+1)*(y^(-1/a)+x^(-1/a))^(a-1)*exp(-(y^(-1/a)+x^(-1/a))^a)*exp(1/y)  )-Unif)^2
  }
  x <- c()
  for (i in 1:N){
    Unif <- runif(1) # generate U
    y <- Y[i]
    x[i] <- optim(par=1,fn=to_opt,lower=0,upper=10^6,method="Brent")$par
  }
  sims <- sims %>% mutate(X=x)
  names(sims)[-1] <- paste0("Y",ncol(sims))
  return(sims)
}

link_norm <- function(sims,dep=1/2) {
  Y <- sims %>% dplyr::select("Y") %>% pull()
  N <- length(Y1)
  r <- dep
  # generate x
  U <- runif(N)
  y_N <- qnorm(exp(-1/Y))
  x_N <- qnorm(1-U)*(1-r^2)^(1/2) + r*y_N
  x <- (-log(pnorm(x_N)))^(-1)
  sims <- sims %>% mutate(Z=x)
  return(sims)
}



# transform from Fréchet to Laplace
frechet_laplace_pit <- function(x) {
  if (exp(-1/x)<0.5) {
    y <-log(2*exp(-1/x))
  }
  else {
    y <--log(2*(1-exp(-1/x)))
  }
  return(y)
}

# transform back from Laplace to Fréchet margins
laplace_frechet_pit <- function(y) {
  if (y<0) {
    x <- (log(2)-y)^(-1)
  }
  else {
    x <- -(log(1-exp(-(y+log(2)))))^(-1)
  }
  return(x)
}

# transform from normal to Laplace margins
norm_laplace_pit <-  function(x) {
  if (pnorm(x)<0.5) {
    y <-log(2*pnorm(x))
  }
  else {
    y <--log(2*(1-pnorm(x)))
  }
  return(y)
}

# conditioning on Y_1 being extreme to model Y_2 (conditional model in bivariate case)
# calculates MLE a_hat, b_hat, mu_hat and sig_hat
Y_likelihood <- function(theta,df=Y_given_1_extreme,given=1,sim=2) {
  a <- theta[1]
  b <- theta[2]
  mu <- theta[3]
  sig <- theta[4]
  Y1 <- df %>% dplyr::select(starts_with("Y") & contains(as.character(given))) %>% pull()
  Y2 <- df %>% dplyr::select(starts_with("Y") & contains(as.character(sim))) %>% pull()
  #lik <-  prod(1/(Y_1^b *sig)*exp(-(Y_2-a*Y_1-mu*Y_1^b)^2/(2*(Y_1^b*sig)^2)) )
  if (a<(-1) | a>1 | b<0 | b>=1) {
    log_lik <- (-10^6)
  }
  else {
  log_lik <- sum(-log(Y1^b *sig) + (-(Y2-a*Y1-mu*Y1^b)^2/(2*(Y1^b*sig)^2))  )
  }
  return(log_lik)
}

Y_likelihood_initial <- function(theta,df=Y_given_1_extreme,given=1,sim=2) {
  a <- theta[1]
  b <- 0
  mu <- theta[2]
  sig <- theta[3]
  Y1 <- df %>% dplyr::select(starts_with("Y") & contains(as.character(given))) %>% pull()
  Y2 <- df %>% dplyr::select(starts_with("Y") & contains(as.character(sim))) %>% pull()
  #lik <-  prod(1/(Y_1^b *sig)*exp(-(Y_2-a*Y_1-mu*Y_1^b)^2/(2*(Y_1^b*sig)^2)) )
  if (a<(-1) | a>1 | b<0 | b>=1) {
    log_lik <- (-10^6)
  }
  else {
    log_lik <- sum(-log(Y1^b *sig) + (-(Y2-a*Y1-mu*Y1^b)^2/(2*(Y1^b*sig)^2))  )
  }
  return(log_lik)
}

Y_likelihood_fix_ab <- function(theta,a=1,b=0,df=Y_given_1_extreme,given=1,sim=2) {
  mu <- theta[1]
  sig <- theta[2]
  Y1 <- df %>% dplyr::select(starts_with("Y") & contains(as.character(given))) %>% pull()
  Y2 <- df %>% dplyr::select(starts_with("Y") & contains(as.character(sim))) %>% pull()
    log_lik <- sum(-log(Y1^b *sig) + (-(Y2-a*Y1-mu*Y1^b)^2/(2*(Y1^b*sig)^2))  )
  return(log_lik)
}

Y_likelihood_constrained <- function(theta,df=Y_given_1_extreme,given=1,sim=2) {
 v <- 0.99
   a <- theta[1]
  b <- theta[2]
  mu <- theta[3]
  sig <- theta[4]
  Y1 <- df %>% dplyr::select(starts_with("Y") & contains(as.character(given))) %>% pull()
  Y2 <- df %>% dplyr::select(starts_with("Y") & contains(as.character(sim))) %>% pull()
  #lik <-  prod(1/(Y_1^b *sig)*exp(-(Y_2-a*Y_1-mu*Y_1^b)^2/(2*(Y_1^b*sig)^2)) )
  # positive residual quantile
  q <- 0.9
  # zp <- quantile(Y2-Y1,q)
  zp <- max(Y2-Y1)
  # residual quantile
  # z <- quantile( (Y1-a*Y1)/(Y1^b),q)
  z <- max((Y1-a*Y1)/(Y1^b))
  # negative residual quantile
  # zn <- quantile(Y2+Y1,q)
  zn <- max(Y2+Y1)
  
  if (  ( (a<=min(1,1-b*z*v^(b-1),1-v^(b-1)*z+v^(-1)*zp) )|(  ((1-b*z*v^(b-1))<a & a<=1) & (( (1-b^(-1))*(b*z)^(1/(1-b)) *(1-a)^(-b/(1-b)) +zp)>0  ) ) )&
        ( (a<=min(1,1+b*z*v^(b-1),1+v^(b-1)*z-v^(-1)*zn) )|(  ((1+b*z*v^(b-1))<(-a) & (-a)<=1) & (( (1-b^(-1))*(-b*z)^(1/(1-b)) *(1+a)^(-b/(1-b)) -zn)>0  ) ) ) 
        ) {
    
    log_lik <- sum(-log(Y1^b *sig) + (-(Y2-a*Y1-mu*Y1^b)^2/(2*(Y1^b*sig)^2))  )
  }
  else {
    log_lik <- (-10^6)
  }
  return(log_lik)
}



# generate from the empirical distribution of the residuals
F_smooth_Z <- function(Z) {
  Z_smooth <- c()
  for (i in seq_along(Z)) {
    z <- Z[i]
   Z_smooth[i] <-  mean(pnorm((z-Z)/density(Z)$bw))
  }
 return(Z_smooth) 
}

# transform from normal back to kernel smoothed by minimising the square difference
norm_to_orig <- function(Z_N,emp_res) {
  Z <- data.frame(matrix(ncol=ncol(Z_N),nrow=nrow(Z_N)))
  s <- seq(0.02,0.98,length.out=49)
  Zs <- data.frame(matrix(ncol=ncol(Z_N),nrow=length(s)))
  # optimise for these 49 values of s
  to_opt <- function(z) {
    return( (mean(pnorm((z-emp_res[,i])/density(emp_res[,i])$bw)) - s[j])^2)
  }
  for (i in 1:ncol(Z_N)) {
  for (j in 1:length(s)) {
    Zs[j,i] <- optim(fn=to_opt,par=1)$par
  }
  }

 to_opt <- function(z) {
 return( (mean(pnorm((z-emp_res[,i])/density(emp_res[,i])$bw)) - pnorm(Z_N[j,i]))^2)
 }
 
for (i in 1:ncol(Z_N)) {
  for (j in 1:nrow(Z_N)) {
    # optimise cdf using 49 linear segments rather than optimising all directly
    if (pnorm(Z_N[j,i])< min(s) | pnorm(Z_N[j,i])>= max(s) ) {
    Z[j,i] <- optim(fn=to_opt,par=1)$par
    }
    else {
      k <- which.min(pnorm(Z_N[j,i])>s)-1
      a <- 0.02/(Zs[k+1,i]-Zs[k,i])
      b <- -a*Zs[k,i]+s[k]
      Z[j,i] <- (pnorm(Z_N[j,i])-b)/(a)
    }
  }
}
  return(Z)
}

# generate a table of a,b,mu,sig,rho estimates given each of the variables
par_summary <- function(sims,v=0.9) {
  df <- sims %>% dplyr::select(starts_with("Y"))
  par_sum <- data.frame(matrix(nrow=6,ncol=0))
  par_sum_init <- data.frame(matrix(nrow=5,ncol=0))
  # Y_not_1_extreme <- df %>% filter(Y_1<quantile(Y_1,v))
  Z_2 <- c()
  Z_3 <- c()
  Z_N_2 <- c()
  Z_N_3 <- c()
  given <- c()
  d <- ncol(df)
  for (j in 1:ncol(df)) {
    Y_given_1_extreme <- df %>% filter(df[,j]>quantile(df[,j],v))
    res <- c(1:d)[-j]
    opt <- list()
    init_par <- list()
    init_lik <- c()
  for (i in 2:d) {
    # get initial parameters
    init_opt <- optim(par=c(0.5,0,1), fn=Y_likelihood_initial,df=Y_given_1_extreme,given=j,sim=res[i-1],control = list(fnscale=-1))
    init_lik <- append(init_lik,init_opt$value)
    init_par[[i-1]] <- c(init_opt$par[1],0.2,init_opt$par[2],init_opt$par[3])
    # optimise using the initial parameters
    opt[[i-1]] <- optim(par=init_par[[i-1]],fn = Y_likelihood,df=Y_given_1_extreme,given=j,sim=res[i-1],control = list(fnscale=-1))
  }
  a_hat <- c(opt[[1]]$par[1],opt[[2]]$par[1])
  b_hat <- c(opt[[1]]$par[2],opt[[2]]$par[2])
  mu_hat <- c(opt[[1]]$par[3],opt[[2]]$par[3])
  sig_hat <- c(opt[[1]]$par[4],opt[[2]]$par[4])
  lik <-c(opt[[1]]$value,opt[[2]]$value)
  
  a_hat_init <- c(init_par[[1]][1],init_par[[2]][1])
  b_hat_init <- c(0.2,0.2)
  mu_hat_init <- c(init_par[[1]][3],init_par[[2]][3])
  sig_hat_init <- c(init_par[[1]][4],init_par[[2]][4])
  
  # generate residual Z ----
  Y1 <- Y_given_1_extreme[,j]
  Y2 <- Y_given_1_extreme[,res[1]]
  Y3 <- Y_given_1_extreme[,res[2]]
  
  tmp_z2 <- (Y2-a_hat[1]*Y1)/(Y1^b_hat[1])
  tmp_z3 <- (Y3-a_hat[2]*Y1)/(Y1^b_hat[2])
  
  Z_2 <- append(Z_2,tmp_z2)
  Z_3 <- append(Z_3,tmp_z3)
  given <- append(given,rep(j,(N*(1-v)) ))
  
  # calculate the normal using the PIT
  tmp_zn2 <- qnorm(F_smooth_Z(tmp_z2))
  tmp_zn3 <- qnorm(F_smooth_Z(tmp_z3))
  Z_N_2 <- append(Z_N_2,tmp_zn2)
  Z_N_3 <- append(Z_N_3,tmp_zn3)
  
  rho_hat <- cor(tmp_zn2,tmp_zn3)
  par_sum <- cbind(par_sum,data.frame(matrix(round(c(lik,a_hat,b_hat,mu_hat,sig_hat,rep(rho_hat,2)),3),nrow=6,ncol=2,byrow=TRUE)))
  par_sum_init <- cbind(par_sum_init,data.frame(matrix(round(c(init_lik,a_hat_init,b_hat_init,mu_hat_init,sig_hat_init),3),nrow=5,ncol=2,byrow=TRUE)))
  
  }
  return(rbind(par_sum_init,par_sum))
}

# generate a table of a,b,mu,sig,rho estimates given each of the variables
par_est <- function(sims,v=0.9) {
  df <- sims %>% dplyr::select(starts_with("Y"))
  #par_sum <- data.frame(lik=numeric(),a=numeric(),b=numeric(),mu=numeric(),sig=numeric(),given=as.character())
  # Y_not_1_extreme <- df %>% filter(Y_1<quantile(Y_1,v))
  given <- c()
  lik <- c()
  a_hat <- c()
  b_hat <- c()
  mu_hat <- c()
  sig_hat <- c()
  res_var <- c()
  d <- ncol(df)
  for (j in 1:d) {
    Y_given_1_extreme <- df %>% filter(df[,j]>quantile(df[,j],v))
    res <- c(1:d)[-j]
    init_par <- c()
    init_lik <- c()
    for (i in 2:d) {
      # optimise using the initial parameters
      init_opt <- optim(par=c(0.5,0,1), fn=Y_likelihood_initial,df=Y_given_1_extreme,given=j,sim=res[i-1],control = list(fnscale=-1))
      init_par <- c(init_opt$par[1],0.2,init_opt$par[2],init_opt$par[3])
      opt <- optim(par=init_par,fn = Y_likelihood,df=Y_given_1_extreme,given=j,sim=res[i-1],control = list(fnscale=-1))
      a_hat <- append(a_hat,opt$par[1])
      b_hat <- append(b_hat,opt$par[2])
      mu_hat <- append(mu_hat,opt$par[3])
      sig_hat <- append(sig_hat,opt$par[4])
      lik <- append(lik,opt$value)
      res_var <- append(res_var,res[i-1])
    }
  }
  par_sum <- data.frame("lik" = lik, "a" = a_hat, "b" = b_hat, "mu" = mu_hat, "sig" = sig_hat, "given" = rep(1:d,each=(d-1)), "res" = res_var)
  return(par_sum)
}

plot_residual <- function(sims,v=0.99) {
  df <- sims %>% dplyr::select(starts_with("Y"))
  par_sum <- data.frame(matrix(nrow=5,ncol=0))
  
  # Y_not_1_extreme <- df %>% filter(Y_1<quantile(Y_1,v))
  Z_2 <- c()
  Z_3 <- c()
  Z_N_2 <- c()
  Z_N_3 <- c()
  given <- c()
  d <- ncol(df)
  for (j in 1:ncol(df)) {
    Y_given_1_extreme <- df %>% filter(df[,j]>quantile(df[,j],v))
    res <- c(1:d)[-j]
    opt <- list()
    for (i in 2:d) {
      # get initial parameters
      init_opt <- optim(par=c(1,0,1), fn=Y_likelihood_initial,df=Y_given_1_extreme,given=j,sim=res[i-1],control = list(fnscale=-1))$par
      init_par <- c(init_opt[1],0.2,init_opt[2],init_opt[3])
      # optimise using the initial parameters
      opt[[i-1]] <- optim(par=init_par,fn = Y_likelihood,df=Y_given_1_extreme,given=j,sim=res[i-1],control = list(fnscale=-1))
    }
    a_hat <- c(opt[[1]]$par[1],opt[[2]]$par[1])
    b_hat <- c(opt[[1]]$par[2],opt[[2]]$par[2])
    mu_hat <- c(opt[[1]]$par[3],opt[[2]]$par[3])
    sig_hat <- c(opt[[1]]$par[4],opt[[2]]$par[4])
    
    # generate residual Z ----
    Y_1 <- Y_given_1_extreme[,j]
    Y_2 <- Y_given_1_extreme[,res[1]]
    Y_3 <- Y_given_1_extreme[,res[2]]
    
   tmp_z2 <- (Y_2-a_hat[1]*Y_1)/(Y_1^b_hat[1])
   tmp_z3 <- (Y_3-a_hat[2]*Y_1)/(Y_1^b_hat[2])
    
   Z_2 <- append(Z_2,tmp_z2)
   Z_3 <- append(Z_3,tmp_z3)
   given <- append(given,rep(j,50))
    
    # calculate the normal using the PIT
    Z_N_2 <- append(Z_N_2,qnorm(F_smooth_Z(tmp_z2)))
    Z_N_3 <- append(Z_N_3,qnorm(F_smooth_Z(tmp_z3)))
    
    rho_hat <- cor(Z_N_2,Z_N_3)
    # par_sum <- cbind(par_sum,data.frame(matrix(round(c(a_hat,b_hat,mu_hat,sig_hat,rep(rho_hat,2)),3),nrow=5,ncol=2,byrow=TRUE)))
 par_sum <- data.frame(Z_2,Z_3,Z_N_2,Z_N_3,given)
     }
  return(par_sum)
}

plot_residual_trueab <- function(sims,v=0.9) {
  df <- sims %>% dplyr::select(starts_with("Y"))
  par_sum <- data.frame(matrix(nrow=6,ncol=0))
  
  # Y_not_1_extreme <- df %>% filter(Y_1<quantile(Y_1,v))
  Z_2 <- c()
  Z_3 <- c()
  Z_N_2 <- c()
  Z_N_3 <- c()
  given <- c()
  d <- ncol(df)
  for (j in 1:ncol(df)) {
    Y_given_1_extreme <- df %>% filter(df[,j]>quantile(df[,j],v))
    res <- c(1:d)[-j]
    opt <- list()
    for (i in 2:d) {
      # get initial parameters
     # init_opt <- optim(par=c(1,0,1), fn=Y_likelihood_initial,df=Y_given_1_extreme,given=j,sim=res[i-1],control = list(fnscale=-1))$par
      init_par <- c(0,1)
      # optimise using the initial parameters
      opt[[i-1]] <- optim(par=init_par,fn = Y_likelihood_fix_ab,a=1,b=0,df=Y_given_1_extreme,given=j,sim=res[i-1],control = list(fnscale=-1))
    }
    a_hat <- c(1,1)
    b_hat <- c(0,0)
    mu_hat <- c(opt[[1]]$par[1],opt[[2]]$par[1])
    sig_hat <- c(opt[[1]]$par[2],opt[[2]]$par[2])
    lik <- c(opt[[1]]$value,opt[[2]]$value)
    
    # generate residual Z ----
    Y1 <- Y_given_1_extreme[,j]
    Y2 <- Y_given_1_extreme[,res[1]]
    Y3 <- Y_given_1_extreme[,res[2]]
    
    tmp_z2 <- (Y2-a_hat[1]*Y1)/(Y1^b_hat[1])
    tmp_z3 <- (Y3-a_hat[2]*Y1)/(Y1^b_hat[2])
    
    Z_2 <- append(Z_2,tmp_z2)
    Z_3 <- append(Z_3,tmp_z3)
    given <- append(given,rep(j,(N*(1-v))))
    
    # calculate the normal using the PIT
    Z_N_2 <- append(Z_N_2,qnorm(F_smooth_Z(tmp_z2)))
    Z_N_3 <- append(Z_N_3,qnorm(F_smooth_Z(tmp_z3)))
    
    rho_hat <- cor(Z_N_2,Z_N_3)
     par_sum <- cbind(par_sum,data.frame(matrix(round(c(lik,a_hat,b_hat,mu_hat,sig_hat,rep(rho_hat,2)),3),nrow=6,ncol=2,byrow=TRUE)))
     #par_sum <- data.frame(Z_2,Z_3,Z_N_2,Z_N_3,given)
  }
  return(par_sum)
}

plot_simulated <- function(sims=sims,v=0.99,sim_threshold=0.999,given=1) {
  df <- sims %>% dplyr::select(starts_with("Y"))
  par_sum <- data.frame(matrix(nrow=5,ncol=0))
  # Y_not_1_extreme <- df %>% filter(Y_1<quantile(Y_1,v))
  Z_2 <- c()
  Z_3 <- c()
  Z_N_2 <- c()
  Z_N_3 <- c()
  j <- given
  cond_colours <- c("#C11432","#66A64F","#009ADA")
  d <- ncol(df)
 
    Y_given_1_extreme <- df %>% filter(df[,j]>quantile(df[,j],v))
    res <- c(1:d)[-j]
    opt <- list()
    for (i in 2:d) {
      # get initial parameters
      init_opt <- optim(par=c(1,0,1), fn=Y_likelihood_initial,df=Y_given_1_extreme,given=j,sim=res[i-1],control = list(fnscale=-1))$par
      init_par <- c(init_opt[1],0.2,init_opt[2],init_opt[3])
      # optimise using the initial parameters
      opt[[i-1]] <- optim(par=init_par,fn = Y_likelihood,df=Y_given_1_extreme,given=j,sim=res[i-1],control = list(fnscale=-1))
    }
    a_hat <- c(opt[[1]]$par[1],opt[[2]]$par[1])
    b_hat <- c(opt[[1]]$par[2],opt[[2]]$par[2])
    mu_hat <- c(opt[[1]]$par[3],opt[[2]]$par[3])
    sig_hat <- c(opt[[1]]$par[4],opt[[2]]$par[4])
    
    # generate residual Z ----
    Y_1 <- Y_given_1_extreme[,j]
    Y_2 <- Y_given_1_extreme[,res[1]]
    Y_3 <- Y_given_1_extreme[,res[2]]
    
    tmp_z2 <- (Y_2-a_hat[1]*Y_1)/(Y_1^b_hat[1])
    tmp_z3 <- (Y_3-a_hat[2]*Y_1)/(Y_1^b_hat[2])
    
    Z_2 <- append(Z_2,tmp_z2)
    Z_3 <- append(Z_3,tmp_z3)
    given <- append(given,rep(j,(N*(1-v))))
    
    # calculate the normal using the PIT
    Z_N_2 <- append(Z_N_2,qnorm(F_smooth_Z(tmp_z2)))
    Z_N_3 <- append(Z_N_3,qnorm(F_smooth_Z(tmp_z3)))
    
    rho_hat <- cor(Z_N_2,Z_N_3)
    # generate from normal
    if (j==5) {
      Z_N <- mvrnorm(n=1000,mu=c(0,0),Sigma=matrix(c(1,0,0,1),2,2))
    }
    else {
    Z_N <- mvrnorm(n=1000,mu=c(0,0),Sigma=matrix(c(1,rho_hat,rho_hat,1),2,2))
    }
    Z <- data.frame(Z_2,Z_3)
    
    # transform back to original margins
    Z_star <- norm_to_orig(Z_N=Z_N,emp_res = Z)
    
    U <- runif(1000)
    Y_1_gen <- -log(2*(1-0.999)) + rexp(1000)
    Gen_Y_1 <- data.frame(Y_1=Y_1_gen)
    
    # for each Y, generate a residual and calculate Y_2
    Y_1 <- Gen_Y_1$Y_1
    # Y_2 <- a_hat*Y_1 + Y_1^b_hat *x
    # Y_3 <-  a_hat*Y_1 + Y_1^b_hat *y
    Y_2 <- a_hat*Y_1 + Y_1^b_hat *Z_star[,1]
    Y_3 <-  a_hat*Y_1 + Y_1^b_hat *Z_star[,2]
    
    Gen_Y_1 <- Gen_Y_1 %>% mutate(Y_2=Y_2,Y_3=Y_3) %>% mutate(sim=rep("model",1000))
    names(Gen_Y_1) <- c(paste0("Y",j),paste0("Y",res[1]),paste0("Y",res[2]),"sim")
    # generate Y_1 (extrapolate so above largest observed value)
    
    #plot
    Y1 <- Y_given_1_extreme[,j]
    Y2 <- Y_given_1_extreme[,res[1]]
    Y3 <- Y_given_1_extreme[,res[2]]
    tmp <- data.frame(Y1,Y2,Y3) %>% mutate(sim=rep("data",500))
    names(tmp) <- c(paste0("Y",j),paste0("Y",res[1]),paste0("Y",res[2]),"sim")
    thres <- frechet_laplace_pit( qfrechet(0.999))
    v <- 0.99
    l <- min(Gen_Y_1 %>% dplyr::select(-sim),Y1,Y2,Y3)
    u <- max(Gen_Y_1 %>% dplyr::select(-sim),Y1,Y2,Y3)
    Gen_orig <- rbind(Gen_Y_1,tmp)
    p1 <- ggplot(Gen_orig) +  annotate("rect",xmin=thres, ymin=thres, xmax=Inf,ymax=Inf, alpha=0.25,fill="#C11432") + geom_point(aes(x=Y1,y=Y2,col=sim),alpha=0.5) + geom_vline(xintercept=thres,col=cond_colours[j],linetype="dashed") +geom_abline(slope=0,intercept=thres,col=cond_colours[j],linetype="dashed") +
      scale_color_manual(values = c("data"="black","model" = cond_colours[j])) + xlab(TeX("$Y_1$")) +ylab(TeX("$Y_2$")) +
      xlim(c(l,u)) + ylim(c(l,u)) 
    p2 <- ggplot(Gen_orig)+  annotate("rect",xmin=thres, ymin=thres, xmax=Inf,ymax=Inf, alpha=0.25,fill="#C11432")  + geom_point(aes(x=Y2,y=Y3,col=sim),alpha=0.5) + 
      scale_color_manual(values = c("data"="black","model" = cond_colours[j])) + xlab(TeX("$Y_2$")) +ylab(TeX("$Y_3$")) + geom_vline(xintercept=thres,col=cond_colours[j],linetype="dashed") +geom_abline(slope=0,intercept=thres,col=cond_colours[j],linetype="dashed") +
    xlim(c(l,u)) + ylim(c(l,u))
    p3 <- ggplot(Gen_orig)+  annotate("rect",xmin=thres, ymin=thres, xmax=Inf,ymax=Inf, alpha=0.25,fill="#C11432")  + geom_point(aes(x=Y1,y=Y3,col=sim),alpha=0.5) + 
      scale_color_manual(values = c("data"="black","model" = cond_colours[j]))+ xlab(TeX("$Y_1$")) +ylab(TeX("$Y_3$"))  + geom_vline(xintercept=thres,col=cond_colours[j],linetype="dashed") +geom_abline(slope=0,intercept=thres,col=cond_colours[j],linetype="dashed") +
      xlim(c(l,u)) + ylim(c(l,u))
    p <- grid.arrange(p1,p2,p3,ncol=3)
    
    # par_sum <- cbind(par_sum,data.frame(matrix(round(c(a_hat,b_hat,mu_hat,sig_hat,rep(rho_hat,2)),3),nrow=5,ncol=2,byrow=TRUE)))
    #par_sum <- data.frame(Z_2,Z_3,Z_N_2,Z_N_3,given)

  return(p1)
}

simulated <- function(sims=sims,v=0.9,sim_threshold=0.999,given=1) {
  df <- sims %>% dplyr::select(starts_with("Y"))
  par_sum <- data.frame(matrix(nrow=5,ncol=0))
  
  # Y_not_1_extreme <- df %>% filter(Y_1<quantile(Y_1,v))
  Z_2 <- c()
  Z_3 <- c()
  Z_N_2 <- c()
  Z_N_3 <- c()
  j <- given
  cond_colours <- c("#C11432","#66A64F","#009ADA")
  d <- ncol(df)
  
  Y_given_1_extreme <- df %>% filter(df[,j]>quantile(df[,j],v))
  res <- c(1:d)[-j]
  opt <- list()
  for (i in 2:d) {
    # get initial parameters
    init_opt <- optim(par=c(1,0,1), fn=Y_likelihood_initial,df=Y_given_1_extreme,given=j,sim=res[i-1],control = list(fnscale=-1))$par
    init_par <- c(init_opt[1],0.2,init_opt[2],init_opt[3])
    # optimise using the initial parameters
    opt[[i-1]] <- optim(par=init_par,fn = Y_likelihood,df=Y_given_1_extreme,given=j,sim=res[i-1],control = list(fnscale=-1))
  }
  a_hat <- c(opt[[1]]$par[1],opt[[2]]$par[1])
  b_hat <- c(opt[[1]]$par[2],opt[[2]]$par[2])
  mu_hat <- c(opt[[1]]$par[3],opt[[2]]$par[3])
  sig_hat <- c(opt[[1]]$par[4],opt[[2]]$par[4])
  
  # generate residual Z ----
  Y_1 <- Y_given_1_extreme[,j]
  Y_2 <- Y_given_1_extreme[,res[1]]
  Y_3 <- Y_given_1_extreme[,res[2]]
  
  tmp_z2 <- (Y_2-a_hat[1]*Y_1)/(Y_1^b_hat[1])
  tmp_z3 <- (Y_3-a_hat[2]*Y_1)/(Y_1^b_hat[2])
  
  Z_2 <- append(Z_2,tmp_z2)
  Z_3 <- append(Z_3,tmp_z3)
  given <- append(given,rep(j,500))
  
  # calculate the normal using the PIT
  Z_N_2 <- append(Z_N_2,qnorm(F_smooth_Z(tmp_z2)))
  Z_N_3 <- append(Z_N_3,qnorm(F_smooth_Z(tmp_z3)))
  
  rho_hat <- cor(Z_N_2,Z_N_3)
  # generate from normal
  
  Z_N <- mvrnorm(n=1000,mu=c(0,0),Sigma=matrix(c(1,rho_hat,rho_hat,1),2,2))
  Z <- data.frame(Z_2,Z_3)
  
  # transform back to original margins
  Z_star <- norm_to_orig(Z_N=Z_N,emp_res = Z)
  
  U <- runif(1000)
  Y_1_gen <- -log(2*(1-0.999)) + rexp(1000)
  Gen_Y_1 <- data.frame(Y_1=Y_1_gen)
  
  # for each Y, generate a residual and calculate Y_2
  Y_1 <- Y_given_1_extreme[,j]
  # Y_2 <- a_hat*Y_1 + Y_1^b_hat *x
  # Y_3 <-  a_hat*Y_1 + Y_1^b_hat *y
  Y_2 <- a_hat*Y_1 + Y_1^b_hat *Z_star[,1]
  Y_3 <-  a_hat*Y_1 + Y_1^b_hat *Z_star[,2]
  Gen_Y_1 <- Gen_Y_1 %>% mutate(Y_2=Y_2,Y_3=Y_3) %>% mutate(sim=rep("conditional_model",1000))
  # generate Y_1 (extrapolate so above largest observed value)
  
  #plot
  Y_2 <- Y_given_1_extreme[,res[1]]
  Y_3 <- Y_given_1_extreme[,res[2]]
  thres <- frechet_laplace_pit( qfrechet(0.999))
  
  Gen_orig <- rbind(Gen_Y_1,data.frame(Y_1,Y_2,Y_3) %>% mutate(sim=rep("original_laplace",500)))
  p1 <- ggplot(Gen_orig) + geom_point(aes(x=Y_1,y=Y_2,col=sim),alpha=0.5) + geom_vline(xintercept=thres,col=cond_colours[j],linetype="dashed") +geom_abline(slope=0,intercept=thres,col=cond_colours[j],linetype="dashed") +
    scale_color_manual(values = c("original_laplace"="black","conditional_model" = cond_colours[j])) 
  p2 <- ggplot(Gen_orig) + geom_point(aes(x=Y_2,y=Y_3,col=sim),alpha=0.5) + 
    scale_color_manual(values = c("original_laplace"="black","conditional_model" = cond_colours[j]))  + geom_vline(xintercept=thres,col=cond_colours[j],linetype="dashed") +geom_abline(slope=0,intercept=thres,col=cond_colours[j],linetype="dashed")
  p3 <- ggplot(Gen_orig) + geom_point(aes(x=Y_1,y=Y_3,col=sim),alpha=0.5) + 
    scale_color_manual(values = c("original_laplace"="black","conditional_model" = cond_colours[j]))  + geom_vline(xintercept=thres,col=cond_colours[j],linetype="dashed") +geom_abline(slope=0,intercept=thres,col=cond_colours[j],linetype="dashed")
  p <- grid.arrange(p1,p2,p3,ncol=3)
  # par_sum <- cbind(par_sum,data.frame(matrix(round(c(a_hat,b_hat,mu_hat,sig_hat,rep(rho_hat,2)),3),nrow=5,ncol=2,byrow=TRUE)))
  #par_sum <- data.frame(Z_2,Z_3,Z_N_2,Z_N_3,given)
  
  return(Gen_orig)
}

plot_cond_quantile <- function(sims=sims,v=0.99) {
df <- sims %>% dplyr::select(starts_with("Y"))
par_sum <- data.frame(matrix(nrow=6,ncol=0))
par_sum_init <- data.frame(matrix(nrow=5,ncol=0))
# Y_not_1_extreme <- df %>% filter(Y_1<quantile(Y_1,v))
Z_2 <- c()
Z_3 <- c()

given <- c()
d <- ncol(df)
for (j in 1:ncol(df)) {
  Y_given_1_extreme <- df %>% filter(df[,j]>quantile(df[,j],v))
  res <- c(1:d)[-j]
  opt <- list()
  init_par <- list()
  init_lik <- c()
  for (i in 2:d) {
    # get initial parameters
    init_opt <- optim(par=c(0.5,0,1), fn=Y_likelihood_initial,df=Y_given_1_extreme,given=j,sim=res[i-1],control = list(fnscale=-1))
    init_lik <- append(init_lik,init_opt$value)
    init_par[[i-1]] <- c(init_opt$par[1],0.2,init_opt$par[2],init_opt$par[3])
    # optimise using the initial parameters
    opt[[i-1]] <- optim(par=init_par[[i-1]],fn = Y_likelihood,df=Y_given_1_extreme,given=j,sim=res[i-1],control = list(fnscale=-1))
  }
  a_hat <- c(opt[[1]]$par[1],opt[[2]]$par[1])
  b_hat <- c(opt[[1]]$par[2],opt[[2]]$par[2])
  mu_hat <- c(opt[[1]]$par[3],opt[[2]]$par[3])
  sig_hat <- c(opt[[1]]$par[4],opt[[2]]$par[4])
  lik <-c(opt[[1]]$value,opt[[2]]$value)
  
  a_hat_init <- c(init_par[[1]][1],init_par[[2]][1])
  b_hat_init <- c(0.2,0.2)
  mu_hat_init <- c(init_par[[1]][3],init_par[[2]][3])
  sig_hat_init <- c(init_par[[1]][4],init_par[[2]][4])
  
  # generate residual Z ----
  Y_1 <- Y_given_1_extreme[,j]
  Y_2 <- Y_given_1_extreme[,res[1]]
  Y_3 <- Y_given_1_extreme[,res[2]]
  
  tmp_z2 <- (Y_2-a_hat[1]*Y_1)/(Y_1^b_hat[1])
  tmp_z3 <- (Y_3-a_hat[2]*Y_1)/(Y_1^b_hat[2])
  
  Z_2 <- append(Z_2,tmp_z2)
  Z_3 <- append(Z_3,tmp_z3)
  given <- append(given,rep(j,50))
  

  
  rho_hat <- cor(tmp_zn2,tmp_zn3)
  par_sum <- cbind(par_sum,data.frame(matrix(round(c(lik,a_hat,b_hat,mu_hat,sig_hat,rep(rho_hat,2)),3),nrow=6,ncol=2,byrow=TRUE)))
  par_sum_init <- cbind(par_sum_init,data.frame(matrix(round(c(init_lik,a_hat_init,b_hat_init,mu_hat_init,sig_hat_init),3),nrow=5,ncol=2,byrow=TRUE)))

  
}
p <- ggplot() + geom_point(data=Y_given_1_extreme,aes(x=Y_1,y=Y_2),alpha=0.5) + 
  geom_line(data=data.frame(x=x,y=yl),aes(x=x,y=y),linetype="dashed",col="#C11432") +
  geom_line(data=data.frame(x=x,y=ym),aes(x=x,y=y),linetype="dashed",col="#C11432") +
  geom_line(data=data.frame(x=x,y=yp),aes(x=x,y=y),linetype="dashed",col="#C11432") +
  geom_line(data=data.frame(x=x,y=ylt),aes(x=x,y=y),linetype="dashed",col="black",alpha=0.5) +
  geom_line(data=data.frame(x=x,y=ymt),aes(x=x,y=y),linetype="dashed",col="black",alpha=0.5) +
  geom_line(data=data.frame(x=x,y=ypt),aes(x=x,y=y),linetype="dashed",col="black",alpha=0.5) +
  geom_line(data=data.frame(x=x,y=ylb),aes(x=x,y=y),linetype="dashed",col="#009ada",alpha=0.5) +
  geom_line(data=data.frame(x=x,y=ymb),aes(x=x,y=y),linetype="dashed",col="#009ada",alpha=0.5) +
  geom_line(data=data.frame(x=x,y=ypb),aes(x=x,y=y),linetype="dashed",col="#009ada",alpha=0.5)
return(p)
}

# Laplace margin density
g_laplace <- function(z,a) {
 # -(a-1)*(1/a)*exp(z/a)*(1+(exp(z/a)))^(a-2)
 # -2^(-1/a)*(a-1)*(1/a)*exp(z/a)*(1+(exp(z/a)))^(a-2)
  (a-1)*(-1/a)*exp(-z/a)*(1+(exp(-z/a)))^(a-2)
  
}

G_laplace <- function(z,a) {
(1+exp(-z/a))^(a-1)
}

ggplot() + xlim(-10,20) + geom_function(fun=G_laplace,args=list(a=0.25),col="#C11432")+ geom_function(fun=G_laplace,args=list(a=1/2),col="#009ada")+ geom_function(fun=G_laplace,args=list(a=0.75),col="#66A64F")+ xlab(TeX("$z$")) + ylab(TeX("$G(z)$"))
ggplot() + xlim(-10,20) + geom_function(fun=g_laplace,args=list(a=0.25),col="#C11432")+ geom_function(fun=g_laplace,args=list(a=1/2),col="#009ada")+ geom_function(fun=g_laplace,args=list(a=0.75),col="#66A64F") + ylim(0,1) + xlab(TeX("$z$")) + ylab(TeX("$g(z)$"))
p2 <- ggplot() + xlim(-10,20) + geom_function(fun=g_laplace,args=list(a=1/2)) + ylim(0,1)+ xlab(TeX("$z$")) + ylab(TeX("$g(z)$"))
p3 <- ggplot() + xlim(-10,20) + geom_function(fun=g_laplace,args=list(a=0.75))  + ylim(0,1)+ xlab(TeX("$z$")) + ylab(TeX("$g(z)$"))
grid.arrange(p1,p2,p3,ncol=1)

# function(sims=sims,v=0.99,sim_threshold=0.999,given=1) {
  df <- (sims %>% dplyr::select(starts_with("Y")))[,1:2]
  par_sum <- data.frame(matrix(nrow=5,ncol=0))
  # Y_not_1_extreme <- df %>% filter(Y_1<quantile(Y_1,v))
  Z <- c()
  j <- given
  cond_colours <- c("#C11432","#66A64F","#009ADA")
  d <- ncol(df)
  
  Y_given_1_extreme <- df %>% filter(df[,j]>quantile(df[,j],v))
  res <- c(1:d)[-j]
  opt <- list()
  for (i in 2:d) {
    # get initial parameters
    init_opt <- optim(par=c(1,0,1), fn=Y_likelihood_initial,df=Y_given_1_extreme,given=j,sim=res[i-1],control = list(fnscale=-1))$par
    init_par <- c(init_opt[1],0.2,init_opt[2],init_opt[3])
    # optimise using the initial parameters
    opt[[i-1]] <- optim(par=init_par,fn = Y_likelihood,df=Y_given_1_extreme,given=j,sim=res[i-1],control = list(fnscale=-1))
  }
  a_hat <- c(opt[[1]]$par[1])
  b_hat <- c(opt[[1]]$par[2])
  mu_hat <- c(opt[[1]]$par[3])
  sig_hat <- c(opt[[1]]$par[4])
  
  # generate residual Z ----
  Y_1 <- Y_given_1_extreme[,j]
  Y_2 <- Y_given_1_extreme[,res[1]]

  
  tmp_z2 <- (Y_2-a_hat[1]*Y_1)/(Y_1^b_hat[1])
# fit DL to the residuals
  fg <- function(z) {
    integrate(function(x) exp(-x) * x^(z - 1), 0, Inf)
  }
  DLLL <- function(x,theta) {
    mu <- theta[1]
  sig <- theta[2]
  delta <- theta[3]
  a <- theta[4]
  b <- theta[5]
  x <- (Y_2-a*Y_1)/(Y_1^b)
    if(sig<=0 | delta<=0 | a<(-1) | a>1 | b<0 | b>=1){return(-10e10)}
    -sum(dgnorm(x,mu=mu,alpha=sig,beta=delta,log=T))
  }
  # optimise over the parameters
  opt <- optim(DLLL,x=tmp_z2,par=c(0,1,1,0.8,0.6))
  mu <- opt$par[1]
  sig <- opt$par[2]
  delta <-  opt$par[3]

  # transform back to original margins
  Z_star <- rgnorm(n=1000,mu=mu,alpha=sig,beta=delta)
  
  U <- runif(1000)
  Y_1_gen <- -log(2*(1-0.999)) + rexp(1000)
  Gen_Y_1 <- data.frame(Y_1=Y_1_gen)
  
  # for each Y, generate a residual and calculate Y_2
  Y_1 <- Gen_Y_1$Y_1
   Y_2 <- a_hat*Y_1 + Y_1^b_hat *Z_star
  Gen_Y_1 <- Gen_Y_1 %>% mutate(Y_2=Y_2) %>% mutate(sim=rep("model",1000))
  names(Gen_Y_1) <- c(paste0("Y",j),paste0("Y",res[1]),"sim")
  # generate Y_1 (extrapolate so above largest observed value)
  
  #plot
  Y1 <- Y_given_1_extreme[,j]
  Y2 <- Y_given_1_extreme[,res[1]]
  tmp <- data.frame(Y1,Y2) %>% mutate(sim=rep("data",50))
  names(tmp) <- c(paste0("Y",j),paste0("Y",res[1]),"sim")
  thres <- frechet_laplace_pit( qfrechet(0.999))
  v <- 0.99
  l <- min(Gen_Y_1 %>% dplyr::select(-sim),Y1,Y2)
  u <- max(Gen_Y_1 %>% dplyr::select(-sim),Y1,Y2)
  Gen_orig <- rbind(Gen_Y_1,tmp)
  p1 <- ggplot(Gen_orig) +  annotate("rect",xmin=thres, ymin=thres, xmax=Inf,ymax=Inf, alpha=0.25,fill="#C11432") + geom_point(aes(x=Y1,y=Y2,col=sim),alpha=0.5) + geom_vline(xintercept=thres,col=cond_colours[j],linetype="dashed") +geom_abline(slope=0,intercept=thres,col=cond_colours[j],linetype="dashed") +
    scale_color_manual(values = c("data"="black","model" = cond_colours[j])) + xlab(TeX("$Y_1$")) +ylab(TeX("$Y_2$")) +
    xlim(c(l,u)) + ylim(c(l,u)) 
  p2 <- ggplot(Gen_orig)+  annotate("rect",xmin=thres, ymin=thres, xmax=Inf,ymax=Inf, alpha=0.25,fill="#C11432")  + geom_point(aes(x=Y2,y=Y3,col=sim),alpha=0.5) + 
    scale_color_manual(values = c("data"="black","model" = cond_colours[j])) + xlab(TeX("$Y_2$")) +ylab(TeX("$Y_3$")) + geom_vline(xintercept=thres,col=cond_colours[j],linetype="dashed") +geom_abline(slope=0,intercept=thres,col=cond_colours[j],linetype="dashed") +
    xlim(c(l,u)) + ylim(c(l,u))
  p3 <- ggplot(Gen_orig)+  annotate("rect",xmin=thres, ymin=thres, xmax=Inf,ymax=Inf, alpha=0.25,fill="#C11432")  + geom_point(aes(x=Y1,y=Y3,col=sim),alpha=0.5) + 
    scale_color_manual(values = c("data"="black","model" = cond_colours[j]))+ xlab(TeX("$Y_1$")) +ylab(TeX("$Y_3$"))  + geom_vline(xintercept=thres,col=cond_colours[j],linetype="dashed") +geom_abline(slope=0,intercept=thres,col=cond_colours[j],linetype="dashed") +
    xlim(c(l,u)) + ylim(c(l,u))
  p <- grid.arrange(p1,p2,p3,ncol=3)
  
  # par_sum <- cbind(par_sum,data.frame(matrix(round(c(a_hat,b_hat,mu_hat,sig_hat,rep(rho_hat,2)),3),nrow=5,ncol=2,byrow=TRUE)))
  #par_sum <- data.frame(Z_2,Z_3,Z_N_2,Z_N_3,given)
  

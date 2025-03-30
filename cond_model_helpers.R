# generate from the smoothed empirical distribution of the residuals
F_smooth_Z <- function(Z) {
  Z_smooth <- c()
  for (i in seq_along(Z)) {
    z <- Z[i]
   Z_smooth[i] <-  mean(pnorm((z-Z)/density(Z)$bw))
  }
 return(Z_smooth) 
}

# transform from normal back to kernel smoothed by minimising the square difference
norm_to_orig <- function(ZN,emp_res) {
  Z <- data.frame(matrix(ncol=ncol(ZN),nrow=nrow(ZN)))
  s <- seq(0.05,0.95,length.out=49)
  Zs <- data.frame(matrix(ncol=ncol(ZN),nrow=length(s)))
  # optimise for these 49 values of s
  to_opts <- function(z) {
    return( (mean(pnorm((z-emp_res[,l])/density(emp_res[,l])$bw)) - s[m])^2)
  }
  for (l in 1:ncol(ZN)) {
  for (m in 1:length(s)) {
    Zs[m,l] <- optim(fn=to_opts,par=1,lower = -12, upper = 12, method="Brent")$par
  }
  }

 to_optZ <- function(z) {
 return( (mean(pnorm((z-emp_res[,i])/density(emp_res[,i])$bw)) - pnorm(ZN[j,i]))^2)
 }
 
for (i in 1:ncol(ZN)) {
  for (j in 1:nrow(ZN)) {
    # optimise cdf using 49 linear segments rather than optimising all directly
    if (pnorm(ZN[j,i])< min(s) | pnorm(ZN[j,i])>= max(s) ) {
    Z[j,i] <- optim(fn=to_optZ,par=1,lower = -12,upper = 12,method="Brent")$par
    }
    else {
      k <- which.min(pnorm(ZN[j,i])>s)-1
      a <- (s[2]-s[1])/(Zs[k+1,i]-Zs[k,i])
      b <- -a*Zs[k,i]+s[k]
      Z[j,i] <- (pnorm(ZN[j,i])-b)/(a)
    }
  }
}
  return(Z)
}

# keef constraints for beta
keef_constraint1 <- function(b,a,Y1,Y2) {
  if (b<0 | b>1) {return(10^6)}
  v <- max(Y1)
  ZmAI <- max((Y2-a*Y1)/(Y1^b))
  ZmAD <- max(Y2-Y1)
  return((1-b*ZmAI*v^(b-1)  -a)^2)
}

keef_constraint2 <- function(b,a,Y1,Y2) {
  if (b<0 | b>1) {return(10^6)}
  v <- max(Y1)
  ZmAI <- max((Y2-a*Y1)/(Y1^b))
  ZmAD <- max(Y2-Y1)
  return((1-v^(b-1)*ZmAI+v^(-1)*ZmAD  -a)^2)
}


# generate a table of parameter estimates conditional on (given) each of the specified vector of variables
par_est <- function(df=sims,v=0.99,given=c(1),margin="AGG",method="two_step", a=NULL, keef_constraints=0) {
  lik <- lika <- likb <- lik2 <- a_hat <- b_hat <- mu_hat <- mu_agg_hat <- sig_hat <- sig_agg_hat <- sigl_hat <- sigu_hat <- delta_hat <- deltal_hat <- deltau_hat <- res_var <- c()
  names(df) <- paste0("Y",1:ncol(df))
  d <- ncol(df)
  for (j in given) {
    Y_given1extreme <- df %>% filter(df[,j]>quantile(df[,j],v))
    res <- c(1:d)[-j]
    init_par <- c()
    init_lik <- c()
    if (is.vector(a) & method=="sequential3") {
      a_hat <- a
    }
    for (i in 2:d) {
      # optimise using the initial parameters
      Y1 <- Y_given1extreme[,j]
      Y2 <- Y_given1extreme[,res[i-1]]
      if (method=="sequential") {
        init_para <- c(0.8,0,1)
        opta <- optim(par=init_para,fn = Y_likelihood,df=Y_given1extreme,given=j,sim=res[i-1],b_hat=0,control = list(fnscale=-1,maxit=2000))
        a_hat <- append(a_hat,opta$par[1])
        lika <- append(lika,-opta$value)
        init_parb <- c(0.2,0,1)
        optb <- optim(par=init_parb,fn = Y_likelihood,df=Y_given1extreme,given=j,sim=res[i-1],a_hat=opta$par[1],control = list(fnscale=-1,maxit=2000))
        b_hat <- append(b_hat,optb$par[length(optb$par)-2])
        Z2 <- (Y2-opta$par[1]*Y1)/(Y1^optb$par[1])
        mu_hat <- append(mu_hat,optb$par[length(optb$par)-1])
        sig_hat <- append(sig_hat,optb$par[length(optb$par)])         
        likb <- append(likb,-optb$value)
      }
      if (method=="sequential2") {
        init_para <- c(0.8,0,1)
        opta <- optim(par=init_para,fn = Y_likelihood,df=Y_given1extreme,given=j,sim=res[i-1],b_hat=0,control = list(fnscale=-1,maxit=2000))
        a_hat <- append(a_hat,opta$par[1])
        lika <- append(lika,-opta$value)
        if (1 %in% keef_constraints) {
        b_max1 <- optim(par=0.8,fn = keef_constraint1,a=a_hat[length(a_hat)],Y1=Y1,Y2=Y2,control = list(maxit=2000),lower=0,upper=1, method = "Brent")$par
        } else {b_max1 <- 1}
        if (2 %in% keef_constraints) {
          b_max2 <- optim(par=0.8,fn = keef_constraint2,a=a_hat[length(a_hat)],Y1=Y1,Y2=Y2,control = list(maxit=2000),lower=0, upper=1, method = "Brent")$par
        } else {b_max2 <- 1} 
        b_max <- min(b_max1,b_max2)
        init_parb <- c(b_max/2,0,1)
        #optb <- optim(par=init_parb,fn = Y_likelihood,df=Y_given1extreme,given=j,sim=res[i-1],a_hat=opta$par[1],lower=c(0,-Inf,0),upper = c(b_max,Inf,4),control = list(fnscale=-1,maxit=2000), method = "L-BFGS-B")
        optb <- optim(par=init_parb,fn = Y_likelihood,df=Y_given1extreme,given=j,sim=res[i-1],a_hat=opta$par[1],b_max=b_max,control = list(fnscale=-1,maxit=2000), method = "Nelder-Mead")
        b_hat <- append(b_hat,optb$par[length(optb$par)-2])
        Z2 <- (Y2-opta$par[1]*Y1)/(Y1^optb$par[1])
       optmusig <- optim(par=init_parb,fn = Y_likelihood,df=Y_given1extreme,given=j,sim=res[i-1],a_hat=opta$par[1],b_hat=optb$par[1],control = list(fnscale=-1,maxit=2000))
         mu_hat <- append(mu_hat,optmusig$par[length(optmusig$par)-1])
        sig_hat <- append(sig_hat,optmusig$par[length(optmusig$par)])  
        lik <- append(lik,-optmusig$value)
        likb <- append(likb,-optb$value)
      }
      
      if (method=="sequentialGG") {
        init_para <- c(0.8,0,1,2)
        opta <- optim(par=init_para,fn = Y_likelihoodGG,df=Y_given1extreme,given=j,sim=res[i-1],b_hat=0,control = list(fnscale=-1,maxit=2000))
        a_hat <- append(a_hat,opta$par[1])
        lika <- append(lika,-opta$value)
        if (1 %in% keef_constraints) {
          b_max1 <- optim(par=0.8,fn = keef_constraint1,a=a_hat[length(a_hat)],Y1=Y1,Y2=Y2,control = list(maxit=2000),lower=0,upper=1, method = "Brent")$par
        } else {b_max1 <- 1}
        if (2 %in% keef_constraints) {
          b_max2 <- optim(par=0.8,fn = keef_constraint2,a=a_hat[length(a_hat)],Y1=Y1,Y2=Y2,control = list(maxit=2000),lower=0, upper=1, method = "Brent")$par
        } else {b_max2 <- 1} 
        b_max <- min(b_max1,b_max2)
        init_parb <- c(b_max/2,0,1,2)
        #optb <- optim(par=init_parb,fn = Y_likelihood,df=Y_given1extreme,given=j,sim=res[i-1],a_hat=opta$par[1],lower=c(0,-Inf,0),upper = c(b_max,Inf,4),control = list(fnscale=-1,maxit=2000), method = "L-BFGS-B")
        optb <- optim(par=init_parb,fn = Y_likelihoodGG,df=Y_given1extreme,given=j,sim=res[i-1],a_hat=opta$par[1],b_max=b_max,control = list(fnscale=-1,maxit=2000), method = "Nelder-Mead")
        b_hat <- append(b_hat,optb$par[length(optb$par)-3])
        Z2 <- (Y2-opta$par[1]*Y1)/(Y1^optb$par[1])
        init_parGG <- c(0,1,2)
        opt3 <- optim(par=init_parGG,fn = Y_likelihoodGG,df=Y_given1extreme,given=j,sim=res[i-1],a_hat=opta$par[1],b_hat=optb$par[1],control = list(fnscale=-1,maxit=2000))
        mu_hat <- append(mu_hat,opt3$par[length(opt3$par)-2])
        sig_hat <- append(sig_hat,opt3$par[length(opt3$par)-1]) 
        delta_hat <- append(delta_hat,opt3$par[length(opt3$par)])        
        lik <- append(lik,-opt3$value)
        likb <- append(likb,-optb$value)
      }
      if (method=="sequential3") {
        init_parb <- c(0.2,0,1)
        optb <- optim(par=init_parb,fn = Y_likelihood,df=Y_given1extreme,given=j,sim=res[i-1],a_hat=a_hat[i-1],control = list(fnscale=-1,maxit=2000))
        b_hat <- append(b_hat,optb$par[length(optb$par)-2])
        Z2 <- (Y2-a_hat[i-1]*Y1)/(Y1^optb$par[1])
        optmusig <- optim(par=init_parb,fn = Y_likelihood,df=Y_given1extreme,given=j,sim=res[i-1],a_hat=a_hat[i-1],b_hat=optb$par[2],control = list(fnscale=-1,maxit=2000))
        mu_hat <- append(mu_hat,optmusig$par[length(optmusig$par)-1])
        sig_hat <- append(sig_hat,optmusig$par[length(optmusig$par)])         
        likb <- append(likb,optb$value)
      }      
      if (method=="two_step" | (method=="one_step" & margin=="Normal")) {
        init_par <- c(0.8,0.2,0,1)
        opt <- optim(par=init_par,fn = Y_likelihood,df=Y_given1extreme,given=j,sim=res[i-1],control = list(fnscale=-1,maxit=2000))
        a_hat <- append(a_hat,opt$par[1])
        b_hat <- append(b_hat,opt$par[2])
        mu_hat <- append(mu_hat,opt$par[3])
        sig_hat <- append(sig_hat,opt$par[4])
        lik <- append(lik,opt$value)
        }
      if (margin=="AGGdelta" & method=="one_step") {
        opt <- optim(fn=AGGdelta_onestep,x=data.frame(Y1,Y2),par=c(0,1,1.5,1.5,0.8,0.3),control=list(maxit=2000))
        a_hat <- append(a_hat,opt$par[5])
        b_hat <- append(b_hat,opt$par[6])
        mu_agg_hat <- append(mu_agg_hat,opt$par[1])
        sig_agg_hat <- append(sig_agg_hat,opt$par[2])
        deltal_hat <- append(deltal_hat,opt$par[3])
        deltau_hat <- append(deltau_hat,opt$par[4])
        lik <- append(lik,opt$value)
      }
      
      if (margin=="GenGaus" & method!="one_step") {
        opt <- optim(fn=NLL_GenGaus,x=Z2,par=c(mean(Z2),sd(Z2),1.5),control=list(maxit=2000),method = "SANN")
        mu_agg_hat <- append(mu_agg_hat,opt$par[1])
        sig_agg_hat <- append(sig_agg_hat,opt$par[2])
        delta_hat <- append(delta_hat,opt$par[3])
        lik2 <- append(lik2,opt$value)
      }
      
      if (margin=="AGGdelta" & method!="one_step") {
        opt <- optim(fn=NLL_AGGdelta,x=Z2,par=c(mean(Z2),sd(Z2),1.2,1.8),control=list(maxit=2000),method = "SANN")
        mu_agg_hat <- append(mu_agg_hat,opt$par[1])
        sig_agg_hat <- append(sig_agg_hat,opt$par[2])
        deltal_hat <- append(deltal_hat,opt$par[3])
        deltau_hat <- append(deltau_hat,opt$par[4])
        lik2 <- append(lik2,opt$value)
      }
      
      if (margin=="AGGsig" & method!="one_step") {
        opt <- optim(fn=NLL_AGGsig,x=Z2,par=c(mean(Z2),sd(Z2),sd(Z2),1.5),control=list(maxit=2000),method = "SANN")
        mu_agg_hat <- append(mu_agg_hat,opt$par[1])
        sigl_hat <- append(sigl_hat,opt$par[2])
        sigu_hat <- append(sigu_hat,opt$par[3])
        delta_hat <- append(delta_hat,opt$par[4])
        lik2 <- append(lik2,opt$value)
      }
      if (margin=="AGG" & method!="one_step") {
        opt <- optim(fn=NLL_AGG,x=Z2,par=c(mean(Z2),sd(Z2),sd(Z2),1.2,1.8),control=list(maxit=2000),method = "SANN")
        mu_agg_hat <- append(mu_agg_hat,opt$par[1])
        sigl_hat <- append(sigl_hat,opt$par[2])
        sigu_hat <- append(sigu_hat,opt$par[3])
        deltal_hat <- append(deltal_hat,opt$par[4])
        deltau_hat <- append(deltau_hat,opt$par[5])
        lik2 <- append(lik2,opt$value)
      }
      res_var <- append(res_var,res[i-1])
    }
  }
  nas <- rep(NA,length(a_hat)) # NA values for parameters not used by a given method
  if (margin=="Normal") {
    if (method=="one_step" | method=="two_step") {
  par_sum <- data.frame("lik" = lik,"lika" = nas,"likb"=nas,"lik2"=nas,
                        "a" = a_hat, "b" = b_hat,
                        "mu" = mu_hat,"mu_agg"=nas,
                        "sig" = sig_hat,"sig_agg"=nas,"sigl"=nas,"sigu"=nas,
                        "delta"=nas,"deltal"=nas,"deltau"=nas,
                        "given" = rep(given,each=(d-1)), "res" = res_var)}
    if (method %in% c("sequential","sequential2")) {
  par_sum <- data.frame("lik"=nas, "lika" = lika,"likb"=likb,"lik2"=lik,
                        "a" = a_hat, "b" = b_hat,
                        "mu" = mu_hat,"mu_agg"=nas,
                        "sig" = sig_hat,"sig_agg"=nas,"sigl"=nas,"sigu"=nas,
                        "delta"=nas,"deltal"=nas,"deltau"=nas,
                        "given" = rep(given,each=(d-1)), "res" = res_var)  }
  }
  if (margin=="AGGdelta" & method=="one_step") {
    par_sum <- data.frame("lik" = lik,"lika"=nas,"likb"=nas,"lik2"=nas,
                          "a" = a_hat, "b" = b_hat,
                          "mu" = nas,"mu_agg"= mu_agg_hat,
                          "sig" = nas,"sig_agg"= sig_agg_hat,"sigl"=nas,"sigu"=nas,
                          "delta"=nas,"deltal" = deltal_hat, "deltau" = deltau_hat,
                          "given" = rep(given,each=(d-1)), "res" = res_var)
  }
  if (margin=="AGGdelta" & method=="two_step") {
    par_sum <- data.frame("lik" = lik, "lika"=nas,"likb"=nas,"lik2"=lik2,
                          "a" = a_hat, "b" = b_hat,
                          "mu" = mu_hat,"mu_agg"=mu_agg_hat,
                          "sig" = sig_hat,"sig_agg"=sig_agg_hat,"sigl"=nas,"sigu"=nas,
                          "delta"=nas,"deltal" = deltal_hat, "deltau" = deltau_hat,
                          "given" = rep(given,each=(d-1)), "res" = res_var)
  }
  
  if (margin=="GenGaus" & method=="two_step") {
    par_sum <- data.frame("lik" = lik, "lika"=nas,"likb"=nas,"lik2"=lik2,
                          "a" = a_hat, "b" = b_hat,
                          "mu" = mu_hat,"mu_agg"=mu_agg_hat,
                          "sig" = sig_hat,"sig_agg" = sig_agg_hat,"sigl" = nas,"sigu" = nas,
                          "delta"= delta_hat,"deltal" = nas, "deltau" = nas,
                          "given" = rep(given,each=(d-1)), "res" = res_var)
  }
  
  if (margin=="AGGsig" & method=="two_step") {
    par_sum <- data.frame("lik" = lik, "lika"=nas,"likb"=nas,"lik2"=lik2,
                          "a" = a_hat, "b" = b_hat,
                          "mu" = mu_hat,"mu_agg"=mu_agg_hat,
                          "sig" = sig_hat,"sig_agg"=nas,"sigl"=sigl_hat,"sigu"=sigu_hat,
                          "delta"= delta_hat,"deltal" = nas, "deltau" = nas,
                          "given" = rep(given,each=(d-1)), "res" = res_var)
  }
  
  if (margin=="AGG" & method=="two_step") {
    par_sum <- data.frame("lik" = lik, "lika"=nas,"likb"=nas,"lik2"=lik2,
                          "a" = a_hat, "b" = b_hat,
                          "mu" = mu_hat,"mu_agg"=mu_agg_hat,
                          "sig" = sig_hat,"sig_agg"=nas,"sigl"=sigl_hat,"sigu"=sigu_hat,
                          "delta"=nas,"deltal" = deltal_hat, "deltau" = deltau_hat,
                          "given" = rep(given,each=(d-1)), "res" = res_var)
  }
  
  if (margin=="AGG" & method=="sequentialGG") {
    par_sum <- data.frame("lik" = lik, "lika"=lika,"likb"=likb,"lik2"=lik2,
                          "a" = a_hat, "b" = b_hat,
                          "mu" = mu_hat,"mu_agg"=mu_agg_hat,
                          "sig" = sig_hat,"sig_agg"=nas,"sigl"=sigl_hat,"sigu"=sigu_hat,
                          "delta"=delta_hat,"deltal" = deltal_hat, "deltau" = deltau_hat,
                          "given" = rep(given,each=(d-1)), "res" = res_var)
  }
  
  if (margin=="AGG" & method %in% c("sequential","sequential2")) {
    par_sum <- data.frame("lik" = nas, "lika"= lika,"likb"= likb,"lik2"=lik2,
                          "a" = a_hat, "b" = b_hat,
                          "mu" = mu_hat,"mu_agg"=mu_agg_hat,
                          "sig" = sig_hat,"sig_agg"=nas,"sigl"=sigl_hat,"sigu"=sigu_hat,
                          "delta"=nas,"deltal" = deltal_hat, "deltau" = deltau_hat,
                          "given" = rep(given,each=(d-1)), "res" = res_var)
  }
  if (margin=="AGGsig" & method %in% c("sequential","sequential2")) {
    par_sum <- data.frame("lik" = nas, "lika"= lika,"likb"= likb,"lik2"=lik2,
                          "a" = a_hat, "b" = b_hat,
                          "mu" = mu_hat,"mu_agg"=mu_agg_hat,
                          "sig" = sig_hat,"sig_agg"=nas,"sigl"=sigl_hat,"sigu"=sigu_hat,
                          "delta"= delta_hat,"deltal" = nas, "deltau" = nas,
                          "given" = rep(given,each=(d-1)), "res" = res_var)
  }
  if (margin=="AGGdelta" & method %in% c("sequential","sequential2")) {
    par_sum <- data.frame("lik" = nas, "lika"= lika,"likb"= likb,"lik2"=lik2,
                          "a" = a_hat, "b" = b_hat,
                          "mu" = mu_hat,"mu_agg"=mu_agg_hat,
                          "sig" = sig_hat,"sig_agg"=sig_agg_hat,"sigl"=nas,"sigu"=nas,
                          "delta"= nas,"deltal" = deltal_hat, "deltau" = deltau_hat,
                          "given" = rep(given,each=(d-1)), "res" = res_var)
  }
  if (margin=="GenGaus" & method %in% c("sequential","sequential2")) {
    par_sum <- data.frame("lik" = lik, "lika"= lika,"likb"= likb,"lik2"=lik2,
                          "a" = a_hat, "b" = b_hat,
                          "mu" = mu_hat,"mu_agg"=mu_agg_hat,
                          "sig" = sig_hat,"sig_agg"=sig_agg_hat,"sigl"= nas,"sigu"=nas,
                          "delta"= delta_hat,"deltal" = nas, "deltau" = nas,
                          "given" = rep(given,each=(d-1)), "res" = res_var)
  }
  
  if (margin=="AGG" & method %in% c("sequential3")) {
    par_sum <- data.frame("lik" = nas, "lika"= nas,"likb"= likb,"lik2"=lik2,
                          "a" = a_hat, "b" = b_hat,
                          "mu" = mu_hat,"mu_agg"=mu_agg_hat,
                          "sig" = sig_hat,"sig_agg"=nas,"sigl"=sigl_hat,"sigu"=sigu_hat,
                          "delta"=nas,"deltal" = deltal_hat, "deltau" = deltau_hat,
                          "given" = rep(given,each=(d-1)), "res" = res_var)
  }
  
  if (margin=="AGGdelta" & method %in% c("sequential","sequential2")) {
    par_sum <- data.frame("lik" = nas, "lika"=lika,"likb"=likb,"lik2"=lik2,
                          "a" = a_hat, "b" = b_hat,
                          "mu" = mu_hat,"mu_agg"=mu_agg_hat,
                          "sig" = sig_hat,"sig_agg"=sig_agg_hat,"sigl"=nas,"sigu"=nas,
                          "delta"=nas,"deltal" = deltal_hat, "deltau" = deltau_hat,
                          "given" = rep(given,each=(d-1)), "res" = res_var)
  }
  
  return(par_sum)
}

par_est_ite <- function(df=sims,d1j = d1j, v=0.9, given=c(1),N=100, show_ite=FALSE,mu_init=NULL,sig_init=NULL,b_init=NULL,method="onephi",SN=NULL, b_inc=FALSE)  {
  names(df) <- paste0("Y",1:ncol(df))
  d <- ncol(df)
  Y_given1extreme <- df %>% filter(df[,given]>quantile(df[,given],v))
  nv <- nrow(Y_given1extreme)
  res <- c(1:d)[-given]
  d1j <- d1j/1000000
  a <- b <- mu <- sig <- data.frame(matrix(ncol=(N+1),nrow = (d-1)))
  if (method=="onephi") {
  phi. <- c()
  }
  if (method=="twophi") {
    SN <- as.numeric(SN[-given])*as.numeric(SN[given])+as.numeric(!SN[-given])*as.numeric(!SN[given])
    phi1. <- c()
    phi0. <- c()
  }
  # calculate a with initial values for mu and sigma
  if (is.numeric(mu_init) & is.numeric(sig_init)) {
   mu[,1] <- mu_init
   sig[,1] <- sig_init
  } else {
  mu[,1] <- 0
  sig[,1] <- 1

  }
  if (is.numeric(b_init)) {
    b[,1] <- b_init
  } else {
    b[,1] <- 0
  }
  if (method=="onephi") {
    phi_init <- 1
    opt <- optim(fn=NLL_expalpha_HT,par=phi_init,df = Y_given1extreme,d1j. = d1j,mu1=as.numeric(mu[,1]),sig1=as.numeric(sig[,1]),beta1=as.numeric(b[,1]),d.=d,given.=given,res.=res,control=list(maxit=2000),method = "BFGS")
    phi <- opt$par
    phi. <- append(phi.,phi)
    a[,1] <- exp(-phi*d1j)
  }
  if (method=="twophi") {
    phi_init <- c(1,1)
    opt <- optim(fn=NLL_expalpha_twophi,par=phi_init,df=Y_given1extreme,d1j.=d1j,SN.=SN,mu1=as.numeric(mu[,1]),sig1=as.numeric(sig[,1]),beta1=as.numeric(b[,1]),d.=d,given.=given,res.=res,control=list(maxit=2000),method = "BFGS")
    phi1 <- opt$par[1]
    phi0 <- opt$par[2]
    phi1. <- append(phi1.,phi1)
    phi0. <- append(phi0.,phi0)
    a[,1] <- exp(-(phi1*as.numeric(SN)+phi0*as.numeric(!SN))*d1j)
    }
  for (i in 1:N) {
    for (j in 1:(d-1)) {
      if (b_inc==FALSE) {
   mu[j,i+1] <- 1/nv*sum(as.numeric(Y_given1extreme[,res[j]])-as.numeric(a[j,i])*as.numeric(Y_given1extreme[,given]))
   sig[j,i+1] <- sqrt(1/nv*sum((as.numeric(Y_given1extreme[,res[j]])-as.numeric(a[j,i])*as.numeric(Y_given1extreme[,given])-as.numeric(mu[j,i+1]))^2))
   b[j,i+1] <- 0
      }
      if (b_inc==TRUE) {
        init_parb <- c(b[j,i],mu[j,i],sig[j,i])
        optb <- optim(par=init_parb,fn = Y_likelihood,df=Y_given1extreme,given=given,sim=res[j],a_hat=as.numeric(a[j,i]),control = list(fnscale=-1,maxit=2000))
        b[j,i+1] <- optb$par[1]
        mu[j,i+1] <- optb$par[2]
        sig[j,i+1] <- optb$par[3]
      }
   }
   # calculate a 
    if (method=="onephi") {
    opt <- optim(fn=NLL_expalpha_HT,df = Y_given1extreme, d1j. = d1j, mu1=as.numeric(mu[,i+1]),sig1=as.numeric(sig[,i+1]),beta1=as.numeric(b[,i+1]),d.=d,given.=given,res.=res,par=phi_init,control=list(maxit=2000),method = "BFGS")
    phi <- opt$par
    phi. <- append(phi.,phi)
    a[,i+1] <- exp(-phi*d1j)
    }
    if (method=="twophi") {
      opt <- optim(fn=NLL_expalpha_twophi,par=phi_init,df=Y_given1extreme,d1j.=d1j,SN.=SN,mu1=as.numeric(mu[,i+1]),sig1=as.numeric(sig[,i+1]),beta1=as.numeric(b[,i+1]),d.=d,given.=given,res.=res,control=list(maxit=2000),method = "BFGS")
      phi1 <- opt$par[1]
      phi0 <- opt$par[2]
      phi1. <- append(phi1.,phi1)
      phi0. <- append(phi0.,phi0)
      a[,i+1] <- exp(-(phi1*as.numeric(SN)+phi0*as.numeric(!SN))*d1j)
    }
  }
    par_sum <- data.frame("a" = as.numeric(a[,N+1]),"mu" = as.numeric(mu[,N+1]), "sig" = as.numeric(sig[,N+1]))
  if (show_ite == TRUE) {
    if (method=="onephi") {
    return(list(a,b,mu,sig,par_sum,phi.))
    }
    if (method=="twophi") {
      return(list(a,b,mu,sig,par_sum,phi1.,phi0.))
    }
  } else {return(par_sum)}
}

# calculate the observed residuals from estimated alpha and beta
observed_residuals <- function(df=sims,given=1,v=0.99,a=NULL,b=NULL) {
  j <- given
  a_hat <- b_hat <- res_var <- c()
  tmp_z <- tmp_z1 <- c()
  df_orig <- df
  names(df) <- paste0("Y",1:ncol(df))
  d <- ncol(df)
  Y_given1extreme <- df %>% filter(df[,j]>quantile(df[,j],v))
  nv <- nrow(Y_given1extreme)
  res <- c(1:d)[-j]
  init_par <- c()
  for (i in 2:d) {
    if (is.numeric(a) & is.numeric(b)) {
      a_hat <- a[i-1]
      b_hat <- b[i-1]
      res_var <- append(res_var,rep(paste0("Z",res[i-1]),nv))
      Y1 <- Y_given1extreme[,j]
      Y2 <- Y_given1extreme[,res[i-1]]
      tmp_z <- append(tmp_z,(Y2-a_hat*Y1/(Y1^b_hat)))
    }
  }
  Z <- data.frame(res_var,tmp_z) %>% mutate(res_var=factor(res_var,levels=paste0("Z",res))) %>% group_by(res_var) %>% 
    mutate(row = row_number()) %>%
    tidyr::pivot_wider(names_from = res_var, values_from = tmp_z) %>% 
    dplyr::select(-row)
  if (sum(names(df_orig)!=paste0("Y",1:ncol(df_orig)))==ncol(df_orig)) {
    names(Z) <- names(df_orig)[-j]
  }
  return(Z)
}


plot_simulated <- function(sims=sims,v=0.99,sim_threshold=0.999,given=1) {
  df <- sims %>% dplyr::select(starts_with("Y"))
  par_sum <- data.frame(matrix(nrow=5,ncol=0))
  # Y_not_1_extreme <- df %>% filter(Y_1<quantile(Y_1,v))
  Z2 <- c()
  Z3 <- c()
  ZN2 <- c()
  ZN3 <- c()
  j <- given
  cond_colours <- c("#C11432","#66A64F","#009ADA")
  d <- ncol(df)
 
    Y_given1extreme <- df %>% filter(df[,j]>quantile(df[,j],v))
    res <- c(1:d)[-j]
    opt <- list()
    for (i in 2:d) {
      # get initial parameters
      init_opt <- optim(par=c(1,0,1), fn=Y_likelihood_initial,df=Y_given1extreme,given=j,sim=res[i-1],control = list(fnscale=-1))$par
      init_par <- c(init_opt[1],0.2,init_opt[2],init_opt[3])
      # optimise using the initial parameters
      opt[[i-1]] <- optim(par=init_par,fn = Y_likelihood,df=Y_given1extreme,given=j,sim=res[i-1],control = list(fnscale=-1))
    }
    a_hat <- c(opt[[1]]$par[1],opt[[2]]$par[1])
    b_hat <- c(opt[[1]]$par[2],opt[[2]]$par[2])
    mu_hat <- c(opt[[1]]$par[3],opt[[2]]$par[3])
    sig_hat <- c(opt[[1]]$par[4],opt[[2]]$par[4])
    
    # generate residual Z ----
    Y1 <- Y_given1extreme[,j]
    Y2 <- Y_given1extreme[,res[1]]
    Y3 <- Y_given1extreme[,res[2]]
    
    tmp_z2 <- (Y2-a_hat[1]*Y1)/(Y1^b_hat[1])
    tmp_z3 <- (Y3-a_hat[2]*Y1)/(Y1^b_hat[2])
    
    Z2 <- append(Z2,tmp_z2)
    Z3 <- append(Z3,tmp_z3)
    given <- append(given,rep(j,(N*(1-v))))
    
    # calculate the normal using the PIT
    ZN2 <- append(ZN2,qnorm(F_smooth_Z(tmp_z2)))
    ZN3 <- append(ZN3,qnorm(F_smooth_Z(tmp_z3)))
    
    rho_hat <- cor(ZN2,ZN3)
    # generate from normal
    if (j==5) {
      ZN <- mvrnorm(n=1000,mu=c(0,0),Sigma=matrix(c(1,0,0,1),2,2))
    }
    else {
    ZN <- mvrnorm(n=1000,mu=c(0,0),Sigma=matrix(c(1,rho_hat,rho_hat,1),2,2))
    }
    Z <- data.frame(Z2,Z3)
    
    # transform back to original margins
    Z_star <- norm_to_orig(ZN=ZN,emp_res = Z)
    
    U <- runif(1000)
    Y1gen <- -log(2*(1-0.999)) + rexp(1000)
    GenY1 <- data.frame(Y1=Y1gen)
    
    # for each Y, generate a residual and calculate Y2
    Y1 <- GenY1$Y1
    # Y2 <- a_hat*Y1 + Y1^b_hat *x
    # Y3 <-  a_hat*Y1 + Y1^b_hat *y
    Y2 <- a_hat*Y1 + Y1^b_hat *Z_star[,1]
    Y3 <-  a_hat*Y1 + Y1^b_hat *Z_star[,2]
    
    GenY1 <- GenY1 %>% mutate(Y2=Y2,Y3=Y3) %>% mutate(sim=rep("model",1000))
    names(GenY1) <- c(paste0("Y",j),paste0("Y",res[1]),paste0("Y",res[2]),"sim")
    # generate Y1 (extrapolate so above largest observed value)
    
    #plot
    Y1 <- Ygiven1extreme[,j]
    Y2 <- Ygiven1extreme[,res[1]]
    Y3 <- Ygiven1extreme[,res[2]]
    tmp <- data.frame(Y1,Y2,Y3) %>% mutate(sim=rep("data",500))
    names(tmp) <- c(paste0("Y",j),paste0("Y",res[1]),paste0("Y",res[2]),"sim")
    thres <- frechet_laplace_pit( qfrechet(0.999))
    v <- 0.99
    l <- min(GenY1 %>% dplyr::select(-sim),Y1,Y2,Y3)
    u <- max(GenY1 %>% dplyr::select(-sim),Y1,Y2,Y3)
    Gen_orig <- rbind(GenY1,tmp)
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
  Z2 <- c()
  Z3 <- c()
  ZN2 <- c()
  ZN3 <- c()
  j <- given
  cond_colours <- c("#C11432","#66A64F","#009ADA")
  d <- ncol(df)
  
  Y_given1extreme <- df %>% filter(df[,j]>quantile(df[,j],v))
  res <- c(1:d)[-j]
  opt <- list()
  for (i in 2:d) {
    # get initial parameters
    init_opt <- optim(par=c(1,0,1), fn=Y_likelihood_initial,df=Y_given1extreme,given=j,sim=res[i-1],control = list(fnscale=-1))$par
    init_par <- c(init_opt[1],0.2,init_opt[2],init_opt[3])
    # optimise using the initial parameters
    opt[[i-1]] <- optim(par=init_par,fn = Y_likelihood,df=Y_given1extreme,given=j,sim=res[i-1],control = list(fnscale=-1))
  }
  a_hat <- c(opt[[1]]$par[1],opt[[2]]$par[1])
  b_hat <- c(opt[[1]]$par[2],opt[[2]]$par[2])
  mu_hat <- c(opt[[1]]$par[3],opt[[2]]$par[3])
  sig_hat <- c(opt[[1]]$par[4],opt[[2]]$par[4])
  
  # generate residual Z ----
  Y1 <- Y_given1extreme[,j]
  Y2 <- Y_given1extreme[,res[1]]
  Y3 <- Y_given1extreme[,res[2]]
  
  tmp_z2 <- (Y2-a_hat[1]*Y1)/(Y1^b_hat[1])
  tmp_z3 <- (Y_3-a_hat[2]*Y1)/(Y1^b_hat[2])
  
  Z2 <- append(Z2,tmp_z2)
  Z3 <- append(Z3,tmp_z3)
  given <- append(given,rep(j,500))
  
  # calculate the normal using the PIT
  ZN2 <- append(ZN2,qnorm(F_smooth_Z(tmp_z2)))
  ZN3 <- append(ZN3,qnorm(F_smooth_Z(tmp_z3)))
  
  rho_hat <- cor(ZN2,ZN3)
  # generate from normal
  
  ZN <- mvrnorm(n=1000,mu=c(0,0),Sigma=matrix(c(1,rho_hat,rho_hat,1),2,2))
  Z <- data.frame(Z2,Z3)
  
  # transform back to original margins
  Z_star <- norm_to_orig(ZN=ZN,emp_res = Z)
  
  U <- runif(1000)
  Y1_gen <- -log(2*(1-0.999)) + rexp(1000)
  Gen_Y1 <- data.frame(Y1=Y1_gen)
  
  # for each Y, generate a residual and calculate Y_2
  Y1 <- Y_given1extreme[,j]
  # Y_2 <- a_hat*Y_1 + Y_1^b_hat *x
  # Y_3 <-  a_hat*Y_1 + Y_1^b_hat *y
  Y2 <- a_hat*Y1 + Y1^b_hat *Z_star[,1]
  Y3 <-  a_hat*Y1 + Y1^b_hat *Z_star[,2]
  Gen_Y1 <- Gen_Y1 %>% mutate(Y2=Y2,Y3=Y3) %>% mutate(sim=rep("conditional_model",1000))
  # generate Y_1 (extrapolate so above largest observed value)
  
  #plot
  Y2 <- Y_given1extreme[,res[1]]
  Y3 <- Y_given1extreme[,res[2]]
  thres <- frechet_laplace_pit( qfrechet(0.999))
  
  Gen_orig <- rbind(Gen_Y1,data.frame(Y1,Y2,Y3) %>% mutate(sim=rep("original_laplace",500)))
  p1 <- ggplot(Gen_orig) + geom_point(aes(x=Y1,y=Y2,col=sim),alpha=0.5) + geom_vline(xintercept=thres,col=cond_colours[j],linetype="dashed") +geom_abline(slope=0,intercept=thres,col=cond_colours[j],linetype="dashed") +
    scale_color_manual(values = c("original_laplace"="black","conditional_model" = cond_colours[j])) 
  p2 <- ggplot(Gen_orig) + geom_point(aes(x=Y2,y=Y3,col=sim),alpha=0.5) + 
    scale_color_manual(values = c("original_laplace"="black","conditional_model" = cond_colours[j]))  + geom_vline(xintercept=thres,col=cond_colours[j],linetype="dashed") +geom_abline(slope=0,intercept=thres,col=cond_colours[j],linetype="dashed")
  p3 <- ggplot(Gen_orig) + geom_point(aes(x=Y1,y=Y3,col=sim),alpha=0.5) + 
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
# Y_not_1_extreme <- df %>% filter(Y1<quantile(Y1,v))
Z2 <- c()
Z3 <- c()

given <- c()
d <- ncol(df)
for (j in 1:ncol(df)) {
  Y_given1extreme <- df %>% filter(df[,j]>quantile(df[,j],v))
  res <- c(1:d)[-j]
  opt <- list()
  init_par <- list()
  init_lik <- c()
  for (i in 2:d) {
    # get initial parameters
    init_opt <- optim(par=c(0.5,0,1), fn=Y_likelihood_initial,df=Y_given1extreme,given=j,sim=res[i-1],control = list(fnscale=-1))
    init_lik <- append(init_lik,init_opt$value)
    init_par[[i-1]] <- c(init_opt$par[1],0.2,init_opt$par[2],init_opt$par[3])
    # optimise using the initial parameters
    opt[[i-1]] <- optim(par=init_par[[i-1]],fn = Y_likelihood,df=Y_given1extreme,given=j,sim=res[i-1],control = list(fnscale=-1))
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
  Y1 <- Y_given1extreme[,j]
  Y2 <- Y_given1extreme[,res[1]]
  Y3 <- Y_given1extreme[,res[2]]
  
  tmp_z2 <- (Y2-a_hat[1]*Y1)/(Y1^b_hat[1])
  tmp_z3 <- (Y3-a_hat[2]*Y1)/(Y1^b_hat[2])
  
  Z2 <- append(Z2,tmp_z2)
  Z3 <- append(Z3,tmp_z3)
  given <- append(given,rep(j,50))
  
  rho_hat <- cor(tmp_zn2,tmp_zn3)
  par_sum <- cbind(par_sum,data.frame(matrix(round(c(lik,a_hat,b_hat,mu_hat,sig_hat,rep(rho_hat,2)),3),nrow=6,ncol=2,byrow=TRUE)))
  par_sum_init <- cbind(par_sum_init,data.frame(matrix(round(c(init_lik,a_hat_init,b_hat_init,mu_hat_init,sig_hat_init),3),nrow=5,ncol=2,byrow=TRUE)))
}
p <- ggplot() + geom_point(data=Y_given1extreme,aes(x=Y1,y=Y2),alpha=0.5) + 
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
  (a-1)*(-1/a)*exp(-z/a)*(1+(exp(-z/a)))^(a-2)
}

G_laplace <- function(z,a) {
(1+exp(-z/a))^(a-1)
}
  

# find site closest to a vector of lon lat coordinates
find_site_index <- function(site = Inverness,grid_uk = uk_sf_rot %>% dplyr::select()) {
  # convert site coordinates
  site_sf <- st_sfc(st_point(site),crs=4326) 
  if (nrow(grid_uk)==445) {
  site_sf <- st_transform(site_sf, crs=27700) # set BNG crs
  }
  x <- which.min(as.numeric(st_distance(site_sf,grid_uk)))
  # check it works
  #return(tm_shape(uk_sf_rot) + tm_dots() + tm_shape(site_sf) + tm_dots(col="#C11432") + tm_shape(uk_sf_rot[x,])+ tm_dots(col="#009ADA"))
  return(x)
}

# function for plotting parameter estimates on a map and against distance
map_param <- function(tmp_est,method = "AGG", facet_var = "cond_site",title_map="",grid_uk=uk_temp_sf) {
  misscol <- "aquamarine"
  Nsites <- max(tmp_est$res, tmp_est$given,na.rm=TRUE)
  if (identical(facet_var,"cond_site")) {
    Nfacet <- length(unique(tmp_est$cond_site))
    facet_label <- unique(tmp_est$cond_site)
    nrow_facet <- ceiling(length(unique(tmp_est$cond_site))/3)
    legend_outside_size <- 0.3
  } else if (identical(facet_var,"tau")) {
    Nfacet <- nrow(tmp_est)/Nsites
    facet_label <- levels(tmp_est$tau)
    nrow_facet <- 1
    legend_outside_size <- 0.1
  } else if (identical(facet_var,"q")) {
    Nfacet <- nrow(tmp_est)/Nsites
    facet_label <- levels(tmp_est$q)
    nrow_facet <- 1
    legend_outside_size <- 0.2 
  } else if (identical(facet_var[2], c("tau"))) {
    Nfacet <- length(unique(tmp_est$cond_site))*length(unique(tmp_est$tau))
    facet_label <- list(levels(tmp_est %>% dplyr::select(facet_var[1]) %>% pull()),levels(tmp_est %>% dplyr::select(facet_var[2]) %>% pull()))   
     nrow_facet <- length(unique(tmp_est$cond_site))
    legend_outside_size <- 0.2 
  }
  rep_Nsites <- nrow(tmp_est)/Nsites
  tmp <- tmp_est %>% mutate(site_index = rep(1:Nsites,rep_Nsites))
  if (nrow(grid_uk)==445) {
    uk_temp_sf <- grid_uk %>% mutate(site_index=1:Nsites) %>% dplyr::select(site_index) %>% cbind(ukcp18[,1:8])
    uk_tmp <- tmp %>% left_join(uk_temp_sf,by=c("site_index")) 
  } else {
    uk_temp_sf <- grid_uk %>% mutate(site_index=1:Nsites)
    uk_tmp <- tmp %>% left_join(uk_temp_sf,by=c("site_index"))
  }
  uk_tmp1 <- st_as_sf(uk_tmp)
  if (method=="q") {
    if (identical(facet_var,"q")) {
    pq <- tm_shape(uk_tmp1) + tm_dots(col="value",style="cont",size=0.3,palette="viridis",colorNA=misscol,title=TeX("$q_p$"), textNA = "Conditioning site") + tm_facets(by=facet_var,nrow = nrow_facet) +  tm_layout(panel.labels = facet_label,legend.outside.size=legend_outside_size,asp=0.5,legend.text.size = 1,legend.title.size=1.5, title=title_map) 
    } else if (identical(facet_var, c("q","tau"))) {
      pq <- tm_shape(uk_tmp1) + tm_dots(col="value",style="cont",size=0.3,palette="Blues",colorNA=misscol,title=TeX("$q_p$"), textNA = "Conditioning site") + tm_facets(by=facet_var) +  tm_layout(legend.outside.size=legend_outside_size,asp=0.5,legend.text.size = 1,legend.title.size=1.5, title=title_map) 
    } else if (identical(facet_var, c("rl","tau"))) {
      pq <- tm_shape(uk_tmp1) + tm_dots(col="rl",style="cont",size=0.3,palette="Blues",colorNA=misscol,title="Return level (days)", textNA = "Conditioning site") + tm_facets(by=c("q","tau")) +  tm_layout(legend.outside.size=legend_outside_size,asp=0.5,legend.text.size = 1,legend.title.size=1.5, title=title_map) 
    } else if (identical(facet_var, c("cond_site","tau"))) {
      pq <- tm_shape(uk_tmp1) + tm_dots(col="value",n=8,style="quantile",size=0.3,palette="Blues",colorNA=misscol,title=TeX("$q_p$"), textNA = "Conditioning site") + tm_facets(by=facet_var) +  tm_layout(legend.outside.size=legend_outside_size,asp=0.5,legend.text.size = 1,legend.title.size=1.5, title=title_map) 
    }
    
  }
  if (method=="rl") {
    if (identical(facet_var, c("cond_site","tau"))) {
      pq <- tm_shape(uk_tmp1) + tm_dots(col="rl",n=8,style="quantile",size=0.3,palette="Blues",colorNA=misscol,title="Return level (days)", textNA = "Conditioning site") + tm_facets(by=facet_var) +  tm_layout(legend.outside.size=legend_outside_size,asp=0.5,legend.text.size = 1,legend.title.size=1.5, title=title_map) 
    }
  }
  # return level difference rl075-rl025
  if (method=="rldiff") {
    if (identical(facet_var, c("cond_site","tau"))) {
      pq <- tm_shape(uk_tmp1) + tm_dots(col="rldiff",n=8,style="quantile",size=0.3,palette="Blues",colorNA=misscol,title="Return level difference (days)", textNA = "Conditioning site") + tm_facets(by=facet_var) +  tm_layout(legend.outside.size=legend_outside_size,asp=0.5,legend.text.size = 1,legend.title.size=1.5, title=title_map) 
    }
  }
  # relative return level difference (rl075-rl025)/rl075
  if (method=="rlreldiff") {
    if (identical(facet_var, c("cond_site","tau"))) {
      pq <- tm_shape(uk_tmp1) + tm_dots(col="rlreldiff",n=8,style="quantile",size=0.3,palette="Blues",colorNA=misscol,title="Return level \n relative difference (days)", textNA = "Conditioning site") + tm_facets(by=facet_var) +  tm_layout(legend.outside.size=legend_outside_size,asp=0.5,legend.text.size = 1,legend.title.size=1.5, title=title_map) 
    }
  }
  if (method=="AGG") {
  uk_tmp1 <- uk_tmp1 %>% mutate(sigdiff=sigu-sigl) %>% mutate(deltadiff=deltau-deltal)
  }
  

  if (method=="max_tau") {
    lims <- seq(min(as.numeric(uk_tmp1$a),na.rm = TRUE),max(as.numeric(uk_tmp1$a),na.rm = TRUE),length.out=5)
    pa <- tm_shape(uk_tmp1) + tm_dots(col="a",style="cont",size=0.3,palette="viridis",colorNA=misscol,title=TeX("$\\alpha_{max}$"), breaks=lims,textNA = "Conditioning site") + tm_facets(by=facet_var,nrow = nrow_facet) +  tm_layout(panel.labels = facet_label,legend.outside.size=legend_outside_size,asp=0.5,legend.text.size = 1,legend.title.size=1.5, title=title_map) 
    pa_tau <- tm_shape(uk_tmp1) + tm_dots(col="a_tau",style="cat",size=0.3,palette="-RdBu",colorNA=misscol,title=TeX("$\\tau_{\\alpha_{max}}$"), breaks=lims,textNA = "Conditioning site") + tm_facets(by=facet_var,nrow = nrow_facet) +  tm_layout(panel.labels = facet_label,legend.outside.size=legend_outside_size,asp=0.5,legend.text.size = 1,legend.title.size=1.5, title=title_map) 
  }
  if (method %in% c("Normal","AGG")) {
    lims <- seq(min(as.numeric(uk_tmp1$a),na.rm = TRUE),max(as.numeric(uk_tmp1$a),na.rm = TRUE),length.out=5)
    pa <- tm_shape(uk_tmp1) + tm_dots(col="a",style="cont",size=0.3,palette="viridis",colorNA=misscol,title=TeX("$\\alpha$"), breaks=lims,textNA = "Conditioning site") + tm_facets(by=facet_var,nrow = nrow_facet) +  tm_layout(panel.labels = facet_label,legend.outside.size=legend_outside_size,asp=0.5,legend.text.size = 1,legend.title.size=1.5, title=title_map) 

    lims <- seq(min(as.numeric(uk_tmp1$b),na.rm = TRUE),max(as.numeric(uk_tmp1$b),na.rm = TRUE),length.out=6)
    pb <- tm_shape(uk_tmp1) + tm_dots(col="b",style="cont",size=0.3,palette="viridis",colorNA=misscol,title=TeX("$\\beta$"), breaks=lims,textNA = "Conditioning site") + tm_facets(by=facet_var,nrow = nrow_facet) +  tm_layout(panel.labels = facet_label,legend.outside.size=legend_outside_size,asp=0.5,legend.text.size = 1,legend.title.size=1.5, title=title_map) 
    lims <- seq(min(as.numeric(uk_tmp1$mu),na.rm = TRUE),max(as.numeric(uk_tmp1$mu),na.rm = TRUE),length.out=6)
    pmu <- tm_shape(uk_tmp1) + tm_dots(col="mu",style="cont",size=0.3,palette="-RdBu",colorNA=misscol,midpoint=0,title=TeX("$\\mu$"), breaks=lims,textNA = "Conditioning site") + tm_facets(by=facet_var,nrow = nrow_facet) +  tm_layout(panel.labels = facet_label,legend.outside.size=legend_outside_size,asp=0.5,legend.text.size = 1,legend.title.size=1.5, title=title_map) 
    
    sigmax <- 3
    lims <- seq(min(as.numeric(uk_tmp1$sig),na.rm = TRUE),min(max(as.numeric(uk_tmp1$sig),na.rm = TRUE),sigmax),length.out=6)
    psig <- tm_shape(uk_tmp1 %>% filter(is.na(sig) |sig < sigmax)) + tm_dots(col="sig",style="cont",size=0.3,palette="viridis",colorNA=misscol,title=TeX("$\\sigma$"), breaks=lims,textNA = "Conditioning site") + tm_facets(by=facet_var,nrow = nrow_facet) +  tm_layout(panel.labels = facet_label,legend.outside.size=legend_outside_size,asp=0.5,legend.text.size = 1,legend.title.size=1.5, title=title_map) 
    
      }

    if (method=="AGG") {
    lims <- seq(min(as.numeric(uk_tmp1$mu_agg),na.rm = TRUE),max(as.numeric(uk_tmp1$mu_agg),na.rm = TRUE),length.out=6)
    pmuagg <- tm_shape(uk_tmp1) + tm_dots(col="mu_agg",style="cont",size=0.3,palette="-RdBu",colorNA=misscol,midpoint=0,title=TeX("$\\mu_{AGG}$"), breaks=lims,textNA = "Conditioning site") + tm_facets(by=facet_var,nrow = nrow_facet) +  tm_layout(panel.labels = facet_label,legend.outside.size=legend_outside_size,asp=0.5,legend.text.size = 1,legend.title.size=1.5, title=title_map) 

    lims <- seq(min(as.numeric(uk_tmp1$sigl),na.rm = TRUE),sigmax,length.out=6)
    psigl <- tm_shape(uk_tmp1 %>% filter(is.na(sigl) | sigl < sigmax)) + tm_dots(col="sigl",style="cont",size=0.3,palette="viridis",colorNA=misscol,title=TeX("$\\sigma_l$"), breaks=lims,textNA = "Conditioning site") + tm_facets(by=facet_var,nrow = nrow_facet) +  tm_layout(panel.labels = facet_label,legend.outside.size=legend_outside_size,asp=0.5,legend.text.size = 1,legend.title.size=1.5, title=title_map) 

    lims <- seq(min(as.numeric(uk_tmp1$sigu),na.rm = TRUE),sigmax,length.out=6)
    psigu <- tm_shape(uk_tmp1 %>% filter(is.na(sigu) | sigu < sigmax)) + tm_dots(col="sigu",style="cont",size=0.3,palette="viridis",colorNA=misscol,title=TeX("$\\sigma_u$"), breaks=lims,textNA = "Conditioning site") + tm_facets(by=facet_var,nrow = nrow_facet) +  tm_layout(panel.labels = facet_label,legend.outside.size=legend_outside_size,asp=0.5,legend.text.size = 1,legend.title.size=1.5, title=title_map) 

    lims <- seq(min(as.numeric(uk_tmp1$sigdiff),na.rm = TRUE),sigmax,length.out=6)
    psigdiff <- tm_shape(uk_tmp1 %>% filter(is.na(sigdiff) | sigdiff < sigmax)) + tm_dots(col="sigdiff",style="cont",size=0.3,palette="-RdBu",midpoint=0,colorNA=misscol,title=TeX("$\\sigma_u-\\sigma_l$"), breaks=lims,textNA = "Conditioning site") + tm_facets(by=facet_var,nrow = nrow_facet) +  tm_layout(panel.labels = facet_label,legend.outside.size=legend_outside_size,asp=0.5,legend.text.size = 1,legend.title.size=1.5, title=title_map) 

    deltamax <- 5
    deltamin <- min(as.numeric(uk_tmp1$deltal),as.numeric(uk_tmp1$deltau),na.rm = TRUE)
    lims <- seq(deltamin,deltamax,length.out=6)
    pdeltal <- tm_shape(uk_tmp1 %>% filter(is.na(deltal)| deltal < deltamax)) + tm_dots(col="deltal",style="cont",size=0.3,palette="viridis",colorNA=misscol,title=TeX("$\\delta_l$"), breaks=lims,textNA = "Conditioning site") + tm_facets(by=facet_var,nrow = nrow_facet) +  tm_layout(panel.labels = facet_label,legend.outside.size=legend_outside_size,asp=0.5,legend.text.size = 1,legend.title.size=1.5, title=title_map) 

    lims <- seq(deltamin,deltamax,length.out=6)
    pdeltau <- tm_shape(uk_tmp1 %>% filter(is.na(deltau) | deltau < deltamax)) + tm_dots(col="deltau",style="cont",size=0.3,palette="viridis",colorNA=misscol,title=TeX("$\\delta_u$"), breaks=lims,textNA = "Conditioning site") + tm_facets(by=facet_var,nrow = nrow_facet) +  tm_layout(panel.labels = facet_label,legend.outside.size=legend_outside_size,asp=0.5,legend.text.size = 1,legend.title.size=1.5, title=title_map) 
    
    deltamin <- -5
    lims <- seq(deltamin,deltamax,length.out=6)
    pdeltadiff <- tm_shape(uk_tmp1 %>% filter(is.na(deltadiff ) | (deltadiff < deltamax&deltadiff>deltamin))) + tm_dots(col="deltadiff",style="cont",size=0.3,palette="-RdBu",colorNA=misscol,midpoint=0,title=TeX("$\\delta_u-\\delta_l$"), breaks=lims,textNA = "Conditioning site") + tm_facets(by=facet_var,nrow = nrow_facet) +  tm_layout(panel.labels = facet_label,legend.outside.size=legend_outside_size,asp=0.5,legend.text.size = 1,legend.title.size=1.5, title=title_map) 
    return(list(pa,pb,pmu,psig,pmuagg,psigl,psigu,psigdiff,pdeltal,pdeltau,pdeltadiff)) 
    } else if (method=="Normal") {
    return(list(pa,pb,pmu,psig))
    } else if (method=="max_tau") {return(list(pa,pa_tau))}
  else if (method %in% c("q","rl","rldiff","rlreldiff")) {return(pq)}
  else {return(list(pa,pmu,psig))}
}

shift_time <- function(sims=sims,cond_site=cond_site,tau=0, Ndays_season = 90) {
  if (tau==0) {
    return(sims)
  }
  dayshift <- c(-3:3)
  Nyears <- nrow(sims)/Ndays_season
  dayremove <- c(rep(Ndays_season,3),0,1,2,3)-c(2,1,rep(0,5))
  daysremove_condsite <- daysremove_othersites <-  c()
  daysremove_condsite <- as.numeric(sapply(1:Nyears,function(i){append(daysremove_condsite,dayremove[dayshift %in% (0:-tau)[-1]]+Ndays_season*(i-1))}))
  daysremove_othersites <- as.numeric(sapply(1:Nyears,function(i){append(daysremove_othersites,dayremove[dayshift %in% (0:tau)[-1]]+Ndays_season*(i-1))}))
  
  sims_tau <-  data.frame(matrix(ncol=ncol(sims),nrow=nrow(sims)-abs(tau)*Nyears))
  names(sims_tau) <- names(sims)
  sims_tau[,cond_site] <- sims[-daysremove_condsite,cond_site]
  sims_tau[,-cond_site] <- sims[-daysremove_othersites,-cond_site]
  return(sims_tau)
}

#' Link parameter estimates to spatial locations
#'
#' @param tmp_est Output of par_est function or rbind multiple outputs.
#' @param grid_uk sf object of spatial location points.
#'
#' @return sf object (to be used in map_param function)
#' @export
#'
#' @examples
est_join_spatial <- function(tmp_est,grid_uk=uk_temp_sf) {
Nsites <- max(tmp_est$res, tmp_est$given,na.rm=TRUE)
rep_Nsites <- nrow(tmp_est)/Nsites
tmp <- tmp_est %>% mutate(site_index = rep(1:Nsites,rep_Nsites))
if (nrow(grid_uk)==445) {
  uk_temp_sf <- grid_uk %>% mutate(site_index=1:Nsites) %>% dplyr::select(site_index) %>% cbind(ukcp18[,1:8])
  uk_tmp <- tmp %>% left_join(uk_temp_sf,by=c("site_index")) 
} else {
  uk_temp_sf <- grid_uk %>% mutate(site_index=1:Nsites)
  uk_tmp <- tmp %>% left_join(uk_temp_sf,by=c("site_index"))
}
uk_tmp1 <- st_as_sf(uk_tmp)
return(uk_tmp1)
}

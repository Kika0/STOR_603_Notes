# check Simon's script
# load the data
load("../luna/kristina/MSdata01/ukgd_cpm85_5k_x84y20_MSdata01.RData",verb=TRUE)

load("../luna/kristina/P2q/ukgd_cpm85_5k_x84y20_MSp2q.RData",verb=TRUE)

load("../luna/kristina/MSGpdParam/ukgd_cpm85_5k_x84y20.MSGpdParam.2025-02-26.RData",verb=TRUE)

o.u        <- data01$u
o.q        <- double(length(o.u))
ms.thgpd.u <- 0.94
## assuming here gpdpar is gpdpar$scale gpdpar$shape gpdpar$thresh each of length(data01) ie gpdpar[length(data01),3]

for (i in seq_along(o.u)) {
  
  if(!is.na(o.u[i])) {
    
    if(o.u[i]<=ms.thgpd.u) {
      
      o.q[i] <- qgam.p2q.fn[[i]](o.u[i])
      
    } else {
      
      m <- (1-ms.thgpd.u) / (1-o.u[i])
      
      if(gpdpar[i,2]!=0) {
        o.q[i] <- gpdpar[i,3] + gpdpar[i,1]*(m^(gpdpar[i,2]) - 1)/gpdpar[i,2] # coles p81
      } else  o.q[i] <- gpdpar[i,3] + gpdpar[i,1]*log(m)                            # coles p81
      
    }
    
  }
  
}

plot(data01$x,o.q)
dev.copy(jpeg,filename="../check_P2q_x84y20.jpg");
dev.off ();

# write a function to transforma field on a given date
unif_orig_P2q <- function(data,P2q,gpdpar,ms.thgpd.u=0.94) {
  y <- rep(0,length(data))
  for (i in seq_along(data)) {
    
    if(!is.na(data[i])) {
      
      if(y[i]<=ms.thgpd.u) {
        
        y[i] <- P2q[[i]](o.u[i])
        
      } else {
        
        m <- (1-ms.thgpd.u) / (1-y[i])
        
        if(gpdpar[i,2]!=0) {
          y[i] <- gpdpar[i,3] + gpdpar[i,1]*(m^(gpdpar[i,2]) - 1)/gpdpar[i,2] # coles p81
        } else  y[i] <- gpdpar[i,3] + gpdpar[i,1]*log(m)                            # coles p81
        
      }
      
    }
    
  }
  return(y)
}

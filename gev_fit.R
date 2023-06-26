gev_fit <- function(x,maxit=10000) {
  
    # ismev package has option to predefine
    # initial parameter values
    sig_init <- sqrt(6 * var(x))/pi
    mu_init <- mean(x) - 0.57722 * sig_init
    xi_init <- 0.1
    init <- c(mu_init,sig_init,xi_init)
    gev_lik <- function(par) {
    mu <- par[1]
    sig <- par[2]
    xi <- par[3]
    -(-length(x)*log(sig)-(1+1/xi)*sum(log(1+xi*((x-mu)/sig)))-sum((1+xi*((x-mu)/sig))^(-1/xi)))
    }
    opt <- optim(par=init,fn = gev_lik,hessian = TRUE,control = list(maxit=maxit))
    z <- list()
    z$trans <- FALSE
    z$nllh <- opt$value
    z$data <- x
    z$mle <- opt$par
    z$cov <- solve(opt$hessian)
    z$se <- sqrt(diag(z$cov))
    z$hessian <- opt$hessian
    return(z)
    }
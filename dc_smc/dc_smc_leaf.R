dc_smc_leaf <- function(data, N, p=1, Ntotal) {
  if (length(data)%%p != 0) stop("p should be a divisor of the length of data")
  nt <- length(data)/p
  data_bar <- mean(data)
  
  min_lim <- log(.Machine$double.xmin)
  
  # Storage
  mu <- array(rep(NA, N*Ntotal*nt), dim = c(N, Ntotal, nt))
  mu_re <- array(rep(NA, N*nt), dim = c(N, nt))
  
  w_log <- array(rep(NA, N*nt), dim = c(N, nt))
  W <- array(rep(NA, N*nt), dim = c(N, nt))
  
  
  # Initialization
  mu_0 <- rep(0, N)
  w_log_0 <- rep(0, N)
  W_0 <- rep(1/N, N)
  
  
  
  for (i in 1:(nt-1)) {
    if (i == 1) {
      tim0 <- proc.time()
      
      U <- runif(1, 0, 1)
      A <- Sys_resamp(W=W_0, P=N, U=U)
      mu_0_re <- mu_0[A]
      
      
      ## MCMC kernel
      for (j in 1:(Ntotal-1)) {
        if (j == 1) {
          mu_star <- mvrnorm(1, mu = mu_0_re, Sigma = 10000*diag(N))
          mu[,j,i] <- mu_star
        }
        
        mu_star <- mvrnorm(1, mu = mu[,j,i], Sigma = 10000*diag(N))
        mu[,j+1,i] <- mu_star
      }
      
      ## Update weights
      w_log[,i] <- weight_check(data = data[1:p],
                                mu = mu[,Ntotal,i],
                                min_lim = min_lim)
      
      if (max(w_log[,i]) == min_lim) {
        W[,i] <- rep(1/N, N)
      } else {
        W[,i] <- exp(w_log[,i]) / sum(exp(w_log[,i]))
      }
      
      
      tim1 <- proc.time()
      timT <- tim1 - tim0
      print(paste0('iteration ', i, ', CPU time is ', round((tim1-tim0)[[3]], 2), 's, ', 'total CPU time is ', round(timT[[3]], digits=2), 's'))
    }
    
    
    tim0 <- proc.time()
    
    U <- runif(1, 0, 1)
    A <- Sys_resamp(W=W[,i], P=N, U=U)
    mu_re[,i] <- mu[A,Ntotal,i]
    
    
    ## MCMC kernel
    for (j in 1:(Ntotal-1)) {
      if (j == 1) {
        mu_star <- mvrnorm(1, mu = mu_re[,i], Sigma = 10000*diag(N))
        mu[,j,i+1] <- mu_star
      }
      
      mu_star <- mvrnorm(1, mu = mu[,j,i+1], Sigma = 10000*diag(N))
      mu[,j+1,i+1] <- mu_star
    }
    
    ## Update weights
    w_log[,i+1] <- weight_check(data = data[(i*p+1):((i+1)*p)],
                                 mu = mu[,Ntotal,i+1],
                                 min_lim = min_lim)
    if (max(w_log[,i+1]) == min_lim) {
      W[,i+1] <- rep(1/N, N)
    } else {
      W[,i+1] <- exp(w_log[,i+1]) / sum(exp(w_log[,i+1]))
    }
    
    
    tim1 <- proc.time()
    timT <- timT + tim1 - tim0
    print(paste0('iteration ', i+1, ', CPU time is ', round((tim1-tim0)[[3]], 2), 's, ', 'total CPU time is ', round(timT[[3]], digits=2), 's'))
  }
  
  U <- runif(1, 0, 1)
  A <- Sys_resamp(W=W[,nt], P=N, U=U)
  
  return(list(mu_post=mu, W=W, w_log=w_log, A=A))
}









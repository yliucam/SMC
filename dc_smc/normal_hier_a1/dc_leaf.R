dc_leaf <- function(y, n, prior_mu, prior_sigma, N, Ntotal) {
  
  # Storage
  x <- array(rep(NA, N*n), dim = c(N, n))
  theta_0 <- rep(NA, N)
  w_log <- array(rep(NA, N*n), dim = c(N, n))
  W <- array(rep(NA, N*n), dim = c(N, n))
  w_theta0_log <- rep(NA, N)
  W_theta0 <- rep(NA, N)
  
  # Initialization
  x_0 <- 0
  w_log_0 <- 0
  W_0 <- 1/N
  theta_0_init <- 1
  
  for (i in 1:n) {
    
    #mu_post <- (prior_mu + prior_sigma^2 * y[i]) / (1+prior_sigma^2)
    #sigma_post <- sqrt(prior_sigma^2 * (1 + prior_sigma^2))
    #x[,i] <- rnorm(N, mu_post, sigma_post)
    x[,i] <- rnorm(N, y[i], 1)
    
    ## Weights update
    w_log[,i] <- sapply(x[,i], function(x) dnorm(y[i], x, 1, log = T))
    #w_log[,i] <- weight_check(w_log = w_log[,i], min_lim = min_lim)
    
    #if (sum(exp(w_log[,i])) == Inf) w_log[,i] <- w_log[,i] - max(w_log[,i])
    #W[,i] <- exp(w_log[,i]) / sum(exp(w_log[,i]))
    W[,i] <- exp(w_log[,i] - matrixStats::logSumExp(w_log[,i]))
    
    #U <- runif(1, 0, 1)
    #A <- Sys_resamp(W = W[,i], P = N, U = U)
    #x[,i] <- x[A, i]
    #w_log[,i] <- w_log[A, i]
    #W[,i] <- exp(w_log[,i] - matrixStats::logSumExp(w_log[,i]))
  }
  
  
  return(list(x=x, w_log=w_log, W=W, theta_0=theta_0, w_theta0_log=w_theta0_log, W_theta0=W_theta0))
}

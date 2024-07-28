dc_leaf <- function(y, n, prior_mu, prior_sigma, N, min_lim) {
  
  # Storage
  x <- array(rep(NA, N*n), dim = c(N, n))
  w_log <- array(rep(NA, N*n), dim = c(N, n))
  W <- array(rep(NA, N*n), dim = c(N, n))
  
  # Initialization
  x_0 <- 0
  w_log_0 <- 0
  W_0 <- 1/N
  
  for (i in 1:n) {
    mu_post <- prior_sigma^2 * y[i] / (1+prior_sigma^2)
    sigma_post <- sqrt(prior_sigma^2 * (1+prior_sigma^2))
    x[,i] <- rnorm(N, mu_post, sigma_post)
    
    ## Weights update
    w_log[,i] <- sapply(x[,i], function(x) dnorm(y[i], x, 1, log = T))
    w_log[,i] <- weight_check(w_log = w_log[,i], min_lim = min_lim)
    
    if (sum(exp(w_log[,i])) == Inf) w_log[,i] <- w_log[,i] - max(w_log[,i])
    W[,i] <- exp(w_log[,i]) / sum(exp(w_log[,i]))
  }
  
  
  return(list(x=x, w_log=w_log, W=W))
}

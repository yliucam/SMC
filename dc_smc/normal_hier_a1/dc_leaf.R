dc_leaf <- function(y, n, prior_mu, prior_sigma, N, min_lim, Ntotal) {
  
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
  
  for (i in 1:(n-1)) {
    # x[,1] & theta_0
    if (i == 1) {
      for (j in 1:(Ntotal-1)) {
        if (j == 1) {
          ## Update x[,1]
          x_star <- rnorm(N, x_0, sqrt(2.38))
          
          prior_log <- dnorm(x_0, 0, theta_0_init, log = T)
          prior_star_log <- sapply(x_star, function(x) dnorm(x, 0, theta_0_init, log = T))
          
          like_log <- dnorm(y[i], x_0, 1, log = T)
          like_star_log <- sapply(x_star, function(x) dnorm(y[i], x, 1, log = T))
          
          post_log <- prior_log + like_log
          post_star_log <- prior_star_log + like_star_log
          
          ratio <- post_star_log - post_log
          a <- sapply(ratio, function(x) min(0, x))
          index <- which(log(runif(N, 0, 1)) < a)
          x[index, i] <- x_star[index]
          x[-index, i] <- x_0
          
          
          ## Update theta_0
          theta_0_star <- rlnorm(N, log(theta_0_init), log(2.38))
          
          prior_log <- dexp(theta_0_init, 1, log = T)
          prior_star_log <- sapply(theta_0_star, function(x) dexp(x, 1, log = T))
          
          like_log <- dnorm(cbind(x[,i], theta_0_init), function(x) dnorm(x[1], 0, x[2], log = T))
          like_star_log <- dnorm(cbind(x[,i], theta_0_star), function(x) dnorm(x[1], 0, x[2], log = T))
          
          post_log <- prior_log + like_log
          post_star_log <- prior_star_log + like_star_log
          
          prop_current_to_star <- sapply(theta_0_star, function(x) dlnorm(x, log(theta_0_init), log(2.38), log = T))
          prop_star_to_current <- sapply(theta_0_star, function(x) dlnorm(theta_0_init, log(x), log(2.38), log = T))
          
          ratio <- post_star_log - post_log + prop_star_to_current - prop_current_to_star
          a <- sapply(ratio, function(x) min(0, x))
          index <- which(log(runif(N, 0, 1)) < a)
          theta_0[index] <- theta_0_star
          theta_0[-index] <- theta_0_init
        }
        
        
        ## Update x[,1]
        x_star <- sapply(x[,i], rnorm(1, x, sqrt(2.38)))
        
        prior_log <- sapply(cbind(x[,i], theta_0), function(x) dnorm(x[1], 0, x[2], log = T))
        prior_star_log <- sapply(cbind(x_star, theta_0), function(x) dnorm(x[1], 0, x[2], log = T))
        
        like_log <- sapply(x[,i], function(x) dnorm(y[i], x, 1, log = T))
        like_star_log <- sapply(x_star, function(x) dnorm(y[i], x, 1, log = T))
        
        post_log <- prior_log + like_log
        post_star_log <- prior_star_log + like_star_log
        
        ratio <- post_star_log - post_log
        a <- sapply(ratio, function(x) min(0, x))
        index <- which(log(runif(N, 0, 1)) < a)
        x[index, i] <- x_star[index]
        
        
        ## Update theta_0
        theta_0_star <- sapply(theta_0, function(x) rlnorm(1, log(x), log(2.38)))
        
        prior_log <- sapply(theta_0, function(x) dexp(x, 1, log = T))
        prior_star_log <- sapply(theta_0_star, function(x) dexp(x, 1, log = T))
        
        like_log <- sapply(cbind(x[,i], theta_0), function(x) dnorm(x[1], 0, x[2], log = T))
        like_star_log <- sapply(cbind(x[,i], theta_0_star), function(x) dnorm(x[1], 0, x[2], log = T))
        
        post_log <- prior_log + like_log
        post_star_log <- prior_star_log + like_star_log
        
        prop_current_to_star <- sapply(cbind(theta_0, theta_0_star), function(x) dlnorm(x[2], log(x[1]), log(2.38), log = T))
        prop_star_to_current <- sapply(cbind(theta_0, theta_0_star), function(x) dlnorm(x[1], log(x[2]), log(2.38), log = T))
        
        ratio <- post_log - post_star_log + prop_star_to_current - prop_current_to_star
        a <- sapply(ratio, function(x) min(0, x))
        index <- which(log(runif(N, 0, 1)) < a)
        theta_0[index] <- theta_0_star
      }
      
      # Weights update
      w_log[,i] <- sapply(x[,i], function(x) dnorm(y[i], x, 1, log = T))
      w_log[,i] <- weight_check(w_log = w_log[,i], min_lim = min_lim)
      if (sum(exp(w_log[,i])) == Inf) w_log[,i] <- w_log[,i] - max(w_log[,i])
      W[,i] <- exp(w_log[,i]) / sum(exp(w_log[,i]))
      
      w_theta0_log <- sapply(cbind(x[,i], theta_0), function(x) dnorm(x[1], 0, x[2], log = T))
      w_theta0_log <- weight_check(w_log = w_theta0_log, min_lim = min_lim)
      if (sum(exp(w_theta0_log)) == Inf) w_theta0_log <- w_theta0_log - max(w_theta0_log)
      W_theta0 <- exp(w_theta0_log) / sum(exp(w_theta0_log))
    }
    
    
    mu_post <- prior_sigma^2 * y[i+1] / (1+prior_sigma^2)
    sigma_post <- sqrt(prior_sigma^2 * (1+prior_sigma^2))
    x[,i+1] <- rnorm(N, mu_post, sigma_post)
    
    ## Weights update
    w_log[,i+1] <- sapply(x[,i+1], function(x) dnorm(y[i+1], x, 1, log = T))
    w_log[,i+1] <- weight_check(w_log = w_log[,i+1], min_lim = min_lim)
    
    if (sum(exp(w_log[,i+1])) == Inf) w_log[,i+1] <- w_log[,i+1] - max(w_log[,i+1])
    W[,i+1] <- exp(w_log[,i+1]) / sum(exp(w_log[,i+1]))
  }
  
  
  return(list(x=x, w_log=w_log, W=W, theta_0=theta_0, w_theta0_log=w_theta0_log, W_theta0=W_theta0))
}

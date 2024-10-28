mcmc_kernel <- function(data, x_init, theta_init, prior_lambda, N, Ntotal) {
  
  # Storage
  x <- rep(NA, N)
  theta <- rep(NA, N)
  
  for (i in 1:(Ntotal - 1)) {
    if (i == 1) {
      # Update X
      x_star <- sapply(x_init, function(x) rnorm(1, x, 2.38))
      
      prior_log <- apply(cbind(x_init, theta_init), 1, function(x) dnorm(x[1], 0, x[2], log = T))
      prior_star_log <- apply(cbind(x_star, theta_init), 1, function(x) dnorm(x[1], 0, x[2], log = T))
      
      like_log <- sapply(x_init, function(x) dnorm(data, x, 1, log = T))
      like_star_log <- sapply(x_star, function(x) dnorm(data, x, 1, log = T))
      
      post_log <- prior_log + like_log
      post_star_log <- prior_star_log + like_star_log
      
      ratio <- post_star_log - post_log
      alpha <- sapply(ratio, function(x) min(0, x))
      index <- which(log(runif(N, 0, 1)) < alpha)
      x[index] <- x_star[index]
      x[-index] <- x_init[-index]
      
      # Update theta
      theta_star <- sapply(theta_init, function(x) rlnorm(1, log(x), log(2.38)))
      
      prior_log <- sapply(theta_init, function(x) dexp(x, prior_lambda, log = T))
      prior_star_log <- sapply(theta_star, function(x) dexp(x, prior_lambda, log = T))
      
      like_log <- apply(cbind(x, theta_init), 1, function(x) dnorm(x[1], 0, x[2], log = T))
      like_star_log <- apply(cbind(x, theta_star), 1, function(x) dnorm(x[1], 0, x[2], log = T))
      
      post_log <- prior_log + like_log
      post_star_log <- prior_star_log + like_star_log
      
      prop_current_to_star <- apply(cbind(theta_star, theta_init), 1, function(x) dlnorm(x[1], log(x[2]), log(2.38), log = T))
      prop_star_to_current <- apply(cbind(theta_star, theta_init), 1, function(x) dlnorm(x[2], log(x[1]), log(2.38), log = T))
      
      ratio <- post_star_log - post_log + prop_star_to_current - prop_current_to_star
      alpha <- sapply(ratio, function(x) min(0, x))
      index <- which(log(N, 0, 1) < alpha)
      theta[index] <- theta_star[index]
      theta[-index] <- theta_init[-index]
    }
    
    # Update X
    x_star <- sapply(x, function(x) rnorm(1, x, 2.38))
    
    prior_log <- apply(cbind(x, theta), 1, function(x) dnorm(x[1], 0, x[2], log = T))
    prior_star_log <- apply(cbind(x_star, theta), 1, function(x) dnorm(x[1], 0, x[2], log = T))
    
    like_log <- sapply(x, function(x) dnorm(data, x, 1, log = T))
    like_star_log <- sapply(x_star, function(x) dnorm(data, x, 1, log = T))
    
    post_log <- prior_log + like_log
    post_star_log <- prior_star_log + like_star_log
    
    ratio <- post_star_log - post_log
    alpha <- sapply(ratio, function(x) min(0, x))
    index <- which(log(runif(N, 0, 1)) < alpha)
    x[index] <- x_star[index]
    
    # Update theta
    theta_star <- sapply(theta, function(x) rlnorm(1, log(x), log(2.38)))
    
    prior_log <- sapply(theta, function(x) dexp(x, prior_lambda, log = T))
    prior_star_log <- sapply(theta_star, function(x) dexp(x, prior_lambda, log = T))
    
    like_log <- apply(cbind(x, theta), 1, function(x) dnorm(x[1], 0, x[2], log = T))
    like_star_log <- apply(cbind(x, theta_star), 1, function(x) dnorm(x[1], 0, x[2], log = T))
    
    post_log <- prior_log + like_log
    post_star_log <- prior_star_log + like_star_log
    
    prop_current_to_star <- apply(cbind(theta_star, theta), 1, function(x) dlnorm(x[1], log(x[2]), log(2.38), log = T))
    prop_star_to_current <- apply(cbind(theta_star, theta), 1, function(x) dlnorm(x[2], log(x[1]), log(2.38), log = T))
    
    ratio <- post_star_log - post_log + prop_star_to_current - prop_current_to_star
    alpha <- sapply(ratio, function(x) min(0, x))
    index <- which(log(N, 0, 1) < alpha)
    theta[index] <- theta_star[index]
  }
  
  
  
  
  return(list(x=x, theta=theta))
}

mcmc_kernel <- function(data, node, theta_init, prior_lambda, N, Ntotal) {
  
  # Storage
  theta <- rep(NA, N)
  
  # Initialization
  if (node == 1) {
    theta_init <- sapply(prior_lambda, function(x) rexp(1, x))
  } else {
    theta_init <- theta_init  
  }
  
  
  for (i in 1:(Ntotal-1)) {
    if (i == 1) {
      theta_star <- sapply(theta_init, function(x) rlnorm(1, log(x), log(2.38)))
      
      prior_log <- apply(cbind(theta_init, prior_lambda), 1, function(x) dexp(x[1], 1/x[2], log = T))
      prior_star_log <- apply(cbind(theta_star, prior_lambda), 1, function(x) dexp(x[1], 1/x[2], log = T))
      
      like_log <- apply(cbind(data, theta_init), 1, function(x) dnorm(x[1], 0, x[2], log = T))
      like_star_log <- apply(cbind(data, theta_star), 1, function(x) dnorm(x[1], 0, x[2], log = T))
      
      post_log <- prior_log + like_log
      post_star_log <- prior_star_log + like_star_log
      
      prop_current_to_star <- apply(cbind(theta_init, theta_star), 1, function(x) dlnorm(x[2], log(x[1]), log(2.38), log = T))
      prop_star_to_current <- apply(cbind(theta_init, theta_star), 1, function(x) dlnorm(x[1], log(x[2]), log(2.38), log = T))
      
      ratio <- post_star_log - post_log + prop_star_to_current - prop_current_to_star
      alpha <- sapply(ratio, function(x) min(0, x))
      index <- which(log(runif(N, 0, 1)) < alpha)
      theta[index] <- theta_star[index]
      theta[-index] <- theta_init[-index]
    }
    
    theta_star <- sapply(theta, function(x) rlnorm(1, log(x), log(2.38)))
    
    prior_log <- apply(cbind(theta, prior_lambda), 1, function(x) dexp(x[1], 1/x[2], log = T))
    prior_star_log <- apply(cbind(theta_star, prior_lambda), 1, function(x) dexp(x[1], 1/x[2], log = T))
    
    like_log <- apply(cbind(data, theta), 1, function(x) dnorm(x[1], 0, x[2], log = T))
    like_star_log <- apply(cbind(data, theta_star), 1, function(x) dnorm(x[1], 0, x[2], log = T))
    
    post_log <- prior_log + like_log
    post_star_log <- prior_star_log + like_star_log
    
    prop_current_to_star <- apply(cbind(theta, theta_star), 1, function(x) dlnorm(x[2], log(x[1]), log(2.38), log = T))
    prop_star_to_current <- apply(cbind(theta, theta_star), 1, function(x) dlnorm(x[1], log(x[2]), log(2.38), log = T))
    
    ratio <- post_star_log - post_log + prop_star_to_current - prop_current_to_star
    alpha <- sapply(ratio, function(x) min(0, x))
    index <- which(log(runif(N, 0, 1)) < alpha)
    theta[index] <- theta_star[index]
  }
  
  return(theta)
}

# MCMC/HM kernel help function on the subroot node
sub_mcmc <- function(data, x_0, prior_alpha, prior_beta, like_beta, Ntotal) {
  N <- length(x_0)
  
  # Storage
  x <- rep(NA, N) # No need to store the samples for all iterations!
  
  for (i in 1:(Ntotal-1)) {
    if (i == 1) {
      # Proposal - log-normal random walk
      x_star <- exp(mvrnorm(1, log(x_0), 100*diag(N))) 
      
      # Prior
      prior_log <- sapply(x_0, function(x) dgamma(x, prior_alpha, prior_beta, log = T))
      prior_star_log <- sapply(x_star, function(x) dgamma(x, prior_alpha, prior_beta, log = T))
      
      # Likelihood
      data_and_x <- cbind(data, x_0)
      like_log <- apply(data_and_x, 1, function(x) dbeta(x[1], x[3], like_beta[1], log = T) + dbeta(x[2], x[3], like_beta[2], log = T))
      data_and_x_star <- cbind(data, x_star)
      like_star_log <- apply(data_and_x_star, 1, function(x) dbeta(x[1], x[3], like_beta[1], log = T) + dbeta(x[2], x[3], like_beta[2], log = T))
      
      # Posterior
      post_log <- prior_log + like_log
      post_star_log <- prior_star_log + like_star_log
      
      # Accept/reject - carry out in a matrix format
      ratio <- post_star_log - post_log + log(x_star) - log(x_0) # Note: consider the Jacobian determinant
      alpha <- sapply(ratio, function(x) min(0, x))
      Unif_rand <- runif(N, 0, 1)
      index <- which(log(Unif_rand) < alpha)
      
      x[index] <- x_star[index]
      x[-index] <- x_0[-index]
    }
    
    
    # Proposal - log-normal random walk
    x_star <- exp(mvrnorm(1, log(x), 100*diag(N)))
    
    # Prior
    prior_log <- sapply(x, function(x) dgamma(x, prior_alpha, prior_beta, log = T))
    prior_star_log <- sapply(x_star, function(x) dgamma(x, prior_alpha, prior_beta, log = T))
    
    # Likelihood
    data_and_x <- cbind(data, x)
    like_log <- apply(data_and_x, 1, function(x) dbeta(x[1], x[3], like_beta[1], log = T) + dbeta(x[2], x[3], like_beta[2], log = T))
    data_and_x_star <- cbind(data, x_star)
    like_log <- apply(data_and_x_star, 1, function(x) dbeta(x[1], x[3], like_beta[1], log = T) + dbeta(x[2], x[3], like_beta[2], log = T))
    
    # Accept/reject - carry out in a matrix format
    ratio <- post_star_log - post_log + log(x_star) - log(x) # Note: consider the Jacobian determinant
    alpha <- sapply(ratio, function(x) min(0, x))
    Unif_rand <- runif(N, 0, 1)
    index <- which(log(Unif_rand) < alpha)
    
    x[index] <- x_star[index]
    # No need to do so: x[-index] <- x[-index] !
  }
  
  return(x_post=x)
}



# MCMC/HM kernel help function on the root node
root_mcmc <- function(data, x_0, prior_mu, prior_sigma, like_beta, Ntotal) {
  N <- length(x_0)
  
  # Storage
  x <- rep(NA, N)
  
  for (i in 1:(Ntotal-1)) {
    if (i == 1) {
      # Proposal - log-normal random walk
      x_star <- exp(mvrnorm(1, log(x_0), 100*diag(N)))
      
      
      # Prior
      prior_log <- sapply(x_0, function(x) dlnorm(x, prior_mu, prior_sigma, log = T))
      prior_star_log <- sapply(x_0, function(x) dlnorm(x, prior_mu, prior_sigma, log = T))
      
      # Likelihood
      data_and_x <- cbind(data, x_0)
      like_log <- apply(data_and_x, 1, function(x) dgamma(x[1], x[3], like_beta[1], log = T) + dgamma(x[2], x[3], like_beta[2], log = T))
      data_and_x_star <- cbind(data, x_star)
      like_star_log <- apply(data_and_x_star, 1, function(x) dgamma(x[1], x[3], like_beta[1], log = T) + dgamma(x[2], x[3], like_beta[2], log = T))
      
      # Posterior
      post_log <- prior_log + like_log
      post_star_log <- prior_star_log + like_star_log
      
      # Accept/reject - carry out in a matrix form
      ratio <- post_star_log - post_log + log(x_star) - log(x_0) # Note: consider the Jacobian determinant
      alpha <- sapply(ratio, function(x) min(0, x))
      Unif_rand <- runif(N, 0, 1)
      index <- which(log(Unif_rand) < alpha)
      
      x[index] <- x_star[index]
      x[-index] <- x_0[-index]
    }
    
    
    # Proposal - log-normal random walk
    x_star <- exp(mvrnorm(1, log(x), 100*diag(N)))
    
    # Prior
    prior_log <- sapply(x, function(x) dlnorm(x, prior_mu, prior_sigma, log = T))
    prior_star_log <- sapply(x, function(x) dlnorm(x_star, prior_mu, prior_sigma, log = T))
    
    # Likelihood
    data_and_x <- cbind(data, x)
    like_log <- apply(data_and_x, 1, function(x) dgamma(x[1], x[3], like_beta[1], log = T) + dgamma(x[2], x[3], like_beta[2], log = T))
    data_and_x_star <- cbind(data, x_star)
    like_star_log <- apply(data_and_x_star, 1, function(x) dgamma(x[1], x[3], like_beta[1], log = T) + dgamma(x[2], x[3], like_beta[2], log = T))
    
    # Accept/reject - carry out in a matrix form
    ratio <- post_star_log - post_log + log(x_star) - log(x) # Note: consider the Jacobian determinant
    alpha <- sapply(ratio, function(x) min(0, x))
    Unif_rand <- runif(N, 0, 1)
    index <- which(log(Unif_rand) < alpha)
    
    x[index] <- x_star[index]
    # No need to do so: x[-index] <- x[-index] !
  }
  
  return(x_post=x)
}

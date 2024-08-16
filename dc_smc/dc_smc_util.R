# MCMC/MH kernel help function on the subroot node
sub_mcmc <- function(data, x_0, prior_alpha, prior_beta, like_beta, var_sub, Ntotal) {
  N <- length(x_0)
  d <- ncol(data)
  
  # Storage
  x <- rep(NA, N) # No need to store the samples for all iterations!
  
  for (i in 1:(Ntotal-1)) {
    if (i == 1) {
      # Proposal - log-normal random walk
      x_star <- sapply(x_0, function(x) rlnorm(1, log(x), sdlog = 10))
      #x_star <- exp(mvrnorm(1, log(x_0), 2.38*diag(var_sub, N))) 
      
      # Prior
      #prior_log <- sapply(x_0, function(x) dgamma(x, prior_alpha, prior_beta, log = T))
      #prior_star_log <- sapply(x_star, function(x) dgamma(x, prior_alpha, prior_beta, log = T))
      prior_log <- sapply(x_0, function(x) dunif(x, 0.001, 10000, log = T))
      prior_star_log <- sapply(x_star, function(x) dunif(x, .001, 10000, log = T))
      
      # Likelihood
      like_log <- likelihood_particle_beta(data, matrix(rep(x_0, d), ncol = d), matrix(rep(like_beta, N), ncol = d, byrow = T))
      like_star_log <- likelihood_particle_beta(data, matrix(rep(x_star, d), ncol = d), matrix(rep(like_beta, N), ncol = d, byrow = T))
      #data_and_x <- cbind(data, x_0)
      #like_log <- apply(data_and_x, 1, function(x) dbeta(x[1], x[3], like_beta[1], log = T) + dbeta(x[2], x[3], like_beta[2], log = T))
      #data_and_x_star <- cbind(data, x_star)
      #like_star_log <- apply(data_and_x_star, 1, function(x) dbeta(x[1], x[3], like_beta[1], log = T) + dbeta(x[2], x[3], like_beta[2], log = T))
      
      # Posterior
      post_log <- prior_log + like_log
      post_star_log <- prior_star_log + like_star_log
      
      # Transition/proposal density
      x_and_x_star <- cbind(x_0, x_star)
      prop_x_to_x_star <- apply(x_and_x_star, 1, function(x) dlnorm(x[2], log(x[1]), sdlog = 10, log = T))
      prop_x_star_to_x <- apply(x_and_x_star, 1, function(x) dlnorm(x[1], log(x[2]), sdlog = 10, log = T))
      
      # Accept/reject - carry out in a matrix format
      ratio <- post_star_log - post_log + prop_x_star_to_x - prop_x_to_x_star
      alpha <- sapply(ratio, function(x) min(0, x))
      Unif_rand <- runif(N, 0, 1)
      index <- which(log(Unif_rand) < alpha)
      
      x[index] <- x_star[index]
      x[-index] <- x_0[-index]
    }
    
    
    # Proposal - log-normal random walk
    x_star <- sapply(x, function(x) rlnorm(1, log(x), sdlog = sqrt(2.38*var_sub)))
    #x_star <- exp(mvrnorm(1, log(x), 2.38*diag(var_sub, N)))
    
    # Prior
    #prior_log <- sapply(x, function(x) dgamma(x, prior_alpha, prior_beta, log = T))
    #prior_star_log <- sapply(x_star, function(x) dgamma(x, prior_alpha, prior_beta, log = T))
    prior_log <- sapply(x, function(x) dunif(x, .001, 10000, log = T))
    prior_star_log <- sapply(x_star, function(x) dunif(x, .001, 10000, log = T))
    
    # Likelihood
    like_log <- likelihood_particle_beta(data, matrix(rep(x, d), ncol = d), matrix(rep(like_beta, N), ncol = d, byrow = T))
    like_star_log <- likelihood_particle_beta(data, matrix(rep(x_star, d), ncol = d), matrix(rep(like_beta, N), ncol = d, byrow = T))
    #data_and_x <- cbind(data, x)
    #like_log <- apply(data_and_x, 1, function(x) dbeta(x[1], x[3], like_beta[1], log = T) + dbeta(x[2], x[3], like_beta[2], log = T))
    #data_and_x_star <- cbind(data, x_star)
    #like_star_log <- apply(data_and_x_star, 1, function(x) dbeta(x[1], x[3], like_beta[1], log = T) + dbeta(x[2], x[3], like_beta[2], log = T))
    
    # Posterior
    post_log <- prior_log + like_log
    post_star_log <- prior_star_log + like_star_log
    
    # Transition/proposal density
    x_and_x_star <- cbind(x, x_star)
    prop_x_to_x_star <- apply(x_and_x_star, 1, function(x) dlnorm(x[2], log(x[1]), sdlog = sqrt(2.38*var_sub), log = T))
    prop_x_star_to_x <- apply(x_and_x_star, 1, function(x) dlnorm(x[1], log(x[2]), sdlog = sqrt(2.38*var_sub), log = T))
    
    # Accept/reject - carry out in a matrix format
    ratio <- post_star_log - post_log + prop_x_star_to_x - prop_x_to_x_star
    alpha <- sapply(ratio, function(x) min(0, x))
    Unif_rand <- runif(N, 0, 1)
    index <- which(log(Unif_rand) < alpha)
    
    x[index] <- x_star[index]
    # No need to do so: x[-index] <- x[-index] !
  }
  
  return(x_post=x)
}



# MCMC/MH kernel help function on the root node
root_mcmc <- function(data, x_0, prior_mu, prior_sigma, like_beta, var_root, Ntotal) {
  N <- length(x_0)
  
  # Storage
  x <- rep(NA, N)
  
  for (i in 1:(Ntotal-1)) {
    if (i == 1) {
      # Proposal - log-normal random walk
      x_star <- sapply(x_0, function(x) rlnorm(1, log(x), 2.38*sqrt(var_root)))
      #x_star <- exp(mvrnorm(1, log(x_0), 10*diag(N)))
      
      
      # Prior
      #prior_log <- sapply(x_0, function(x) dlnorm(x, prior_mu, prior_sigma, log = T))
      #prior_star_log <- sapply(x_star, function(x) dlnorm(x, prior_mu, prior_sigma, log = T))
      prior_log <- sapply(x_0, function(x) dunif(x, .001, 10000, log = T))
      prior_star_log <- sapply(x_star, function(x) dunif(x, .001, 10000, log = T))
      
      # Likelihood
      data_and_x <- cbind(data, x_0)
      like_log <- apply(data_and_x, 1, function(x) dgamma(x[1], x[3], like_beta[1], log = T) + dgamma(x[2], x[3], like_beta[2], log = T))
      data_and_x_star <- cbind(data, x_star)
      like_star_log <- apply(data_and_x_star, 1, function(x) dgamma(x[1], x[3], like_beta[1], log = T) + dgamma(x[2], x[3], like_beta[2], log = T))
      
      # Posterior
      post_log <- prior_log + like_log
      post_star_log <- prior_star_log + like_star_log
      
      # Transition/proposal density
      x_and_x_star <- cbind(x_0, x_star)
      prop_x_to_x_star <- apply(x_and_x_star, 1, function(x) dlnorm(x[2], log(x[1]), 2.38*sqrt(var_root), log = T))
      prop_x_star_to_x <- apply(x_and_x_star, 1, function(x) dlnorm(x[1], log(x[2]), 2.38*sqrt(var_root), log = T))
      
      # Accept/reject - carry out in a matrix form
      ratio <- post_star_log - post_log + prop_x_star_to_x - prop_x_to_x_star
      alpha <- sapply(ratio, function(x) min(0, x))
      Unif_rand <- runif(N, 0, 1)
      index <- which(log(Unif_rand) < alpha)
      
      x[index] <- x_star[index]
      x[-index] <- x_0[-index]
    }
    
    
    # Proposal - log-normal random walk
    x_star <- sapply(x, function(x) rlnorm(1, log(x), 2.38*sqrt(var_root)))
    #x_star <- exp(mvrnorm(1, log(x), 10*diag(N)))
    
    # Prior
    #prior_log <- sapply(x, function(x) dlnorm(x, prior_mu, prior_sigma, log = T))
    #prior_star_log <- sapply(x_star, function(x) dlnorm(x, prior_mu, prior_sigma, log = T))
    prior_log <- sapply(x, function(x) dunif(x, .001, 10000, log = T))
    prior_star_log <- sapply(x_star, function(x) dunif(x, .001, 10000, log = T))
    
    # Likelihood
    data_and_x <- cbind(data, x)
    like_log <- apply(data_and_x, 1, function(x) dgamma(x[1], x[3], like_beta[1], log = T) + dgamma(x[2], x[3], like_beta[2], log = T))
    data_and_x_star <- cbind(data, x_star)
    like_star_log <- apply(data_and_x_star, 1, function(x) dgamma(x[1], x[3], like_beta[1], log = T) + dgamma(x[2], x[3], like_beta[2], log = T))
    
    # Posterior
    post_log <- prior_log + like_log
    post_star_log <- prior_star_log + like_star_log
    
    # Transition/proposal density
    x_and_x_star <- cbind(x, x_star)
    prop_x_to_x_star <- apply(x_and_x_star, 1, function(x) dlnorm(x[2], log(x[1]), 2.38*sqrt(var_root), log = T))
    prop_x_star_to_x <- apply(x_and_x_star, 1, function(x) dlnorm(x[1], log(x[2]), 2.38*sqrt(var_root), log = T))
    
    # Accept/reject - carry out in a matrix form
    ratio <- post_star_log - post_log + prop_x_star_to_x - prop_x_to_x_star
    alpha <- sapply(ratio, function(x) min(0, x))
    Unif_rand <- runif(N, 0, 1)
    index <- which(log(Unif_rand) < alpha)
    
    x[index] <- x_star[index]
    # No need to do so: x[-index] <- x[-index] !
  }
  
  return(x_post=x)
}

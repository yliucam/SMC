# MCMC on leaf nodes
dc_melding_leaf_mcmc <- function(data, N, Ntotal) {
  #data_rep <- matrix(rep(data, N), nrow = N, byrow = T)
  
  mu <- array(rep(NA, N*Ntotal), dim = c(N, Ntotal))
  sigma <- rep(NA, N)
  
  mu_0 <- rep(0, N)
  sigma_0 <- rep(1, N)
  
  for (i in 1:(Ntotal-1)) {
    if (i == 1) {
      sigma_star <- rlnorm(N, log(sigma_0[1]), sdlog = sqrt(2.38))
      
      prior_sigma_log <- sapply(sigma_0, function(x) dgamma(x, 1, 1, log = T))
      prior_sigma_star_log <- sapply(sigma_star, function(x) dgamma(x, 1, 1, log = T))
      
      likelihood_sigma_log <- rep(NA, N)
      likelihood_sigma_star_log <- rep(NA, N)
      for (j in 1:N) {
        likelihood_sigma_log[j] <- sum(sapply(data, function(x) dnorm(x, mu_0[j], sigma_0[j], log = T)))
        likelihood_sigma_star_log[j] <- sum(sapply(data, function(x) dnorm(x, mu_0[j], sigma_star[j], log = T)))
      }
      
      post_sigma_log <- prior_sigma_log + likelihood_sigma_log
      post_sigma_star_log <- prior_sigma_star_log + likelihood_sigma_star_log
      
      prop_sigma_to_sigma_star <- apply(cbind(sigma_0, sigma_star), 1, function(x) dlnorm(x[2], log(x[1]), sdlog = sqrt(2.38), log = T))
      prop_sigma_star_to_sigma <- apply(cbind(sigma_0, sigma_star), 1, function(x) dlnorm(x[1], log(x[2]), sdlog = sqrt(2.38), log = T))
      
      ratio <- post_sigma_star_log - post_sigma_log + prop_sigma_star_to_sigma - prop_sigma_to_sigma_star
      alpha <- sapply(ratio, function(x) min(0, x))
      Unif_rand <- runif(N, 0, 1)
      index <- which(log(Unif_rand) < alpha)
      
      sigma[index] <- sigma_star[index]
      sigma[-index] <- sigma_0[-index]
      
      
      mu_star <- rnorm(N, mu_0[1], 10)
      
      prior_mu_log <- sapply(mu_0, function(x) dnorm(x, 0, 1, log = T))
      prior_mu_star_log <- sapply(mu_star, function(x) dnorm(x, 0, 1, log = T))
      
      likelihood_mu_log <- rep(NA, N)
      likelihood_mu_star_log <- rep(NA, N)
      for (j in 1:N) {
        likelihood_mu_log[j] <- sum(sapply(data, function(x) dnorm(x, mu_0[j], sigma[j], log = T)))
        likelihood_mu_star_log[j] <- sum(sapply(data, function(x) dnorm(x, mu_star[j], sigma[j], log = T)))
      }
      
      post_mu_log <- prior_mu_log + likelihood_mu_log
      post_mu_star_log <- prior_mu_star_log + likelihood_mu_star_log
      
      #ratio <- post_mu_star_log - post_mu_log
      ratio <- likelihood_mu_star_log - likelihood_mu_log
      alpha <- sapply(ratio, function(x) min(0, x))
      Unif_rand <- runif(N, 0, 1)
      index <- which(log(Unif_rand) < alpha)
      
      if (length(index) == 0) {
        mu[,i] <- mu_0
      } else {
        mu[index, i] <- mu_star[index]
        mu[-index, i] <- mu_0[-index]  
      }
    }
    
    
    sigma_star <- sapply(sigma, function(x) rlnorm(1, log(x), sdlog = sqrt(2.38)))
    
    prior_sigma_log <- sapply(sigma, function(x) dgamma(x, 1, 1, log = T))
    prior_sigma_star_log <- sapply(sigma_star, function(x) dgamma(x, 1, 1, log = T))
    
    likelihood_sigma_log <- rep(NA, N)
    likelihood_sigma_star_log <- rep(NA, N)
    for (j in 1:N) {
      likelihood_sigma_log[j] <- sum(sapply(data, function(x) dnorm(x, mu[j,i], sigma[j], log = T)))
      likelihood_sigma_star_log[j] <- sum(sapply(data, function(x) dnorm(x, mu[j,i], sigma_star[j], log = T)))
    }
    
    post_sigma_log <- prior_sigma_log + likelihood_sigma_log
    post_sigma_star_log <- prior_sigma_star_log + likelihood_sigma_star_log
    
    prop_sigma_to_sigma_star <- apply(cbind(sigma, sigma_star), 1, function(x) dlnorm(x[2], log(x[1]), sdlog = sqrt(2.38), log = T))
    prop_sigma_star_to_sigma <- apply(cbind(sigma, sigma_star), 1, function(x) dlnorm(x[1], log(x[2]), sdlog = sqrt(2.38), log = T))
    
    ratio <- post_sigma_star_log - post_sigma_log + prop_sigma_star_to_sigma - prop_sigma_to_sigma_star
    alpha <- sapply(ratio, function(x) min(0, x))
    Unif_rand <- runif(N, 0, 1)
    index <- which(log(Unif_rand) < alpha)
    
    sigma[index] <- sigma_star[index]
    
    
    mu_star <- sapply(mu[,i], function(x) rnorm(1, x, 10))
    
    prior_mu_log <- sapply(mu[,i], function(x) dnorm(x, 0, 1, log = T))
    prior_mu_star_log <- sapply(mu_star, function(x) dnorm(x, 0, 1, log = T))
    
    likelihood_mu_log <- rep(NA, N)
    likelihood_mu_star_log <- rep(NA, N)
    for (j in 1:N) {
      likelihood_mu_log[j] <- sum(sapply(data, function(x) dnorm(x, mu[j,i], sigma[j], log = T)))
      likelihood_mu_star_log[j] <- sum(sapply(data, function(x) dnorm(x, mu_star[j], sigma[j], log = T)))
    }
    
    post_mu_log <- prior_mu_log + likelihood_mu_log
    post_mu_star_log <- prior_mu_star_log + likelihood_mu_star_log
    
    #ratio <- post_mu_star_log - post_mu_log
    ratio <- likelihood_mu_star_log - likelihood_mu_log
    alpha <- sapply(ratio, function(x) min(0, x))
    Unif_rand <- runif(N, 0, 1)
    index <- which(log(Unif_rand) < alpha)
    
    if (length(index) == 0) {
      mu[,i+1] <- mu[,i]
    } else {
      mu[index, i+1] <- mu_star[index]
      mu[-index, i+1] <- mu[-index, i]  
    }
  }
  
  w_log <- apply(cbind(data, mu[,Ntotal], sigma), 1, function(x) dnorm(x[1], x[2], x[3], log = T))
  min_lim <- log(.Machine$double.xmin)
  w_log[which(w_log < min_lim)] <- min_lim
  W <- exp(w_log - matrixStats::logSumExp(w_log))
  
  return(list(mu=mu[,Ntotal], sigma=sigma, w_log=w_log, W=W))
}






# MCMC on root
dc_melding_root_mcmc <- function(data, data_phi, Ntotal, mu_init, phi_sigma, alpha_j) {
  N <- dim(data_phi)[1]
  d <- dim(data_phi)[2]
  mu <- array(rep(NA, N*Ntotal), dim = c(N, Ntotal))
  for (i in 1:(Ntotal-1)) {
    if (i == 1) {
      mu_star <- sapply(mu_init, function(x) rlnorm(1, x, sdlog = sqrt(2.38)))
      
      prior_log <- sapply(mu_init, function(x) dgamma(x, 1, 1, log = T))
      prior_star_log <- sapply(mu_star, function(x) dgamma(x, 1, 1, log = T))
      
      likelihood_phi_log <- normal_nodes_log_sum_new(data_phi, mu = matrix(rep(mu_init, d), ncol=d), sigma = matrix(rep(phi_sigma, N), ncol=d, byrow=T))
      likelihood_phi_star_log <- normal_nodes_log_sum_new(data_phi, mu = matrix(rep(mu_star, d), ncol=d), sigma = matrix(rep(phi_sigma, N), ncol=d, byrow=T))
      likelihood_y_log <- exponential_psd_log_sum(data, lambda = mu_init)
      likelihood_y_star_log <- exponential_psd_log_sum(data, lambda = mu_star)
      
      post_log <- prior_log + alpha_j * (likelihood_phi_log + likelihood_y_log)
      post_star_log <- prior_star_log + alpha_j * (likelihood_phi_star_log + likelihood_y_star_log)
      
      prop_mu_to_mu_star <- apply(cbind(mu_init, mu_star), 1, function(x) dlnorm(x[2], log(x[1]), sdlog = sqrt(2.38), log = T))
      prop_mu_star_to_mu <- apply(cbind(mu_init, mu_star), 1, function(x) dlnorm(x[1], log(x[2]), sdlog = sqrt(2.38), log = T))
      
      ratio <- post_star_log - post_log + prop_mu_star_to_mu - prop_mu_to_mu_star
      alpha <- sapply(ratio, function(x) min(0, x))
      Unif_rand <- runif(N, 0, 1)
      index <- which(log(Unif_rand) < alpha)
      
      if (length(index) == 0) {
        mu[,i] <- mu_init
      } else {
        mu[index, i] <- mu_star[index]
        mu[-index, i] <- mu_init[-index]
      }
    }
    
    mu_star <- sapply(mu[,i], function(x) rlnorm(1, x, sdlog = sqrt(2.38)))
    
    prior_log <- sapply(mu[,i], function(x) dgamma(x, 1, 1, log = T))
    prior_star_log <- sapply(mu_star, function(x) dgamma(x, 1, 1, log = T))
    
    likelihood_phi_log <- normal_nodes_log_sum_new(data_phi, mu = matrix(rep(mu[,i], d), ncol=d), sigma = matrix(rep(phi_sigma, N), ncol=d, byrow=T))
    likelihood_phi_star_log <- normal_nodes_log_sum_new(data_phi, mu = matrix(rep(mu_star, d), ncol=d), sigma = matrix(rep(phi_sigma, N), ncol=d, byrow=T))
    likelihood_y_log <- exponential_psd_log_sum(data, lambda = mu[,i])
    likelihood_y_star_log <- exponential_psd_log_sum(data, lambda = mu_star)
    
    post_log <- prior_log + alpha_j * (likelihood_phi_log + likelihood_y_log)
    post_star_log <- prior_star_log + alpha_j * (likelihood_phi_star_log + likelihood_y_star_log)
    
    prop_mu_to_mu_star <- apply(cbind(mu, mu_star), 1, function(x) dlnorm(x[2], log(x[1]), sdlog = sqrt(2.38), log = T))
    prop_mu_star_to_mu <- apply(cbind(mu, mu_star), 1, function(x) dlnorm(x[1], log(x[2]), sdlog = sqrt(2.38), log = T))
    
    ratio <- post_star_log - post_log + prop_mu_star_to_mu - prop_mu_to_mu_star
    alpha <- sapply(ratio, function(x) min(0, x))
    Unif_rand <- runif(N, 0, 1)
    index <- which(log(Unif_rand) < alpha)
    
    if (length(index) == 0) {
      mu[,i+1] <- mu[,i]
    } else {
      mu[index, i+1] <- mu_star[index]
      mu[-index, i+1] <- mu[-index, i]
    }
  }
  
  w_log <- apply(cbind(data_phi, mu[,Ntotal]), 1, function(x) dnorm(x[1], x[3], 1, log = T) + dnorm(x[2], x[3], 2, log = T))
  W <- exp(w_log - matrixStats::logSumExp(w_log))
  
  return(list(mu=mu, w_log=w_log, W=W))
}

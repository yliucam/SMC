dc_melding_smc_mcmc <- function(data,
                                phi,
                                psi_alpha,
                                psi_beta,
                                sigma_init,
                                alpha_j,
                                N,
                                Ntotal) {
  
  sigma <- array(rep(NA, N*Ntotal), dim = c(N, Ntotal))
  
  for (i in 1:(Ntotal-1)) {
    if (i == 1) {
      sigma_star <- sapply(sigma_init, function(x) rlnorm(1, x, sdlog = sqrt(2.38)))
      
      prior_log <- sapply(sigma_init, function(x) dgamma(x, psi_alpha, psi_beta, log = T))
      prior_star_log <- sapply(sigma_star, function(x) dgamma(x, psi_alpha, psi_beta, log = T))
      
      likelihood_log <- normal_uni_log_sum(data = data,
                                           mu = phi,
                                           sigma = sigma_init)
      likelihood_star_log <- normal_uni_log_sum(data = data,
                                                mu = phi,
                                                sigma = sigma_star)
      post_log <- prior_log + alpha_j * likelihood_log
      post_star_log <- prior_star_log + alpha_j * likelihood_star_log
      
      prop_sigma_to_sigma_star <- apply(cbind(sigma_init, sigma_star), 1, function(x) dlnorm(x[2], log(x[1]), sdlog = sqrt(2.38), log = T))
      prop_sigma_star_to_sigma <- apply(cbind(sigma_init, sigma_star), 1, function(x) dlnorm(x[1], log(x[2]), sdlog = sqrt(2.38), log = T))
      
      ratio <- post_star_log - post_log + prop_sigma_star_to_sigma - prop_sigma_to_sigma_star
      alpha <- sapply(ratio, function(x) min(0, x))
      Unif_rand <- runif(N, 0, 1)
      index <- which(log(Unif_rand) < alpha)
      
      if (length(index) == 0) {
        sigma[,i] <- sigma_init
      } else {
        sigma[index, i] <- sigma_star[index]
        sigma[-index, i] <- sigma_init[-index]
      }
    }
    
    sigma_star <- sapply(sigma[,i], function(x) rlnorm(1, x, sdlog = sqrt(2.38)))
    
    prior_log <- sapply(sigma[,i], function(x) dgamma(x, psi_alpha, psi_beta, log = T))
    prior_star_log <- sapply(sigma_star, function(x) dgamma(x, psi_alpha, psi_beta, log = T))
    
    likelihood_log <- normal_uni_log_sum(data = data,
                                         mu = phi,
                                         sigma = sigma[,i])
    likelihood_star_log <- normal_uni_log_sum(data = data,
                                              mu = phi,
                                              sigma = sigma_star)
    
    post_log <- prior_log + alpha_j * likelihood_log
    post_star_log <- prior_star_log + alpha_j * likelihood_star_log
    
    prop_sigma_to_sigma_star <- apply(cbind(sigma[,i], sigma_star), 1, function(x) dlnorm(x[2], log(x[1]), sdlog = sqrt(2.38), log = T))
    prop_sigma_star_to_sigma <- apply(cbind(sigma[,i], sigma_star), 1, function(x) dlnorm(x[1], log(x[2]), sdlog = sqrt(2.38), log = T))
    
    ratio <- post_star_log - post_log + prop_sigma_star_to_sigma - prop_sigma_to_sigma_star
    alpha <- sapply(ratio, function(x) min(0, x))
    Unif_rand <- runif(N, 0, 1)
    index <- which(log(Unif_rand) < alpha)
    
    if (length(index) == 0) {
      sigma[,i+1] <- sigma[,i]
    } else {
      sigma[index, i+1] <- sigma_star[index]
      sigma[-index, i+1] <- sigma[-index, i]
    }
  }
  
  return(list(sigma=sigma))
}







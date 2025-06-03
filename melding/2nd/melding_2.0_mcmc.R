melding_2.0_mcmc <- function(data,
                             Ntotal,
                             burnin,
                             sub,
                             lambda,
                             N) {
  phi12 <- rep(NA, Ntotal)
  phi23 <- rep(NA, Ntotal)
  
  psi2 <- rep(NA, Ntotal)
  
  phi12_prop <- sub[,1]
  phi23_prop <- sub[,2]
  
  phi12_init <- phi12_prop[N]
  phi23_init <- phi23_prop[N]
  
  psi2_init <- 1
  
  for (i in 1:(Ntotal-1)) {
    if (i == 1) {
      # Update phi12
      index1 <- sample(N, 1)
      phi12_star <- phi12_prop[index1]
      
      prior_pool_2_log <- pooling_prior_log_log(lambda = lambda,
                                                phi12 = phi12_init,
                                                phi23 = phi23_init)
      prior_pool_2_star_log <- pooling_prior_log_log(lambda = lambda,
                                                     phi12 = phi12_star,
                                                     phi23 = phi23_init)
      
      likelihood_log <- normal_log_sum(data = data,
                                       mu = phi12_init + phi23_init,
                                       sigma = psi2_init)
      likelihood_star_log <- normal_log_sum(data = data,
                                            mu = phi12_star + phi23_init,
                                            sigma = psi2_init)
      
      post_log <- prior_pool_2_log + likelihood_log
      post_star_log <- prior_pool_2_star_log + likelihood_star_log
      
      ratio_12 <- post_star_log - post_log
      alpha_12 <- min(0, ratio_12)
      unif_rand <- runif(1, 0, 1)
      
      if (log(unif_rand) < alpha_12) {
        phi12[i] <- phi12_star
      } else {
        phi12[i] <- phi12_init
      }
      
      # Update phi23
      index2 <- sample(N, 1)
      phi23_star <- phi23_prop[index2]
      
      prior_pool_2_log <- pooling_prior_log_log(lambda = lambda,
                                                phi12 = phi12[i],
                                                phi23 = phi23_init)
      prior_pool_2_star_log <- pooling_prior_log_log(lambda = lambda,
                                                     phi12 = phi12[i],
                                                     phi23 = phi23_star)
      
      likelihood_log <- normal_log_sum(data = data,
                                       mu = phi12[i] + phi23_init,
                                       sigma = psi2_init)
      likelihood_star_log <- normal_log_sum(data = data,
                                            mu = phi12[i] + phi23_star,
                                            sigma = psi2_init)
      
      post_log <- prior_pool_2_log + likelihood_log
      post_star_log <- prior_pool_2_star_log + likelihood_star_log
      
      ratio_23 <- post_star_log - post_log
      alpha_23 <- min(0, ratio_23)
      unif_rand <- runif(1, 0, 1)
      if (log(unif_rand) < alpha_23) {
        phi23[i] <- phi23_star
      } else {
        phi23[i] <- phi23_init
      }
      
      # Update psi2
      psi2_star <- rlnorm(1, psi2_init, sdlog = sqrt(2.38*10))
      
      prior_log <- dgamma(psi2_init, 1, 1, log = T)
      prior_star_log <- dgamma(psi2_star, 1, 1, log = T)
      
      likelihood_log <- normal_log_sum(data = data,
                                       mu = phi12[i] + phi23[i],
                                       sigma = psi2_init)
      likelihood_star_log <- normal_log_sum(data = data,
                                            mu = phi12[i] + phi23[i],
                                            sigma = psi2_star)
      
      post_log <- prior_log + likelihood_log
      post_star_log <- prior_star_log + likelihood_star_log
      
      prop_psi2_to_psi2_star <- dlnorm(psi2_star, log(psi2_init), sdlog = sqrt(2.38*10), log = T)
      prop_psi2_star_to_psi2 <- dlnorm(psi2_init, log(psi2_star), sdlog = sqrt(2.38*10), log = T)
      
      ratio <- post_star_log - post_log + prop_psi2_star_to_psi2 - prop_psi2_to_psi2_star
      alpha <- min(0, ratio)
      unif_rand <- runif(1, 0, 1)
      if (log(unif_rand) < alpha) {
        psi2[i] <- psi2_star
      } else {
        psi2[i] <- psi2_init
      }
    }
    
    # Update phi12
    index1 <- sample(N, 1)
    phi12_star <- phi12_prop[index1]
    
    prior_pool_2_log <- pooling_prior_log_log(lambda = lambda,
                                              phi12 = phi12[i],
                                              phi23 = phi23[i])
    prior_pool_2_star_log <- pooling_prior_log_log(lambda = lambda,
                                                   phi12 = phi12_star,
                                                   phi23 = phi23[i])
    
    likelihood_log <- normal_log_sum(data = data,
                                     mu = phi12[i] + phi23[i],
                                     sigma = psi2[i])
    likelihood_star_log <- normal_log_sum(data = data,
                                          mu = phi12_star + phi23[i],
                                          sigma = psi2[i])
    
    post_log <- prior_pool_2_log + likelihood_log
    post_star_log <- prior_pool_2_star_log + likelihood_star_log
    
    ratio_12 <- post_star_log - post_log
    alpha_12 <- min(0, ratio_12)
    unif_rand <- runif(1, 0, 1)
    
    if (log(unif_rand) < alpha_12) {
      phi12[i+1] <- phi12_star
    } else {
      phi12[i+1] <- phi12[i]
    }
    
    # Update phi23
    index2 <- sample(N, 1)
    phi23_star <- phi23_prop[index2]
    
    prior_pool_2_log <- pooling_prior_log_log(lambda = lambda,
                                              phi12 = phi12[i+1],
                                              phi23 = phi23[i])
    prior_pool_2_star_log <- pooling_prior_log_log(lambda = lambda,
                                                   phi12 = phi12[i+1],
                                                   phi23 = phi23_star)
    
    likelihood_log <- normal_log_sum(data = data,
                                     mu = phi12[i+1] + phi23[i],
                                     sigma = psi2[i])
    likelihood_star_log <- normal_log_sum(data = data,
                                          mu = phi12[i+1] + phi23_star,
                                          sigma = psi2[i])
    
    post_log <- prior_pool_2_log + likelihood_log
    post_star_log <- prior_pool_2_star_log + likelihood_star_log
    
    ratio_23 <- post_star_log - post_log
    alpha_23 <- min(0, ratio_23)
    unif_rand <- runif(1, 0, 1)
    if (log(unif_rand) < alpha_23) {
      phi23[i+1] <- phi23_star
    } else {
      phi23[i+1] <- phi23[i]
    }
    
    # Update psi2
    psi2_star <- rlnorm(1, psi2[i], sdlog = sqrt(2.38*10))
    
    prior_log <- dgamma(psi2[i], 1, 1, log = T)
    prior_star_log <- dgamma(psi2_star, 1, 1, log = T)
    
    likelihood_log <- normal_log_sum(data = data,
                                     mu = phi12[i+1] + phi23[i+1],
                                     sigma = psi2[i])
    likelihood_star_log <- normal_log_sum(data = data,
                                          mu = phi12[i+1] + phi23[i+1],
                                          sigma = psi2_star)
    
    post_log <- prior_log + likelihood_log
    post_star_log <- prior_star_log + likelihood_star_log
    
    prop_psi2_to_psi2_star <- dlnorm(psi2_star, log(psi2[i]), sdlog = sqrt(2.38*10), log = T)
    prop_psi2_star_to_psi2 <- dlnorm(psi2[i], log(psi2_star), sdlog = sqrt(2.38*10), log = T)
    
    ratio <- post_star_log - post_log + prop_psi2_star_to_psi2 - prop_psi2_to_psi2_star
    alpha <- min(0, ratio)
    unif_rand <- runif(1, 0, 1)
    if (log(unif_rand) < alpha) {
      psi2[i+1] <- psi2_star
    } else {
      psi2[i+1] <- psi2[i]
    }
  }
  
  keep <- (burnin+1):Ntotal
  
  return(list(phi12=phi12[keep], phi23=phi23[keep], psi2=psi2[keep]))
}


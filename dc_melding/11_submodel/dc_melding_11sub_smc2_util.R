pf_sub6_update <- function(data,
                           theta_curr,
                           w_theta_log_curr,
                           W_theta_curr,
                           mu,
                           sigma,
                           rho,
                           Nparam, 
                           Nstate,
                           t) {
  ESS_min <- Nstate / 2
  min_lim <- log(.Machine$double.xmin)
  
  theta_next <- matrix(rep(NA, Nstate*Nparam), nrow = Nstate)
  w_theta_log <- W_theta <- matrix(rep(NA, Nstate*Nparam), nrow = Nstate)
  llike_sum <- rep(NA, Nparam)
  
  for (i in 1:Nparam) {
    ESS <- 1 / sum(W_theta_curr[,i]^2)
    if (ESS < ESS_min) {
      U <- runif(1, 0, 1)
      A <- Sys_resamp(W = W_theta_curr[,i], P = Nstate, U = U)
      theta_curr[,i] <- theta_curr[A,i]
      w_theta_log_curr[,i] <- rep(0, Nstate)
      W_theta_curr[,i] <- rep(1, Nstate)
    }
    
    theta_mean <- mu[i] + rho[i] * (theta_curr[,i] - mu[i])
    theta_next[,i] <- MASS::mvrnorm(1, mu = theta_mean, Sigma = diag(sigma[i]^2, Nstate))
    sigma_sd <- exp(theta_next[,i]/2)
    llike <- dnorm_log_uni(data = rep(data, Nstate),
                           mu = rep(0, Nstate),
                           sigma = sigma_sd)
    w_theta_log[,i] <- w_theta_log_curr[,i] + llike
    if (max(w_theta_log[,i]) < min_lim) w_theta_log[,i] <- w_theta_log[,i]/1e4
    w_theta_log[which(w_theta_log[,i] < min_lim),i] <- min_lim
    W_theta[,i] <- exp(w_theta_log[,i] - matrixStats::logSumExp(w_theta_log[,i]))
    llike_sum[i] <- sum(llike)
  }
  
  return(list(theta_next = theta_next,
              w_theta_log = w_theta_log,
              W_theta = W_theta,
              llike_sum = llike_sum,
              theta_curr = theta_curr))
}


pooled_prior_log_sub6 <- function(phi56,
                                  phi67,
                                  lambda) {
  lambda1 <- lambda[1]
  lambda2 <- lambda[2]
  lambda3 <- lambda[3]
  
  N <- length(phi56)
  
  res <- (lambda1 + lambda2 - 1) * dgamma_log(phi56, rep(2, N), rep(2, N))
  res <- (lambda2 + lambda3 - 1) * dnorm_log_uni(phi67, rep(0, N), rep(2, N))
  
  return(res)
}


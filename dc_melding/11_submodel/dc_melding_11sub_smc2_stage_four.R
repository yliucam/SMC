smc2_sub6 <- function(data,
                      phi56,
                      phi67,
                      Nstate,
                      Nparam,
                      n_sub,
                      lambda,
                      Ntotal,
                      TT,
                      ncore) {
  ESS_min <- Nparam / 2
  min_lim <- log(.Machine$double.xmin)
  
  mu_array <- sigma_array <- rho_array <- array(rep(NA, Nparam*(TT+1)), dim = c(Nparam, TT+1)) 
  w_log_array <- W_array <- array(rep(NA, (TT+1)*Nparam), dim = c(Nparam, TT+1))
  theta_array <- array(rep(NA, (TT+1)*Nstate*Nparam), dim = c(Nstate, Nparam, TT+1))
  w_theta_log_array <- W_theta_array <- array(rep(NA, (TT+1)*Nstate*Nparam), dim = c(Nstate, Nparam, TT+1))
  A5_array <- A7_array <- array(rep(NA, Nparam*(TT+1)), dim = c(Nparam, TT+1))
  
  # Initialize common parameter particles from submodels' samples
  sigma_0 <- phi56
  mu_0 <- phi67
  
  # Store the original submodel trajectory
  A5_0 <- A7_0 <- 1:Nparam
  
  phi_sub <- cbind(sigma_0, mu_0)
  
  ## Particle trajectory
  phi_sub_index <- cbind(A5_0, A7_0)
  
  # Merging
  sigma_merged <- phi56
  sigma_array[,1] <- sigma_merged
  mu_merged <- phi67
  mu_array[,1] <- mu_merged
  A5_array[,1] <- A5_0
  A7_array[,1] <- A7_0
  
  u_pooling_log <- pooled_prior_log_sub6(phi56 = sigma_merged,
                                         phi67 = mu_merged,
                                         lambda)
  
  # Initialize parameter particles from priors
  rho_0 <- rbeta(Nparam, 9, 1)
  rho_array[,1] <- rho_0
  
  # Initialized latent variables
  theta_0 <- matrix(rep(NA, Nstate*Nparam), nrow = Nstate)
  for (i in 1:Nparam) {
    theta_0[,i] <- rnorm(Nstate, mu_merged[i], sigma_merged[i])
  }
  theta_array[,,1] <- theta_0
  
  
  # Initial weights from PF likelihoods
  w_theta_log_0 <- W_theta_0 <- matrix(rep(0, Nstate*Nparam), nrow = Nstate)
  for (i in 1:Nparam) {
    w_theta_log_0[,i] <- dnorm_log_uni(data = theta_0[,i],
                                       mu = rep(mu_merged[i], Nstate),
                                       sigma = rep(sigma_merged[i], Nstate))
    if (max(w_theta_log_0[,i]) < min_lim) w_theta_log_0[,i] <- w_theta_log_0[,i]/1e4
    w_theta_log_0[which(w_theta_log_0[,i] < min_lim),i] <- min_lim
    w_theta_log_array[,i,1] <- w_theta_log_0[,i]
    W_theta_0[,i] <- exp(w_theta_log_0[,i] - matrixStats::logSumExp(w_theta_log_0[,i]))
    W_theta_array[,i,1] <- W_theta_0[,i]
  }
  
  w_log_0 <- colSums(w_theta_log_0) / Nstate
  if (max(w_log_0) < min_lim) w_log_0 <- w_log_0/1e4
  w_log_0[which(w_log_0 < min_lim)] <- min_lim
  w_log_array[,1] <- w_log_0
  W_0 <- exp(w_log_0 - matrixStats::logSumExp(w_log_0))
  W_array[,1] <- W_0
  
  pb <- txtProgressBar(min = 0, max = TT, style = 3)
  
  for (t in 1:(TT-1)) {
    if (t == 1) {
      ESS <- 1 / sum(W_0^2)
      # rejuvenation
      if (ESS < ESS_min) {
        U <- runif(1, 0, 1)
        A <- Sys_resamp(W = W_0, P = Nparam, U = U)
        sigma_resamp <- sigma_merged[A]
        sigma_array[,2] <- sigma_resamp
        mu_resamp <- mu_merged[A]
        mu_array[,2] <- mu_resamp
        A5_array[,2] <- A5_array[A,1]
        A7_array[,2] <- A7_array[A,1]
        u_pooling_log <- u_pooling_log[A]
        rho_resamp <- rho_0[A]
        theta_0 <- theta_0[,A]
        w_theta_log_0 <- w_theta_log_0[,A]
        w_theta_log_array[,,1] <- w_theta_log_0
        W_theta_0 <- apply(w_theta_log_0, 2, function(x) exp(x - matrixStats::logSumExp(x)))
        W_theta_array[,,1] <- W_theta_0
        w_log_0 <- rep(0, Nparam)
        
        # PMCMC kernel
        inits <- list(mu_0 = mu_resamp,
                      sigma_0 = sigma_resamp,
                      rho_0 = rho_resamp,
                      theta_0 = colMeans(theta_0))
        out_pmcmc <- pmcmc_sub6_kernel_init(data = data[1],
                                            inits = inits,
                                            Nparam = Nparam,
                                            Nstate = Nstate,
                                            Ntotal = Ntotal,
                                            ncore = ncore)
        rho_resamp <- out_pmcmc$rho 
        rho_array[,2] <- rho_resamp
        # phi <- out_pmcmc$phi
        # phi_array[,2] <- phi
        # theta_0 <- out_pmcmc$theta
      } else {
        sigma_resamp <- sigma_merged
        sigma_array[,2] <- sigma_resamp
        mu_resamp <- mu_merged
        mu_array[,2] <- mu_resamp
        A5_array[,2] <- A5_array[,1]
        A7_array[,2] <- A7_array[,1]
        rho_resamp <- rho_0
        rho_array[,2] <- rho_resamp
      }
      
      out_pf <- pf_sub6_update(data = data[1],
                               theta_curr = theta_0,
                               w_theta_log_curr = w_theta_log_0,
                               W_theta_curr = W_theta_0,
                               mu = mu_resamp,
                               sigma = sigma_resamp,
                               rho = rho_resamp,
                               Nparam = Nparam, 
                               Nstate = Nstate,
                               t = t)
      
      theta <- out_pf$theta_next
      theta_array[,,2] <- theta
      w_theta_log <- out_pf$w_theta_log
      w_theta_log_array[,,2] <- w_theta_log
      W_theta <- out_pf$W_theta
      W_theta_array[,,2] <- W_theta
      llike_sum <- out_pf$llike_sum
      theta_array[,,1] <- out_pf$theta_curr
      
      w_log_update <- llike_sum / Nstate + u_pooling_log # At t = 1, we target the true posterior directly from \check{p_13}, so this is alike using alpha = 1 in the annealing process
      w_log <- w_log_0 + w_log_update
      if (max(w_log) < min_lim) w_log <- w_log/1e4
      w_log[which(w_log < min_lim)] <- min_lim
      w_log_array[,2] <- w_log
      W <- exp(w_log - matrixStats::logSumExp(w_log))
      W_array[,2] <- W
      
      theta_inits <- colMeans(theta) # used for the next state node -- store here for convenience
      
      setTxtProgressBar(pb, t)
    }
    
    ESS <- 1 / sum(W^2)
    # rejuvenation
    if (ESS < ESS_min) {
      U <- runif(1, 0, 1)
      A <- Sys_resamp(W = W, P = Nparam, U = U)
      sigma_resamp <- sigma_resamp[A]
      sigma_array[,t+2] <- sigma_resamp
      mu_resamp <- mu_resamp[A]
      mu_array[,t+2] <- mu_resamp
      A5_array[,t+2] <- A5_array[A,t+1]
      A7_array[,t+2] <- A7_array[A,t+1]
      rho_resamp <- rho_resamp[A]
      theta <- theta[,A]
      w_theta_log <- w_theta_log[,A]
      w_theta_log_array[,,t+1] <- w_theta_log
      W_theta <- apply(w_theta_log, 2, function(x) exp(x - matrixStats::logSumExp(x)))
      W_theta_array[,,t+1] <- W_theta
      w_log <- rep(0, Nparam)
      
      if (t == 1) {
        theta_inits <- theta_inits[A]
      } else {
        theta_inits <- theta_inits[,A]
      }
      
      # PMCMC kernel
      inits <- list(mu_0 = mu_resamp,
                    sigma_0 = sigma_resamp,
                    rho_0 = rho_resamp,
                    theta_0 = theta_inits)
      out_pmcmc <- pmcmc_sub6_kernel(data = data[1:(t+1)], 
                                     inits = inits,
                                     Nparam = Nparam,
                                     Nstate = Nstate,
                                     Ntotal = Ntotal,
                                     ncore = ncore,
                                     TT = t + 1,
                                     t = t)
      rho_resamp <- out_pmcmc$rho
      rho_array[,t+2] <- rho_resamp
      # phi <- out_pmcmc$phi
      # phi_array[,t+2] <- phi
      # theta <- out_pmcmc$theta
    } else {
      mu_array[,t+2] <- mu_resamp
      sigma_array[,t+2] <- sigma_resamp
      A5_array[,t+2] <- A5_array[,t+1]
      A7_array[,t+2] <- A7_array[,t+1]
      rho_array[,t+2] <- rho_resamp
      # phi_array[,t+2] <- phi
    }
    
    out_pf <- pf_sub6_update(data = data[t+1], 
                             theta_curr = theta,
                             w_theta_log_curr = w_theta_log,
                             W_theta_curr = W_theta,
                             mu = mu_resamp,
                             sigma = sigma_resamp,
                             rho = rho_resamp,
                             Nparam = Nparam, 
                             Nstate = Nstate,
                             t = t + 1)
    
    theta <- out_pf$theta_next
    theta_array[,,t+2] <- theta
    w_theta_log <- out_pf$w_theta_log
    w_theta_log_array[,,t+2] <- w_theta_log
    W_theta <- out_pf$W_theta
    W_theta_array[,,2] <- W_theta
    llike_sum <- out_pf$llike_sum
    theta_array[,,t+1] <- out_pf$theta_curr
    
    # pooling_prior_log <- pooled_prior_log_sub6(phi56 = sigma_resamp,
    #                                            phi67 = mu_resamp,
    #                                            lambda = lambda) # For t > 1, the target posterior involves the pooled prior
    # w_log_update <- llike_sum / Nstate + pooling_prior_log
    w_log_update <- llike_sum / Nstate
    w_log <- w_log + w_log_update
    if (max(w_log) < min_lim) w_log <- w_log/1e4
    w_log[which(w_log < min_lim)] <- min_lim
    w_log_array[,t+2] <- w_log
    W <- exp(w_log - matrixStats::logSumExp(w_log))
    W_array[,t+2] <- W
    
    theta_inits <- t(apply(theta_array[,,2:(t+2)], 3, colMeans)) # used for the next state node -- store here for convenience
    
    setTxtProgressBar(pb, t + 1)
  }
  
  return(list(mu = mu_array,
              sigma = sigma_array,
              rho = rho_array,
              A5 = A5_array,
              A7 = A7_array,
              W = W_array))
}
dc_melding_smc_sub6 <- function(data,
                                out_sub5,
                                out_sub7,
                                N,
                                n_sub,
                                alpha_incre,
                                Ntotal,
                                n_chains,
                                mu_prior,
                                sigma_prior,
                                alpha_prior,
                                beta_prior,
                                balpha_prior,
                                bbeta_prior,
                                lambda) {
  time_begin <- Sys.time()
  
  ESS_min <- N / 2
  min_lim <- log(.Machine$double.xmin)
  
  TT <- length(data)
  nt <- 1 / alpha_incre
  
  # Storage
  ## Common parameters model storage
  phi_sub_merged_array <- array(rep(NA, N*n_sub*(nt+1)), dim = c(N, n_sub, nt+1))
  phi_sub_index_array <- array(rep(NA, N*n_sub*(nt+1)), dim = c(N, n_sub, nt+1))
  
  ## Unique parameter storage
  rho_array <- array(rep(NA, N*(nt+1)), dim = c(N, nt+1))
  
  ## Weights storage
  w_common_log_array <- array(rep(NA, N*(nt+1)), dim = c(N, nt+1))
  W_common_array <- array(rep(NA, N*(nt+1)), dim = c(N, nt+1))
  
  # Initialization
  ## Initialize particle indices for phi12 and phi23 
  phi_sub_index <- cbind(1:N, 1:N)
  phi_sub_index_merged <- phi_sub_index_array[,,1] <- phi_sub_index
  
  ## Initialize common parameters
  phi56_merged <- out_sub5
  phi67_merged <- out_sub7
  phi_sub_merged_array[,,1] <- cbind(phi56_merged, phi67_merged)
  
  ## Pooled prior
  u_pooling_log <- pooled_6_prior_log(phi56 = phi56_merged,
                                      phi67 = phi67_merged,
                                      mu = mu_prior,
                                      sigma = sigma_prior,
                                      alpha = alpha_prior,
                                      beta = beta_prior,
                                      lambda1 = lambda[1],
                                      lambda2 = lambda[2],
                                      lambda3 = lambda[3])
  
  ## Initialize the weights -- only consider the pooled prior as the weights, see (???) in the paper
  w_common_log_0 <- u_pooling_log 
  if (max(w_common_log_0) < min_lim) w_common_log_0 <- w_common_log_0/1e4
  w_common_log_0[which(w_common_log_0 < min_lim)] <- min_lim
  W_common_0 <- exp(w_common_log_0 - matrixStats::logSumExp(w_common_log_0))
  
  ### Resample -- compulsorily, see Line 1(b) in Algorithm 2 in the paper
  U <- runif(1, 0, 1)
  A <- Sys_resamp(W = W_common_0, P = N, U = U)
  phi56_resamp <- phi56_merged[A]
  phi67_resamp <- phi67_merged[A]
  phi_sub_merged_array[,,1] <- cbind(phi56_resamp, phi67_resamp)
  phi_sub_index_resamp <- phi_sub_index_array[,,1] <- phi_sub_index_merged[A,]
  u_pooling_log <- u_pooling_log[A]
  w_common_log_0 <- w_common_log_array[,1] <- rep(0, N)
  W_common_0 <- W_common_array[,1] <- rep(1, N)
  
  ## Initialize the unique parameter
  rho_0 <- rbeta(N, balpha_prior, bbeta_prior)
  rho_array[,1] <- rho_0
  rho_resamp <- rho_0 # redundant, but for convenience
  x_0 <- matrix(rep(phi67_resamp, TT), nrow = N)
  
  # Tempering SMC
  for (i in 1:(nt-1)) {
    if (i == 1) {
      ## Update alpha
      alpha_update <- alpha_incre
      
      ## Update weights
      p_common_llike <- dnorm(data[1], 0, exp(x_0[,1]/2), log = T)
      for (t in 2:TT) {
        p_common_llike <- p_common_llike + dnorm(data[t], 0, exp(x_0[,t]/2), log = T)
      }
      w_common_log <- w_common_log_0 + alpha_incre * u_pooling_log +
        alpha_incre * p_common_llike
      if (max(w_common_log) < min_lim) w_common_log <- w_common_log/1e4
      w_common_log[which(w_common_log < min_lim)] <- min_lim
      w_common_log_array[,2] <- w_common_log
      W_common <- exp(w_common_log - matrixStats::logSumExp(w_common_log))
      W_common_array[,2] <- W_common
      
      ## Resampling -- optionally
      ESS <- 1 / sum(W_common^2)
      if (ESS < ESS_min) {
        U <- runif(1, 0, 1)
        A <- Sys_resamp(W = W_common, P = N, U = U)
        rho_resamp <- rho_resamp[A]
        x_resamp <- x_0[A,]
        phi56_resamp <- phi56_resamp[A]
        phi67_resamp <- phi67_resamp[A]
        phi_sub_merged_array[,,2] <- cbind(phi56_resamp, phi67_resamp)
        phi_sub_index_resamp <- phi_sub_index_resamp[A,]
        phi_sub_index_array[,,2] <- phi_sub_index_resamp
        u_pooling_log <- u_pooling_log[A]
        w_common_log <- w_common_log_array[,2] <- rep(0, N)
        W_common_array[,2] <- rep(1, N)
      } else {
        rho_resamp <- rho_resamp
        x_resamp <- x_0
        phi56_resamp <- phi56_resamp
        phi67_resamp <- phi67_resamp
        phi_sub_merged_array[,,2] <- cbind(phi56_resamp, phi67_resamp)
        phi_sub_index_resamp <- phi_sub_index_resamp
        phi_sub_index_array[,,2] <- phi_sub_index_resamp
      }
      
      ## MCMC kernel
      inits <- list(rho_inits = rho_resamp,
                    x_inits = x_resamp)
      
      out_jags <- dc_melding_jags_sub6(data = data,
                                       mu = phi67_resamp,
                                       sigma = phi56_resamp,
                                       alpha_prior = balpha_prior,
                                       beta_prior = bbeta_prior,
                                       inits = inits,
                                       alpha_j = alpha_update,
                                       N = N,
                                       Ntotal = Ntotal,
                                       n_chains = n_chains)
      rho <- out_jags$rho
      rho_array[,2] <- rho
      x <- out_jags$x
    }
    
    ## Update alpha
    alpha_update <- alpha_update + alpha_incre
    
    ## Update weights
    p_common_llike <- dnorm(data[1], 0, exp(x[,1]/2), log = T)
    for (t in 2:TT) {
      p_common_llike <- p_common_llike + dnorm(data[t], 0, exp(x[,t]/2), log = T)
    }
    w_common_log <- w_common_log + alpha_incre * u_pooling_log +
      alpha_incre * p_common_llike
    if (max(w_common_log) < min_lim) w_common_log <- w_common_log/1e4
    w_common_log[which(w_common_log < min_lim)] <- min_lim
    w_common_log_array[,i+2] <- w_common_log
    W_common <- exp(w_common_log - matrixStats::logSumExp(w_common_log))
    W_common_array[,i+2] <- W_common
    
    ## Resampling -- optionally
    ESS <- 1 / sum(W_common^2)
    if (ESS < ESS_min) {
      U <- runif(1, 0, 1)
      A <- Sys_resamp(W = W_common, P = N, U = U)
      rho_resamp <- rho[A]
      x_resamp <- x[A,]
      phi56_resamp <- phi56_resamp[A]
      phi67_resamp <- phi67_resamp[A]
      phi_sub_merged_array[,,i+2] <- cbind(phi56_resamp, phi67_resamp)
      phi_sub_index_resamp <- phi_sub_index_resamp[A,]
      phi_sub_index_array[,,i+2] <- phi_sub_index_resamp
      u_pooling_log <- u_pooling_log[A]
      w_common_log <- w_common_log_array[,i+2] <- rep(0, N)
      W_common_array[,i+2] <- rep(1, N)
    } else {
      rho_resamp <- rho
      x_resamp <- x
      phi56_resamp <- phi56_resamp
      phi67_resamp <- phi67_resamp
      phi_sub_merged_array[,,i+2] <- cbind(phi56_resamp, phi67_resamp)
      phi_sub_index_resamp <- phi_sub_index_resamp
      phi_sub_index_array[,,i+2] <- phi_sub_index_resamp
    }
    
    ## MCMC kernel
    inits <- list(rho_inits = rho_resamp,
                  x_inits = x_resamp)
    
    out_jags <- dc_melding_jags_sub6(data = data,
                                     mu = phi67_resamp,
                                     sigma = phi56_resamp,
                                     alpha_prior = balpha_prior,
                                     beta_prior = bbeta_prior,
                                     inits = inits,
                                     alpha_j = alpha_update,
                                     N = N,
                                     Ntotal = Ntotal,
                                     n_chains = n_chains)
    rho <- out_jags$rho
    rho_array[,i+2] <- rho
    x <- out_jags$x
  }
  
  time_end <- Sys.time()
  running_time <- difftime(time_end, time_begin)
  
  return(list(psi6 = rho_array,
              phi = phi_sub_merged_array,
              phi_index = phi_sub_index_array,
              W = W_common_array,
              running_time = running_time))
}








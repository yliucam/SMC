dc_melding_smc_sub2 <- function(data,
                                out_sub1,
                                out_sub3,
                                N,
                                n_sub,
                                alpha_incre,
                                Ntotal,
                                n_chains,
                                mu_prior,
                                sigma_prior,
                                alpha_prior,
                                beta_prior,
                                nu_p_prior,
                                lambda) {
  time_begin <- Sys.time()
  
  ESS_min <- N / 2
  min_lim <- log(.Machine$double.xmin)
  
  nt <- 1 / alpha_incre
  
  # Storage
  ## Common parameters model storage
  phi_sub_merged_array <- array(rep(NA, N*n_sub*(nt+1)), dim = c(N, n_sub, nt+1))
  phi_sub_index_array <- array(rep(NA, N*n_sub*(nt+1)), dim = c(N, n_sub, nt+1))
  
  ## Unique parameter storage
  nu_array <- array(rep(NA, N*(nt+1)), dim = c(N, nt+1))
  
  ## Weights storage
  w_common_log_array <- array(rep(NA, N*(nt+1)), dim = c(N, nt+1))
  W_common_array <- array(rep(NA, N*(nt+1)), dim = c(N, nt+1))
  
  # Initialization
  ## Initialize particle indices for phi12 and phi23 
  phi_sub_index <- cbind(1:N, 1:N)
  phi_sub_index_merged <- phi_sub_index_array[,,1] <- phi_sub_index
  
  ## Initialize common parameters
  phi12_merged <- out_sub1
  phi23_merged <- out_sub3
  phi_sub_merged_array[,,1] <- cbind(phi12_merged, phi23_merged)
  
  ## Pooled prior
  u_pooling_log <- pooled_2_prior_log(phi12 = phi12_merged,
                                      phi23 = phi23_merged,
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
  phi12_resamp <- phi12_merged[A]
  phi23_resamp <- phi23_merged[A]
  phi_sub_merged_array[,,1] <- cbind(phi12_resamp, phi23_resamp)
  phi_sub_index_resamp <- phi_sub_index_array[,,1] <- phi_sub_index_merged[A,]
  u_pooling_log <- u_pooling_log[A]
  w_common_log_0 <- w_common_log_array[,1] <- rep(0, N)
  W_common_0 <- W_common_array[,1] <- rep(1, N)
  
  ## Initialize the unique parameter
  nu_0 <- rcat(N, nu_p_prior)
  nu_array[,1] <- nu_0
  nu_resamp <- nu_0 # redundant, but for convenience
  
  # Tempering SMC
  for (i in 1:(nt-1)) {
    if (i == 1) {
      ### Update alpha
      alpha_update <- alpha_incre
      
      ### Update weights
      p_common_llike <- tstudent_log_particle(data = data,
                                              mu = phi12_resamp,
                                              tau = phi23_resamp,
                                              nu = nu_resamp)
      w_common_log <- w_common_log_0 + alpha_incre * u_pooling_log +
        alpha_incre * p_common_llike
      if (max(w_common_log) < min_lim) w_common_log <- w_common_log/1e4
      w_common_log[which(w_common_log < min_lim)] <- min_lim
      w_common_log_array[,2] <- w_common_log
      W_common <- exp(w_common_log - matrixStats::logSumExp(w_common_log))
      W_common_array[,2] <- W_common
      
      ### Resampling -- optionally
      ESS <- 1 / sum(W_common^2)
      if (ESS < ESS_min) {
        U <- runif(1, 0, 1)
        A <- Sys_resamp(W = W_common, P = N, U = U)
        nu_resamp <- nu_resamp[A]
        phi12_resamp <- phi12_resamp[A]
        phi23_resamp <- phi23_resamp[A]
        phi_sub_merged_array[,,2] <- cbind(phi12_resamp, phi23_resamp)
        phi_sub_index_resamp <- phi_sub_index_resamp[A,]
        phi_sub_index_array[,,2] <- phi_sub_index_resamp
        u_pooling_log <- u_pooling_log[A]
        w_common_log <- w_common_log_array[,2] <- rep(0, N)
        W_common_array[,2] <- rep(1, N)
      } else {
        nu_resamp <- nu_resamp
        phi12_resamp <- phi12_resamp
        phi23_resamp <- phi23_resamp
        phi_sub_merged_array[,,2] <- cbind(phi12_resamp, phi23_resamp)
        phi_sub_index_resamp <- phi_sub_index_resamp
        phi_sub_index_array[,,2] <- phi_sub_index_resamp
      }
      
      ### MCMC kernel
      inits <- list(nu_inits = nu_resamp)
      out_jags <- dc_melding_jags_sub2(data = data,
                                       mu = phi12_resamp,
                                       tau = phi23_resamp,
                                       nu_p = nu_p_prior,
                                       alpha_j = alpha_update,
                                       inits = inits,
                                       N = N,
                                       Ntotal = Ntotal,
                                       n_chains = n_chains)
      nu <- out_jags$nu
      nu_array[,2] <- nu
    }
    
    ### Update alpha
    alpha_update <- alpha_update + alpha_incre
    
    ### Update weights
    p_common_llike <- tstudent_log_particle(data = data,
                                            mu = phi12_resamp,
                                            tau = phi23_resamp,
                                            nu = nu)
    w_common_log <- w_common_log + alpha_incre * u_pooling_log +
      alpha_incre * p_common_llike
    if (max(w_common_log) < min_lim) w_common_log <- w_common_log/1e4
    w_common_log[which(w_common_log < min_lim)] <- min_lim
    w_common_log_array[,i+2] <- w_common_log
    W_common <- exp(w_common_log - matrixStats::logSumExp(w_common_log))
    W_common_array[,i+2] <- W_common
    
    ### Resampling -- optionally
    ESS <- 1 / sum(W_common^2)
    if (ESS < ESS_min) {
      U <- runif(1, 0, 1)
      A <- Sys_resamp(W = W_common, P = N, U = U)
      nu_resamp <- nu[A]
      phi12_resamp <- phi12_resamp[A]
      phi23_resamp <- phi23_resamp[A]
      phi_sub_merged_array[,,i+2] <- cbind(phi12_resamp, phi23_resamp)
      phi_sub_index_resamp <- phi_sub_index_resamp[A,]
      phi_sub_index_array[,,i+2] <- phi_sub_index_resamp
      u_pooling_log <- u_pooling_log[A]
      w_common_log <- w_common_log_array[,i+2] <- rep(0, N)
      W_common_array[,i+2] <- rep(1, N)
    } else {
      nu_resamp <- nu
      phi12_resamp <- phi12_resamp
      phi23_resamp <- phi23_resamp
      phi_sub_merged_array[,,i+2] <- cbind(phi12_resamp, phi23_resamp)
      phi_sub_index_resamp <- phi_sub_index_resamp
      phi_sub_index_array[,,i+2] <- phi_sub_index_resamp
    }
    
    ### MCMC kernel
    inits <- list(nu_inits = nu_resamp)
    out_jags <- dc_melding_jags_sub2(data = data,
                                     mu = phi12_resamp,
                                     tau = phi23_resamp,
                                     nu_p = nu_p_prior,
                                     alpha_j = alpha_update,
                                     inits = inits,
                                     N = N,
                                     Ntotal = Ntotal,
                                     n_chains = n_chains)
    nu <- out_jags$nu
    nu_array[,i+2] <- nu
  }
  
  time_end <- Sys.time()
  running_time <- difftime(time_end, time_begin)
  
  return(list(psi2 = nu_array,
              phi = phi_sub_merged_array,
              phi_index = phi_sub_index_array,
              W = W_common_array,
              running_time = running_time))
}









dc_melding_smc_sub10 <- function(data,
                                 out_sub9,
                                 out_sub11,
                                 N,
                                 n_sub,
                                 alpha_incre,
                                 Ntotal,
                                 n_chains,
                                 mu_prior,
                                 sigma_prior,
                                 alpha_prior,
                                 beta_prior,
                                 nu_p_prior,
                                 lambda) {
  time_begin <- Sys.time()
  
  ESS_min <- N / 2
  min_lim <- log(.Machine$double.xmin)
  
  nt <- 1 / alpha_incre
  
  # Storage
  ## Common parameters model storage
  phi_sub_merged_array <- array(rep(NA, N*n_sub*(nt+1)), dim = c(N, n_sub, nt+1))
  phi_sub_index_array <- array(rep(NA, N*n_sub*(nt+1)), dim = c(N, n_sub, nt+1))
  
  ## Unique parameter storage
  nu_array <- array(rep(NA, N*(nt+1)), dim = c(N, nt+1))
  
  ## Weights storage
  w_common_log_array <- array(rep(NA, N*(nt+1)), dim = c(N, nt+1))
  W_common_array <- array(rep(NA, N*(nt+1)), dim = c(N, nt+1))
  
  # Initialization
  ## Initialize particle indices for phi9_10 and phi10_11 
  phi_sub_index <- cbind(1:N, 1:N)
  phi_sub_index_merged <- phi_sub_index_array[,,1] <- phi_sub_index
  
  ## Initialize common parameters
  phi9_10_merged <- out_sub9
  phi10_11_merged <- out_sub11
  phi_sub_merged_array[,,1] <- cbind(phi9_10_merged, phi10_11_merged)
  
  ## Pooled prior
  u_pooling_log <- pooled_10_prior_log(phi9_10 = phi9_10_merged,
                                       phi10_11 = phi10_11_merged,
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
  phi9_10_resamp <- phi9_10_merged[A]
  phi10_11_resamp <- phi10_11_merged[A]
  phi_sub_merged_array[,,1] <- cbind(phi9_10_resamp, phi10_11_resamp)
  phi_sub_index_resamp <- phi_sub_index_array[,,1] <- phi_sub_index_merged[A,]
  u_pooling_log <- u_pooling_log[A]
  w_common_log_0 <- w_common_log_array[,1] <- rep(0, N)
  W_common_0 <- W_common_array[,1] <- rep(1, N)
  
  ## Initialize the unique parameter
  nu_0 <- rcat(N, nu_p_prior)
  nu_array[,1] <- nu_0
  nu_resamp <- nu_0 # redundant but for convenience
  
  # Tempering SMC
  for (i in 1:(nt-1)) {
    if (i == 1) {
      ### Update alpha
      alpha_update <- alpha_incre
      
      ### Update weights
      p_common_llike <- tstudent_log_particle(data = data,
                                              mu = phi9_10_resamp,
                                              tau = phi10_11_resamp,
                                              nu = nu_resamp)
      w_common_log <- w_common_log_0 + alpha_incre * u_pooling_log +
        alpha_incre * p_common_llike
      if (max(w_common_log) < min_lim) w_common_log <- w_common_log/1e4
      w_common_log[which(w_common_log < min_lim)] <- min_lim
      w_common_log_array[,2] <- w_common_log
      W_common <- exp(w_common_log - matrixStats::logSumExp(w_common_log))
      W_common_array[,2] <- W_common
      
      ### Resampling -- optionally
      ESS <- 1 / sum(W_common^2)
      if (ESS < ESS_min) {
        U <- runif(1, 0, 1)
        A <- Sys_resamp(W = W_common, P = N, U = U)
        nu_resamp <- nu_resamp[A]
        phi9_10_resamp <- phi9_10_resamp[A]
        phi10_11_resamp <- phi10_11_resamp[A]
        phi_sub_merged_array[,,2] <- cbind(phi9_10_resamp, phi10_11_resamp)
        phi_sub_index_resamp <- phi_sub_index_resamp[A,]
        phi_sub_index_array[,,2] <- phi_sub_index_resamp
        u_pooling_log <- u_pooling_log[A]
        w_common_log <- w_common_log_array[,2] <- rep(0, N)
        W_common_array[,2] <- rep(1, N)
      } else {
        nu_resamp <- nu_resamp
        phi9_10_resamp <- phi9_10_resamp
        phi10_11_resamp <- phi10_11_resamp
        phi_sub_merged_array[,,2] <- cbind(phi9_10_resamp, phi10_11_resamp)
        phi_sub_index_resamp <- phi_sub_index_resamp
        phi_sub_index_array[,,2] <- phi_sub_index_resamp
      }
      
      ### MCMC kernel
      inits <- list(nu_inits = nu_resamp)
      out_jags <- dc_melding_jags_sub2(data = data,
                                       mu = phi9_10_resamp,
                                       tau = phi10_11_resamp,
                                       nu_p = nu_p_prior,
                                       alpha_j = alpha_update,
                                       inits = inits,
                                       N = N,
                                       Ntotal = Ntotal,
                                       n_chains = n_chains)
      nu <- out_jags$nu
      nu_array[,2] <- nu
    }
    
    ### Update alpha
    alpha_update <- alpha_update + alpha_incre
    
    ### Update weights
    p_common_llike <- tstudent_log_particle(data = data,
                                            mu = phi9_10_resamp,
                                            tau = phi10_11_resamp,
                                            nu = nu)
    w_common_log <- w_common_log + alpha_incre * u_pooling_log +
      alpha_incre * p_common_llike
    if (max(w_common_log) < min_lim) w_common_log <- w_common_log/1e4
    w_common_log[which(w_common_log < min_lim)] <- min_lim
    w_common_log_array[,i+2] <- w_common_log
    W_common <- exp(w_common_log - matrixStats::logSumExp(w_common_log))
    W_common_array[,i+2] <- W_common
    
    ### Resampling -- optionally
    ESS <- 1 / sum(W_common^2)
    if (ESS < ESS_min) {
      U <- runif(1, 0, 1)
      A <- Sys_resamp(W = W_common, P = N, U = U)
      nu_resamp <- nu[A]
      phi9_10_resamp <- phi9_10_resamp[A]
      phi10_11_resamp <- phi10_11_resamp[A]
      phi_sub_merged_array[,,i+2] <- cbind(phi9_10_resamp, phi10_11_resamp)
      phi_sub_index_resamp <- phi_sub_index_resamp[A,]
      phi_sub_index_array[,,i+2] <- phi_sub_index_resamp
      u_pooling_log <- u_pooling_log[A]
      w_common_log <- w_common_log_array[,i+2] <- rep(0, N)
      W_common_array[,i+2] <- rep(1, N)
    } else {
      nu_resamp <- nu
      phi9_10_resamp <- phi9_10_resamp
      phi10_11_resamp <- phi10_11_resamp
      phi_sub_merged_array[,,i+2] <- cbind(phi9_10_resamp, phi10_11_resamp)
      phi_sub_index_resamp <- phi_sub_index_resamp
      phi_sub_index_array[,,i+2] <- phi_sub_index_resamp
    }
    
    ### MCMC kernel
    inits <- list(nu_inits = nu_resamp)
    out_jags <- dc_melding_jags_sub2(data = data,
                                     mu = phi9_10_resamp,
                                     tau = phi10_11_resamp,
                                     nu_p = nu_p_prior,
                                     alpha_j = alpha_update,
                                     inits = inits,
                                     N = N,
                                     Ntotal = Ntotal,
                                     n_chains = n_chains)
    nu <- out_jags$nu
    nu_array[,i+2] <- nu
  }
  
  time_end <- Sys.time()
  running_time <- difftime(time_end, time_begin)
  
  return(list(psi2 = nu_array,
              phi = phi_sub_merged_array,
              phi_index = phi_sub_index_array,
              W = W_common_array,
              running_time = running_time))
}



















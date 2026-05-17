dc_melding_smc_sub4 <- function(data,
                                out_sub3,
                                out_sub5,
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
  phi34_merged <- out_sub3
  phi45_merged <- out_sub5
  phi_sub_merged_array[,,1] <- cbind(phi34_merged, phi45_merged)
  
  ## Pooled prior
  u_pooling_log <- pooled_4_prior_log(phi34 = phi34_merged,
                                      phi45 = phi45_merged,
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
  phi34_resamp <- phi34_merged[A]
  phi45_resamp <- phi45_merged[A]
  phi_sub_merged_array[,,1] <- cbind(phi34_resamp, phi45_resamp)
  phi_sub_index_resamp <- phi_sub_index_array[,,1] <- phi_sub_index_merged[A,]
  u_pooling_log <- u_pooling_log[A]
  w_common_log_0 <- w_common_log_array[,1] <- rep(0, N)
  W_common_0 <- W_common_array[,1] <- rep(1, N)
  
  ## Initialize the unique parameter
  rho_0 <- rbeta(N, balpha_prior, bbeta_prior)
  rho_array[,1] <- rho_0
  rho_resamp <- rho_0 # redundant, but for convenience
  x_0 <- matrix(rep(phi34_resamp, TT), nrow = N)
  
  # Tempering SMC
  for (i in 1:(nt-1)) {
    if (i == 1) {
      ## Update alpha
      alpha_update <- alpha_incre
      
      ## Update weights
      p_common_llike <- dnorm(data[1], x_0[,1], 1, log = T)
      for (t in 2:TT) {
        p_common_llike <- p_common_llike + dnorm(data[t], x_0[,t], 1, log = T)
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
        phi34_resamp <- phi34_resamp[A]
        phi45_resamp <- phi45_resamp[A]
        phi_sub_merged_array[,,2] <- cbind(phi34_resamp, phi45_resamp)
        phi_sub_index_resamp <- phi_sub_index_resamp[A,]
        phi_sub_index_array[,,2] <- phi_sub_index_resamp
        u_pooling_log <- u_pooling_log[A]
        w_common_log <- w_common_log_array[,2] <- rep(0, N)
        W_common_array[,2] <- rep(1, N)
      } else {
        rho_resamp <- rho_resamp
        x_resamp <- x_0
        phi34_resamp <- phi34_resamp
        phi45_resamp <- phi45_resamp
        phi_sub_merged_array[,,2] <- cbind(phi34_resamp, phi45_resamp)
        phi_sub_index_resamp <- phi_sub_index_resamp
        phi_sub_index_array[,,2] <- phi_sub_index_resamp
      }
      
      ## MCMC kernel
      inits <- list(rho_inits = rho_resamp,
                    x_inits = x_resamp)
      out_jags <- dc_melding_jags_sub4(data = data,
                                       mu = phi34_resamp,
                                       sigma = phi45_resamp,
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
    p_common_llike <- dnorm(data[1], x[,1], 1, log = T)
    for (t in 2:TT) {
      p_common_llike <- p_common_llike + dnorm(data[t], x[,t], 1, log = T)
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
      phi34_resamp <- phi34_resamp[A]
      phi45_resamp <- phi45_resamp[A]
      phi_sub_merged_array[,,i+2] <- cbind(phi34_resamp, phi45_resamp)
      phi_sub_index_resamp <- phi_sub_index_resamp[A,]
      phi_sub_index_array[,,i+2] <- phi_sub_index_resamp
      u_pooling_log <- u_pooling_log[A]
      w_common_log <- w_common_log_array[,i+2] <- rep(0, N)
      W_common_array[,i+2] <- rep(1, N)
    } else {
      rho_resamp <- rho
      x_resamp <- x
      phi34_resamp <- phi34_resamp
      phi45_resamp <- phi45_resamp
      phi_sub_merged_array[,,i+2] <- cbind(phi34_resamp, phi45_resamp)
      phi_sub_index_resamp <- phi_sub_index_resamp
      phi_sub_index_array[,,i+2] <- phi_sub_index_resamp
    }
    
    ## MCMC kernel
    inits <- list(rho_inits = rho_resamp,
                  x_inits = x_resamp)
    out_jags <- dc_melding_jags_sub4(data = data,
                                     mu = phi34_resamp,
                                     sigma = phi45_resamp,
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
  
  return(list(psi4 = rho_array,
              phi = phi_sub_merged_array,
              phi_index = phi_sub_index_array,
              W = W_common_array,
              running_time = running_time))
}










dc_melding_smc_sub8 <- function(data,
                                out_sub7,
                                out_sub9,
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
  phi78_merged <- out_sub7
  phi89_merged <- out_sub9
  phi_sub_merged_array[,,1] <- cbind(phi78_merged, phi89_merged)
  
  ## Pooled prior
  u_pooling_log <- pooled_8_prior_log(phi78 = phi78_merged,
                                      phi89 = phi89_merged,
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
  phi78_resamp <- phi78_merged[A]
  phi89_resamp <- phi89_merged[A]
  phi_sub_merged_array[,,1] <- cbind(phi78_resamp, phi89_resamp)
  phi_sub_index_resamp <- phi_sub_index_array[,,1] <- phi_sub_index_merged[A,]
  u_pooling_log <- u_pooling_log[A]
  w_common_log_0 <- w_common_log_array[,1] <- rep(0, N)
  W_common_0 <- W_common_array[,1] <- rep(1, N)
  
  ## Initialize the unique parameter
  rho_0 <- rbeta(N, balpha_prior, bbeta_prior)
  rho_array[,1] <- rho_0
  rho_resamp <- rho_0 # redundant, but for convenience
  # x_0 <- matrix(rep(phi78_resamp, TT), nrow = N)
  x_0 <- matrix(rep(4, N*TT), nrow = N)
  
  # Tempering SMC
  for (i in 1:(nt-1)) {
    if (i == 1) {
      ## Update alpha
      alpha_update <- alpha_incre
      
      ## Update weights
      p_common_llike <- dnorm(data[1], x_0[,1], 1, log = T)
      for (t in 2:TT) {
        p_common_llike <- p_common_llike + dnorm(data[t], x_0[,t], 1, log = T)
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
        phi78_resamp <- phi78_resamp[A]
        phi89_resamp <- phi89_resamp[A]
        phi_sub_merged_array[,,2] <- cbind(phi78_resamp, phi89_resamp)
        phi_sub_index_resamp <- phi_sub_index_resamp[A,]
        phi_sub_index_array[,,2] <- phi_sub_index_resamp
        u_pooling_log <- u_pooling_log[A]
        w_common_log <- w_common_log_array[,2] <- rep(0, N)
        W_common_array[,2] <- rep(1, N)
      } else {
        rho_resamp <- rho_resamp
        x_resamp <- x_0
        phi78_resamp <- phi78_resamp
        phi89_resamp <- phi89_resamp
        phi_sub_merged_array[,,2] <- cbind(phi78_resamp, phi89_resamp)
        phi_sub_index_resamp <- phi_sub_index_resamp
        phi_sub_index_array[,,2] <- phi_sub_index_resamp
      }
      
      ## MCMC kernel
      inits <- list(rho_inits = rho_resamp,
                    x_inits = x_resamp)
      out_jags <- dc_melding_jags_sub4(data = data,
                                       mu = phi78_resamp,
                                       sigma = phi89_resamp,
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
    p_common_llike <- dnorm(data[1], x[,1], 1, log = T)
    for (t in 2:TT) {
      p_common_llike <- p_common_llike + dnorm(data[t], x[,t], 1, log = T)
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
      phi78_resamp <- phi78_resamp[A]
      phi89_resamp <- phi89_resamp[A]
      phi_sub_merged_array[,,i+2] <- cbind(phi78_resamp, phi89_resamp)
      phi_sub_index_resamp <- phi_sub_index_resamp[A,]
      phi_sub_index_array[,,i+2] <- phi_sub_index_resamp
      u_pooling_log <- u_pooling_log[A]
      w_common_log <- w_common_log_array[,i+2] <- rep(0, N)
      W_common_array[,i+2] <- rep(1, N)
    } else {
      rho_resamp <- rho
      x_resamp <- x
      phi78_resamp <- phi78_resamp
      phi89_resamp <- phi89_resamp
      phi_sub_merged_array[,,i+2] <- cbind(phi78_resamp, phi89_resamp)
      phi_sub_index_resamp <- phi_sub_index_resamp
      phi_sub_index_array[,,i+2] <- phi_sub_index_resamp
    }
    
    ## MCMC kernel
    inits <- list(rho_inits = rho_resamp,
                  x_inits = x_resamp)
    out_jags <- dc_melding_jags_sub4(data = data,
                                     mu = phi78_resamp,
                                     sigma = phi89_resamp,
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
  
  return(list(psi8 = rho_array,
              phi = phi_sub_merged_array,
              phi_index = phi_sub_index_array,
              W = W_common_array,
              running_time = running_time))
}









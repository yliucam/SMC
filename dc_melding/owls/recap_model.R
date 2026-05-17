recap_model <- function(data,
                        N,
                        Ntotal,
                        n_chains,
                        alpha_incre) {
  ESS_min <- N / 2
  min_lim <- log(.Machine$double.xmin)
  
  data_f_j <- data$recap_f_j 
  data_f_a <- data$recap_f_a 
  data_m_j <- data$recap_m_j 
  data_m_a <- data$recap_m_a
  
  ti <- dim(data_f_j)[2]
  nt <- 1 / alpha
  
  # Storage
  alpha0_array <- alpha2_array <- array(rep(NA, (nt+1)*N), dim = c(N, nt+1))
  alpha1_array <- alpha4_array <- array(rep(NA, (nt+1)*N), dim = c(N, nt+1))
  alpha5_array <- array(rep(NA, (nt+1)*(ti-1)*N), dim = c(N, ti-1, nt+1))
  w_log_array <- array(rep(NA, (nt+1)*N), dim = c(N, nt+1))
  W_array <- array(rep(NA, (nt+1)*N), dim = c(N, nt+1))
  
  # Initialization
  alpha0_0 <- rtruncnorm(N, -10, 10, 0, 2)
  alpha2_0 <- rtruncnorm(N, -10, 10, 0, 2)
  alpha1_0 <- rtruncnorm(N, -10, 10, 0, 2)
  alpha4_0 <- rtruncnorm(N, -10, 10, 0, 2)
  alpha5_0 <- alpha5_resamp <- matrix(nrow = N, ncol = ti-1)
  for (u in 1:(ti-1)) {
    alpha5_0[,u] ~ rtruncnorm(N, -10, 10, 0, 2)
  }
  
  w_log_0 <- dtruncnorm(alpha0_0, -10, 10, 0, 2) + dtruncnorm(alpha2_0, -10, 10, 0, 2) 
    dtruncnorm(alpha1_0, -10, 10, 0, 2) + dtruncnorm(alpha4_0, -10, 10, 0, 2)
  for (u in 1:(ti-1)) {
    w_log_0 <- w_log_0 + dtruncnorm(alpha5_0[,u], -10, 10, 0, 2)
  }
  w_log_array[,1] <- w_log_0
  W_0 <- exp(w_log_0 - matrixStats::logSumExp(w_log_0))
  W_array[,1] <- W_0
  
  for (i in 1:(nt-1)) {
    if (i == 1) {
      # Resampling -- optionally
      ESS <- 1 / sum(W_0^2)
      if (ESS < ESS_min) {
        U <- runif(1, 0, 1)
        A <- Sys_resamp(W = W_0, P = N, U = U)
        alpha0_resamp <- alpha0_0[A]
        alpha2_resamp <- alpha2_0[A]
        alpha1_resamp <- alpha1_0[A]
        alpha4_resamp <- alpha4_0[A]
        for (u in 1:(ti-1)) {
          alpha5_resamp[,u] <- alpha5_0[A,u]
        }
        w_log_0 <- w_log_array[,1] <- rep(0, N)
        W_0 <- W_array[,1] <- rep(1, N)
      } else {
        alpha0_resamp <- alpha0_0
        alpha2_resamp <- alpha2_0
        alpha1_resamp <- alpha1_0
        alpha4_resamp <- alpha4_0
        alpha5_resamp <- alpha5_0
      }
      
      ## Update alpha
      alpha_update <- alpha_incre
      
      ## MCMC kernel
      inits <- list(alpha0_0 = alpha0_resamp,
                    alpha2_0 = alpha2_resamp,
                    alpha1_0 = alpha1_resamp,
                    alpha4_0 = alpha4_resamp,
                    alpha5_0 = alpha5_resamp)
      out_jags <- recap_jags(data = data,
                             alpha_j = alpha_incre,
                             inits = inits,
                             N = N,
                             Ntotal = Ntotal,
                             n_chains = n_chains)
      alpha0 <- out_jags$alpha0_res
      alpha0_array[,2] <- alpha0
      alpha2 <- out_jags$alpha2_res
      alpha2_array[,2] <- alpha2
      alpha1 <- out_jags$alpha1_res
      
    }
  }
}




count_model <- function(data2,
                        N,
                        alpha_incre,
                        Ntotal,
                        n_chains) {
  tt <- length(data2)
  nt <- 1/alpha_incre
  
  ESS_min <- N / 2
  min_lim <- log(.Machine$double.xmin)
  
  alpha0_array <- alpha2_array <- rho_array <- array(rep(NA, N*(nt+1)), dim = c(N, nt+1))
  alpha6_array <- array(rep(NA, N*(nt+1)), dim = c(N, nt+1))
  xJ_array <- array(rep(NA, tt*N*(nt+1)), dim = c(N, tt, nt+1))
  sur_array <- array(rep(NA, tt*N*(nt+1)), dim = c(N, tt, nt+1))
  imm_array <- array(rep(NA, tt*N*(nt+1)), dim = c(N, tt, nt+1))
  w_log_array <- array(rep(NA, N*(nt+1)), dim = c(N, nt+1))
  W_array <- array(rep(NA, N*(nt+1)), dim = c(N, nt+1))
  
  alpha0_0 <- rtruncnorm(N, -10, 10, 0, 2)
  alpha2_0 <- rtruncnorm(N, -10, 10, 0, 2)
  rho_0 <- rep(5, N)
  alpha6_0 <- rtruncnorm(N, -10, 10, 0, 2)
  # We use the same initializations for the latent variables at every tempering node
  xJ_0 <- matrix(rep(10, tt*N), nrow = N)
  sur_0 <- matrix(rep(10, tt*N), nrow = N)
  imm_0 <- matrix(rep(10, tt*N), nrow = N)
  
  alpha0_array[,1] <- alpha0_0
  alpha2_array[,1] <- alpha2_0
  rho_array[,1] <- rho_0
  alpha6_array[,1] <- alpha6_0
  xJ_array[,,1] <- xJ_0
  sur_array[,,1] <- sur_0
  imm_array[,,1] <- imm_0
  
  w_log_0 <- dtruncnorm(alpha0_0, -10, 10, 0, 2) + dtruncnorm(alpha2_0, -10, 10, 0, 2) +
    dtruncnorm(alpha6_0, -10, 10, 0, 2)
  w_log_array[,1] <- w_log_0
  W_0 <- exp(w_log_0 - matrixStats::logSumExp(w_log_0))
  W_array[,1] <- W_0
  
  pb <- txtProgressBar(min = 1, max = nt, style = 3)
  
  for (i in 1:(nt-1)) {
    if (i == 1) {
      ## Resampling -- optionally
      ESS <- 1 / sum(W_0^2)
      if (ESS < ESS_min) {
        U <- runif(1, 0, 1)
        A <- Sys_resamp(W = W_0, P = N, U = U)
        alpha0_resamp <- alpha0_0[A]
        alpha2_resamp <- alpha2_0[A]
        alpha6_resamp <- alpha6_0[A]
        rho_resamp <- rho_0
        w_log_0 <- w_log_array[,1] <- rep(0, N)
        W_0 <- W_array[,1] <- rep(1, N)
      } else {
        alpha0_resamp <- alpha0_0
        alpha2_resamp <- alpha2_0
        alpha6_resamp <- alpha6_0
        rho_resamp <- rho_0
      }
      
      ## Update alpha
      alpha_update <- alpha_incre
      
      ## MCMC kernel
      inits <- list(alpha0_0 = alpha0_resamp,
                    alpha2_0 = alpha2_resamp,
                    rho_0 = rho_resamp,
                    alpha6_0 = alpha6_resamp,
                    xJ_0 = xJ_0,
                    sur_0 = sur_0,
                    imm_0 = imm_0)
      out_jags <- count_jags(data = data2,
                             alpha_j = alpha_update,
                             inits = inits,
                             N = N,
                             Ntotal = Ntotal,
                             n_chains = n_chains)
      alpha0 <- out_jags$alpha0_res
      alpha0_array[,2] <- alpha0
      alpha2 <- out_jags$alpha2_res
      alpha2_array[,2] <- alpha2
      rho <- out_jags$rho
      rho_array[,2] <- rho
      alpha6 <- out_jags$alpha6_res
      alpha6_array[,2] <- alpha6
      xJ <- out_jags$xJ_res
      xJ_array[,,2] <- xJ
      sur <- out_jags$sur_res
      sur_array[,,2] <- sur
      imm <- out_jags$imm_res
      imm_array[,,2] <- imm
      
      ## Update weights
      llike <- p_common_likelihood(data = data2,
                                   xJ = xJ,
                                   sur = sur,
                                   imm = imm)
      w_log <- w_log_0 + alpha_incre * llike
      if (max(w_log) < min_lim) w_log <- w_log/1e4
      w_log[which(w_log < min_lim)] <- min_lim
      w_log_array[,2] <- w_log
      W <- exp(w_log - matrixStats::logSumExp(w_log))
      W_array[,2] <- W
      
      setTxtProgressBar(pb, i)
    }
    
    ## Resampling -- optionally
    ESS <- 1 / sum(W^2)
    if (ESS < ESS_min) {
      U <- runif(1, 0, 1)
      A <- Sys_resamp(W = W, P = N, U = U)
      alpha0_resamp <- alpha0[A]
      alpha2_resamp <- alpha2[A]
      alpha6_resamp <- alpha6[A]
      rho_resamp <- rho[A]
      alpha0_array[,i+1] <- alpha0_resamp
      alpha2_array[,i+1] <- alpha2_resamp
      alpha6_array[,i+1] <- alpha6_resamp
      rho_array[,i+1] <- rho_resamp
      w_log <- w_log_array[,i+1] <- rep(0, N)
      W <- W_array[,i+1] <- rep(1, N)
    } else {
      alpha0_resamp <- alpha0
      alpha2_resamp <- alpha2
      alpha6_resamp <- alpha6
      rho_resamp <- rho
    }
    
    ## Update alpha
    alpha_update <- alpha_update + alpha_incre
    
    ## MCMC kernel
    inits <- list(alpha0_0 = alpha0_resamp,
                  alpha2_0 = alpha2_resamp,
                  rho_0 = rho_resamp,
                  alpha6_0 = alpha6_resamp,
                  xJ_0 = xJ_0,
                  sur_0 = sur_0,
                  imm_0 = imm_0)
    out_jags <- count_jags(data = data2,
                           alpha_j = alpha_update,
                           inits = inits,
                           N = N,
                           Ntotal = Ntotal,
                           n_chains = n_chains)
    alpha0 <- out_jags$alpha0_res
    alpha0_array[,i+2] <- alpha0
    alpha2 <- out_jags$alpha2_res
    alpha2_array[,i+2] <- alpha2
    rho <- out_jags$rho
    rho_array[,i+2] <- rho
    alpha6 <- out_jags$alpha6_res
    alpha6_array[,i+2] <- alpha6
    xJ <- out_jags$xJ_res
    xJ_array[,,i+2] <- xJ
    sur <- out_jags$sur_res
    sur_array[,,i+2] <- sur
    imm <- out_jags$imm_res
    imm_array[,,i+2] <- imm
    
    ## Update weights
    llike <- p_common_likelihood(data = data2,
                                 xJ = xJ,
                                 sur = sur,
                                 imm = imm)
    w_log <- w_log + alpha_incre * llike
    if (max(w_log) < min_lim) w_log <- w_log/1e4
    w_log[which(w_log < min_lim)] <- min_lim
    w_log_array[,i+2] <- w_log
    W <- exp(w_log - matrixStats::logSumExp(w_log))
    W_array[,i+2] <- W
    
    setTxtProgressBar(pb, i+1)
  }
  
  return(list(alpha0 = alpha0_array,
              alpha2 = alpha2_array,
              rho = rho_array,
              alpha6 = alpha6_array,
              xJ = xJ_array,
              sur = sur_array,
              imm = imm_array,
              W = W_array))
}





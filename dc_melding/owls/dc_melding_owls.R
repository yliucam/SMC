# Systematic resampling - Algorithm 9.6 in Chopin & Papaspiliopoulous (2020)

Sys_resamp <- function(W, P, U) {
  N <- length(W)
  A <- rep(NA, P)
  v <- N * cumsum(W)
  s <- U
  m <- 1
  for (i in 1:P) {
    while (v[m] < s) {
      m <- m + 1
      #s <- s + 1
    }
    A[i] <- m
    s <- s + 1
  }
  
  return(A)
}



dc_melding_owls_common <- function(data2,
                                   N,
                                   n_sub_var,
                                   alpha_incre,
                                   alpha0,
                                   alpha2,
                                   rho,
                                   lambda,
                                   Ntotal,
                                   n_chains) {
  tt <- length(data2)
  nt <- 1/alpha_incre
  
  ESS_min <- N / 2
  min_lim <- log(.Machine$double.xmin)
  
  para_sub <- cbind(alpha0, alpha2, rho)
  lambda <- c(1, 1, 1)
  
  
  ## Common model storage
  phi_sub_merged_array <- array(rep(NA, N*n_sub_var*(nt+1)), dim = c(N, n_sub_var, nt))
  alpha6_array <- array(rep(NA, N*(nt+1)), dim = c(N, nt+1))
  xJ_array <- array(rep(NA, tt*N*(nt+1)), dim = c(N, tt, nt+1))
  sur_array <- array(rep(NA, tt*N*(nt+1)), dim = c(N, tt, nt+1))
  imm_array <- array(rep(NA, tt*N*(nt+1)), dim = c(N, tt, nt+1))
  w_common_log_array <- array(rep(NA, N*(nt+1)), dim = c(N, nt+1))
  W_common_array <- array(rep(NA, N*(nt+1)), dim = c(N, nt+1))
  A_common <- array(rep(NA, N*nt), dim = c(N, nt))
  
  phi_sub_merged_array[,,1] <- cbind(alpha0, alpha2, rho)
  
  u_pooling_log <- pooling_prior_log_log(lambda = lambda,
                                         alpha = cbind(alpha0, alpha2),
                                         rho = rho,
                                         mu = 0,
                                         sigma = 2,
                                         trun = c(-10, 10),
                                         a = 0,
                                         b = 10)
  
  
  alpha6_0 <- rtruncnorm(N, -10, 10, 0, 2)
  xJ_0 <- matrix(rep(10, tt*N), nrow = N)
  sur_0 <- matrix(rep(10, tt*N), nrow = N)
  imm_0 <- matrix(rep(10, tt*N), nrow = N)
  llike_0 <- dtruncnorm(alpha6_0, -10, 10, 0, 2)
  
  alpha6_array[,1] <- alpha6_0
  xJ_array[,,1] <- xJ_0
  sur_array[,,1] <- sur_0
  imm_array[,,1] <- imm_0
  
  w_common_log_0 <- llike_0
  w_common_log_array[,1] <- w_common_log_0
  W_common_0 <- exp(w_common_log_0 - matrixStats::logSumExp(w_common_log_0))
  W_common_array[,1] <- W_common_0
  
  ESS <- 1 / sum(W_common_0^2)
  if ((ESS < ESS_min) && (!all(W_common_0 == 1))) {
    U <- runif(1, 0, 1)
    A <- Sys_resamp(W=W_common_0, P=N, U=U)
    alpha6_resamp <- alpha6_0[A]
    xJ_resamp <- xJ_0[A,]
    sur_resamp <- sur_0[A,]
    imm_resamp <- imm_0[A,]
    alpha0_resamp <- alpha0[A]
    alpha2_resamp <- alpha2[A]
    rho_resamp <- rho[A]
    u_pooling_log <- u_pooling_log[A]
    w_common_log_0 <- rep(0, N)
    W_common_0 <- rep(1, N)
  } else {
    alpha6_resamp <- alpha6_0
    xJ_resamp <- xJ_0
    sur_resamp <- sur_0
    imm_resamp <- imm_0
    alpha0_resamp <- alpha0
    alpha2_resamp <- alpha2
    rho_resamp <- rho
  }
  
  
  for (i in 1:(nt-1)) {
    if (i == 1) {
      ### Update alpha
      alpha_update <- alpha_incre
      
      ### Update weights
      p_common_like_log <- p_common_likelihood(data = data2,
                                               xJ = xJ_resamp,
                                               sur = sur_resamp,
                                               imm = imm_resamp)
      
      w_common_log <- w_common_log_0 + alpha_incre * (u_pooling_log + p_common_like_log)
      if (max(w_common_log) < min_lim) w_common_log <- w_common_log/1e4
      w_common_log[which(w_common_log < min_lim)] <- min_lim
      w_common_log_array[,2] <- w_common_log
      W_common <- exp(w_common_log - matrixStats::logSumExp(w_common_log))
      W_common_array[,i] <- W_common
      
      ### Resampling -- optionally
      ESS <- 1 / sum(W_common^2)
      if (ESS < ESS_min) {
        U <- runif(1, 0, 1)
        A <- Sys_resamp(W=W_common, P=N, U=U)
        alpha6_resamp <- alpha6_0[A]
        xJ_resamp <- xJ_0[A,]
        sur_resamp <- sur_0[A,]
        imm_resamp <- imm_0[A,]
        alpha0_resamp <- alpha0_resamp[A]
        alpha2_resamp <- alpha2_resamp[A]
        rho_resamp <- rho_resamp[A]
        phi_sub_merged_array[,,1] <- cbind(alpha0_resamp, alpha2_resamp, rho_resamp)
        u_pooling_log <- u_pooling_log[A]
        w_common_log <- w_common_log_array[,i+1] <- rep(0, N)
        W_common <- W_common_array[,i+1] <- rep(1, N)
        A_common[,i] <- A
      } else {
        alpha6_resamp <- alpha6_0
        xJ_resamp <- xJ_0
        sur_resamp <- sur_0
        imm_resamp <- imm_0
        phi_sub_merged_array[,,1] <- cbind(alpha0_resamp, alpha2_resamp, rho_resamp)
        A_common[,i] <- A_common[,i]
      }
      
      
      
      ### MCMC kernel
      psi2_inits <- list(alpha6_inits = alpha6_resamp,
                         xJ_inits = xJ_resamp,
                         sur_inits = sur_resamp,
                         imm_inits = imm_resamp)
      out_common <- dc_melding_smc_jags(data = data2,
                                        alpha0 = alpha0_resamp,
                                        alpha2 = alpha2_resamp,
                                        rho = rho_resamp,
                                        psi2_inits = psi2_inits,
                                        alpha_j = alpha_update,
                                        N = N,
                                        Ntotal = Ntotal,
                                        n_chains = n_chains)
      
      alpha6 <- out_common$alpha6_res
      alpha6_array[,2] <- alpha6
      xJ <- out_common$xJ_res
      xJ_array[,,2] <- xJ
      sur <- out_common$sur_res
      sur_array[,,2] <- sur
      imm <- out_common$imm_res
      imm_array[,,2] <- imm
    }
    
    ### Update alpha
    alpha_update <- alpha_update + alpha_incre
    
    ### Update weights
    p_common_like_log <- p_common_likelihood(data = data2,
                                             xJ = xJ,
                                             sur = sur,
                                             imm = imm)
    
    w_common_log <- w_common_log + alpha_incre * (u_pooling_log + p_common_like_log)
    if (max(w_common_log) < min_lim) w_common_log <- w_common_log/1e4
    w_common_log[which(w_common_log < min_lim)] <- min_lim
    w_common_log_array[,i+2] <- w_common_log
    W_common <- exp(w_common_log - matrixStats::logSumExp(w_common_log))
    W_common_array[,i+2] <- W_common
    
    ### Resampling -- optionally
    ESS <- 1 / sum(W_common^2)
    if (ESS < ESS_min) {
      U <- runif(1, 0, 1)
      A <- Sys_resamp(W=W_common, P=N, U=U)
      alpha6_resamp <- alpha6[A]
      xJ_resamp <- xJ[A,]
      sur_resamp <- sur[A,]
      imm_resamp <- imm[A,]
      alpha0_resamp <- alpha0_resamp[A]
      alpha2_resamp <- alpha2_resamp[A]
      rho_resamp <- rho_resamp[A]
      phi_sub_merged_array[,,i+1] <- cbind(alpha0_resamp, alpha2_resamp, rho_resamp)
      u_pooling_log <- u_pooling_log[A]
      w_common_log <- w_common_log_array[,i+2] <- rep(0, N)
      W_common <- W_common_array[,i+2] <- rep(1, N)
      A_common[,i] <- A
    } else {
      alpha6_resamp <- alpha6
      xJ_resamp <- xJ
      sur_resamp <- sur
      imm_resamp <- imm
      phi_sub_merged_array[,,i+1] <- cbind(alpha0_resamp, alpha2_resamp, rho_resamp)
      A_common[,i] <- A_common[,i]
    }
    
    
    ### MCMC kernel
    psi2_inits <- list(alpha6_inits = alpha6_resamp,
                       xJ_inits = xJ_resamp,
                       sur_inits = sur_resamp,
                       imm_inits = imm_resamp)
    out_common <- dc_melding_smc_jags(data = data2,
                                      alpha0 = alpha0_resamp,
                                      alpha2 = alpha2_resamp,
                                      rho = rho_resamp,
                                      psi2_inits = psi2_inits,
                                      alpha_j = alpha_update,
                                      N = N,
                                      Ntotal = Ntotal,
                                      n_chains = n_chains)
    
    alpha6 <- out_common$alpha6_res
    alpha6_array[,i+2] <- alpha6
    xJ <- out_common$xJ_res
    xJ_array[,,i+2] <- xJ
    sur <- out_common$sur_res
    sur_array[,,i+2] <- sur
    imm <- out_common$imm_res
    imm_array[,,i+2] <- imm
  }
  
  psi2_res <- list(alpha6 = alpha6_array,
                   xJ = xJ_array,
                   sur = sur_array,
                   imm = imm_array)
  
  return(list(psi2=psi2_res,
              W=W_common_array,
              phi=phi_sub_merged_array))
}

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




dc_melding_2.0 <- function(data,
                           N,
                           n_sub,
                           m,
                           alpha,
                           lambda, #powers of sub-priors
                           mu_phi_uni,
                           sigma_phi_uni,
                           mu_phi_mult,
                           Sigma_phi_mult,
                           psi_alpha,
                           psi_beta,
                           pooling="log",
                           Ntotal_sub,
                           Ntotal) {
  ESS_min <- N / 2
  min_lim <- log(.Machine$double.xmin)
  
  data1 <- data$y1
  data3 <- data$y3
  data2 <- data$y2
  
  nt <- 1 / alpha
  
  # Storage
  ## Merging results
  phi_sub_merged_array <- array(rep(NA, N*n_sub*(nt+1)), dim = c(N, n_sub, nt+1))
  
  ## Common model storage
  psi2_array <- array(rep(NA, N*nt), dim = c(N, nt))
  w_common_log <- array(rep(NA, N*nt), dim = c(N, nt))
  W_common <- rep(NA, N)
  W_common_array <- array(rep(NA, N*nt), dim = c(N, nt))
  A_common <- array(rep(NA, N*nt), dim = c(N, nt))
  
  
  invisible(capture.output(out_sub <- dc_melding_sub_jags(data = cbind(data1, data3),
                                                          mu_phi = mu_phi_uni,
                                                          sigma_phi = sigma_phi_uni,
                                                          N = N,
                                                          Ntotal = Ntotal_sub)))
  phi_sub <- out_sub$mu[,Ntotal_sub,]
  w_sub_log <- out_sub$w_log[,Ntotal_sub,]
  W_sub <- out_sub$W[,Ntotal_sub,]
  psi_sub <- out_sub$sigma[,Ntotal_sub,]
  
  
  # Merging
  ## Storage for mN matchings
  phi_sub_mN <- array(rep(NA, m*N*n_sub), dim = c(m*N, n_sub))
  w_sub_log_mN <- array(rep(NA, m*N*n_sub), dim = c(m*N, n_sub))
  W_sub_mN <- array(rep(NA, m*N*n_sub), dim = c(m*N, n_sub))
  
  ## Resampling mN particles
  for (sub_i in 1:n_sub) {
    A <- sample(1:N, m*N, replace = T)
    phi_sub_mN[,sub_i] <- phi_sub[A, sub_i]
    w_sub_log_mN[,sub_i] <- w_sub_log[A, sub_i]
    W_sub_mN[,sub_i] <- W_sub[A, sub_i]
  }
  ### Pooled prior decomposition:
  ### P_pool1(phi12)=1, P_pool3(phi23)=1 --- noninformative priors
  ### P_pool2(phi12, phi23)=P_pool(phi12, phi23)
  u_pooling_log_mN <- pooling_prior_log_log(lambda = lambda,
                                            phi = phi_sub_mN,
                                            mu_uni = matrix(rep(mu_phi_uni, m*N), nrow=m*N, byrow=T),
                                            sigma_uni = matrix(rep(sigma_phi_uni, m*N), nrow=m*N, byrow=T),
                                            mu_mult = mu_phi_mult,
                                            Sigma_mult = Sigma_phi_mult)
  p_approx_log <- normal_uni_log_sum(data = data2,
                                     mu = rowSums(phi_sub_mN),
                                     sigma = rep(psi_alpha/psi_beta, m*N))
  V_log <- alpha * (u_pooling_log_mN + p_approx_log)
  V_log[which(V_log < min_lim)] <- min_lim
  V <- exp(V_log)
  V <- V / sum(V)
  
  A <- sample(1:(m*N), N, replace = T, prob = V)
  phi_sub_merged <- phi_sub_mN[A,]
  phi_sub_merged_array[,,1] <- phi_sub_merged
  W_sub_merged <- W_sub_mN[A,]
  u_pooling_log <- u_pooling_log_mN[A]
  p_approx_log <- p_approx_log[A]
  
  
  # Common model
  psi2_0 <- rgamma(N, psi_alpha, psi_beta)
  w_common_log_0 <- rep(0, N)
  
  ## SMC sampler
  for (i in 1:(nt-1)) {
    if (i == 1) {
      ### Update weights
      p_common_log <- normal_uni_log_sum(data = data2,
                                         mu = rowSums(phi_sub_merged),
                                         sigma = psi2_0)
      p_common_approx_log <- normal_uni_log_sum(data = data2,
                                                mu = rowSums(phi_sub_merged),
                                                sigma = rep(psi_alpha/psi_beta, N))
      w_common_log[,i] <- w_common_log_0 + alpha * (1-alpha) * u_pooling_log + 
        alpha * p_common_log - alpha^2 * p_common_approx_log
      if (max(w_common_log[,i]) < min_lim) w_common_log[,i] <- w_common_log[,i]/1e4
      w_common_log[which(w_common_log[,i] < min_lim), i] <- min_lim
      W_common <- exp(w_common_log[,i] - matrixStats::logSumExp(w_common_log[,i]))
      W_common_array[,i] <- W_common
      
      ### Resampling -- optionally
      ESS <- 1 / sum(W_common^2)
      if (ESS < ESS_min) {
        U <- runif(1, 0, 1)
        A <- Sys_resamp(W=W_common, P=N, U=U)
        psi2_resamp <- psi2_0[A]
        phi_sub_merged <- phi_sub_merged[A,]
        phi_sub_merged_array[,,i+1] <- phi_sub_merged
        u_pooling_log <- u_pooling_log[A]
        w_common_log[,i] <- rep(0, N)
        A_common[,i] <- A
      } else {
        psi2_resamp <- psi2_0
        phi_sub_merged_array[,,i+1] <- phi_sub_merged
      }
      
      ### Update alpha
      alpha_update <- alpha
      
      ### MCMC kernel
      #invisible(capture.output(out_common <- dc_melding_SMC_jags(data = data2,
      #                                                           phi = rowSums(phi_sub_merged),
      #                                                           psi_alpha=psi_alpha,
      #                                                           psi_beta=psi_beta,
      #                                                           psi2_init = psi2_resamp,
      #                                                           alpha = alpha_update,
      #                                                           N = N,
      #                                                           Ntotal = Ntotal)))
      
      out_common <- dc_melding_smc_mcmc(data = data2,
                                        phi = rowSums(phi_sub_merged),
                                        psi_alpha=psi_alpha,
                                        psi_beta=psi_beta,
                                        sigma_init = psi2_resamp,
                                        alpha_j = alpha_update,
                                        N = N,
                                        Ntotal = Ntotal)
      
      psi2 <- out_common$sigma[,Ntotal]
      psi2_array[,i] <- psi2
    }
    
    ### Update weights
    p_common_log <- normal_uni_log_sum(data = data2,
                                       mu = rowSums(phi_sub_merged),
                                       sigma = psi2)
    p_common_approx_log <- normal_uni_log_sum(data = data2,
                                              mu = rowSums(phi_sub_merged),
                                              sigma = rep(psi_alpha/psi_beta, N))
    w_common_log[,i+1] <- w_common_log[,i] + alpha * (1-alpha) * u_pooling_log + 
      alpha * p_common_log - alpha^2 * p_common_approx_log
    if (max(w_common_log[,i+1]) < min_lim) w_common_log[,i+1] <- w_common_log[,i+1]/1e4
    w_common_log[which(w_common_log[,i+1] < min_lim), i+1] <- min_lim
    W_common <- exp(w_common_log[,i+1] - matrixStats::logSumExp(w_common_log[,i+1]))
    W_common_array[,i+1] <- W_common
    
    ### Resampling -- optionally
    ESS <- 1 / sum(W_common^2)
    if (ESS < ESS_min) {
      U <- runif(1, 0, 1)
      A <- Sys_resamp(W=W_common, P=N, U=U)
      psi2_resamp <- psi2[A]
      phi_sub_merged <- phi_sub_merged[A,]
      phi_sub_merged_array[,,i+2] <- phi_sub_merged
      u_pooling_log <- u_pooling_log[A]
      w_common_log[,i+1] <- rep(0, N)
      A_common[,i+1] <- A
    } else {
      psi2_resamp <- psi2
      phi_sub_merged_array[,,i+2] <- phi_sub_merged
    }
    
    ### Update alpha
    alpha_update <- alpha_update + alpha
    
    ### MCMC kernel
    #invisible(capture.output(out_common <- dc_melding_SMC_jags(data = data2,
    #                                                           phi = rowSums(phi_sub_merged),
    #                                                           psi_alpha=psi_alpha,
    #                                                           psi_beta=psi_beta,
    #                                                           psi2_init = psi2_resamp,
    #                                                           alpha = alpha_update,
    #                                                           N = N,
    #                                                           Ntotal = Ntotal)))
    out_common <- dc_melding_smc_mcmc(data = data2,
                                      phi = rowSums(phi_sub_merged),
                                      psi_alpha=psi_alpha,
                                      psi_beta=psi_beta,
                                      sigma_init = psi2_resamp,
                                      alpha_j = alpha_update,
                                      N = N,
                                      Ntotal = Ntotal)
    psi2 <- out_common$sigma[,Ntotal]
    psi2_array[,i+1] <- psi2
  }
  
  return(list(psi2=psi2_array,
              W=W_common_array,
              phi_merged=phi_sub_merged_array,
              V=V))
}








dc_melding_smc_sub1 <- function(data,
                                mu_prior,
                                sigma_prior,
                                alpha_prior,
                                beta_prior,
                                N,
                                Ntotal,
                                n_chains,
                                alpha_incre) {
  time_begin <- Sys.time()
  
  ESS_min <- N / 2
  min_lim <- log(.Machine$double.xmin)
  
  nt <- 1 / alpha_incre
  
  # Storage
  mu_array <- sigma_array <- array(rep(NA, (nt+1)*N), dim = c(N, nt+1))
  w_log_array <- W_array <- array(rep(NA, (nt+1)*N), dim = c(N, nt+1))
  
  # Initialization
  mu_0 <- rnorm(N, mu_prior, sigma_prior)
  mu_array[,1] <- mu_0
  sigma_0 <- rep(sd(data), N)
  sigma_array[,1] <- sigma_0
  
  w_log_0 <- dnorm(mu_0, mu_prior, sigma_prior, T) + 
    dgamma(sigma_0, shape = alpha_prior, rate = beta_prior, log = T)
  w_log_0[which(w_log_0 < min_lim)] <- min_lim
  w_log_array[,1] <- w_log_0
  W_0 <- exp(w_log_0 - matrixStats::logSumExp(w_log_0))
  W_array[,1] <- W_0
  
  for (i in 1:(nt-1)) {
    if (i == 1) {
      ### Resampling -- optionally
      ESS <- 1 / sum(W_0^2)
      if (ESS < ESS_min) {
        U <- runif(1, 0, 1)
        A <- Sys_resamp(W = W_0, P = N, U = U)
        mu_resamp <- mu_0[A]
        sigma_resamp <- sigma_0[A]
        w_log_0 <- w_log_array[,1] <- rep(0, N)
      } else {
        mu_resamp <- mu_0
        sigma_resamp <- sigma_0
      }
      
      ## Update alpha
      alpha_update <- alpha_incre
      
      ## MCMC kernel
      inits <- list(mu_0 = mu_resamp,
                    sigma_0 = sigma_resamp)
      out_jags <- dc_melding_jags_sub1(data = data,
                                       mu_mean = mu_prior,
                                       mu_sd = sigma_prior,
                                       sigma_alpha = alpha_prior,
                                       sigma_beta = beta_prior,
                                       inits = inits,
                                       alpha_j = alpha_update,
                                       N = N,
                                       Ntotal = Ntotal,
                                       n_chains = n_chains)
      mu <- out_jags$mu
      mu_array[,2] <- mu
      sigma <- out_jags$sigma
      sigma_array[,2] <- sigma
      
      llike <- normal_uni_log_particle(data, mu, sigma)
      w_log <- w_log_0 + alpha_incre * llike
      if (max(w_log) < min_lim) w_log <- w_log/1e4
      w_log[which(w_log < min_lim)] <- min_lim
      w_log_array[,2] <- w_log
      W <- exp(w_log - matrixStats::logSumExp(w_log))
      W_array[,2] <- W
    }
    
    ### Resampling -- optionally
    ESS <- 1 / sum(W^2)
    if (ESS < ESS_min) {
      U <- runif(1, 0, 1)
      A <- Sys_resamp(W = W, P = N, U = U)
      mu_resamp <- mu[A]
      sigma_resamp <- sigma[A]
      w_log <- w_log_array[,i+1] <- rep(0, N)
    } else {
      mu_resamp <- mu
      sigma_resamp <- sigma
    }
    
    ## Update alpha
    alpha_update <- alpha_update + alpha_incre
    
    ## MCMC kernel
    inits <- list(mu_0 = mu_resamp,
                  sigma_0 = sigma_resamp)
    out_jags <- dc_melding_jags_sub1(data = data,
                                     mu_mean = mu_prior,
                                     mu_sd = sigma_prior,
                                     sigma_alpha = alpha_prior,
                                     sigma_beta = beta_prior,
                                     inits = inits,
                                     alpha_j = alpha_update,
                                     N = N,
                                     Ntotal = Ntotal,
                                     n_chains = n_chains)
    mu <- out_jags$mu
    mu_array[,i+2] <- mu
    sigma <- out_jags$sigma
    sigma_array[,i+2] <- sigma
    
    llike <- normal_uni_log_particle(data, mu, sigma)
    w_log <- w_log_array[,i+1] + alpha_incre * llike
    if (max(w_log) < min_lim) w_log <- w_log/1e4
    w_log[which(w_log < min_lim)] <- min_lim
    w_log_array[,i+2] <- w_log
    W <- exp(w_log - matrixStats::logSumExp(w_log))
    W_array[,i+2] <- W
  }
  
  time_end <- Sys.time()
  running_time <- difftime(time_end, time_begin)
  
  return(list(mu = mu_array,
              sigma = sigma_array,
              W = W_array,
              running_time = running_time))
}






dc_melding_smc_sub3 <- function(data,
                                mu_mean_prior,
                                mu_sd_prior,
                                tau_alpha_prior,
                                tau_beta_prior,
                                nu_p_prior,
                                N,
                                Ntotal,
                                n_chains,
                                alpha_incre) {
  time_begin <- Sys.time()
  
  ESS_min <- N / 2
  min_lim <- log(.Machine$double.xmin)
  
  nt <- 1 / alpha_incre
  
  # Storage
  mu_array <- tau_array <- nu_array <- array(rep(NA, (nt+1)*N), dim = c(N, nt+1))
  w_log_array <- W_array <- array(rep(NA, (nt+1)*N), dim = c(N, nt+1))
  
  # Initialization
  mu_0 <- rep(0, N) #rnorm(N, mu_prior, sigma_prior)
  mu_array[,1] <- mu_0
  tau_0 <- rep(1/sd(data)^2, N)  #rgamma(N, shape = alpha_prior, rate = beta_prior)
  tau_array[,1] <- tau_0
  nu_0 <- rcat(N, nu_p_prior)
  nu_array[,1] <- nu_0
  
  w_log_0 <- normal_uni_log(mu_0, mu_mean_prior, mu_sd_prior) +
    gamma_log(tau_0, tau_alpha_prior, tau_beta_prior)
  w_log_0 <- rep(0, N)
  if (max(w_log_0) < min_lim) w_log_0 <- w_log_0/1e4
  w_log_0[which(w_log_0 < min_lim)] <- min_lim
  w_log_array[,1] <- w_log_0
  W_0 <- exp(w_log_0 - matrixStats::logSumExp(w_log_0))
  W_array[,1] <- W_0
  
  for (i in 1:(nt-1)) {
    if (i == 1) {
      ### Resampling -- optionally
      ESS <- 1 / sum(W_0^2)
      if (ESS < ESS_min) {
        U <- runif(1, 0, 1)
        A <- Sys_resamp(W = W_0, P = N, U = U)
        mu_resamp <- mu_0[A]
        tau_resamp <- tau_0[A]
        nu_resamp <- nu_0[A]
        w_log_0 <- w_log_array[,1] <- rep(0, N)
      } else {
        mu_resamp <- mu_0
        tau_resamp <- tau_0
        nu_resamp <- nu_0
      }
      
      ## Update alpha
      alpha_update <- alpha_incre
      
      ## MCMC kernel
      inits <- list(mu_0 = mu_resamp,
                    tau_0 = tau_resamp,
                    nu_0 = nu_resamp)
      out_jags <- dc_melding_jags_sub3(data = data,
                                       mu_mean = mu_mean_prior,
                                       mu_sd = mu_sd_prior,
                                       tau_alpha = tau_alpha_prior,
                                       tau_beta = tau_beta_prior,
                                       nu_p = nu_p_prior,
                                       alpha_j = alpha_update,
                                       inits = inits,
                                       N = N,
                                       Ntotal = Ntotal,
                                       n_chains = n_chains)
      mu <- out_jags$mu
      mu_array[,2] <- mu
      tau <- out_jags$tau
      tau_array[,2] <- tau
      nu <- out_jags$nu
      nu_array[,2] <- nu
      
      llike <- tstudent_log_particle(data, mu, tau, nu)
      w_log <- w_log_0 + alpha_incre * llike
      if (max(w_log) < min_lim) w_log <- w_log/1e4
      w_log[which(w_log < min_lim)] <- min_lim
      w_log_array[,2] <- w_log
      W <- exp(w_log - matrixStats::logSumExp(w_log))
      W_array[,2] <- W
    }
    
    ### Resampling -- optionally
    ESS <- 1 / sum(W^2)
    if (ESS < ESS_min) {
      U <- runif(1, 0, 1)
      A <- Sys_resamp(W = W, P = N, U = U)
      mu_resamp <- mu[A]
      tau_resamp <- tau[A]
      nu_resamp <- nu[A]
      w_log <- w_log_array[,i+1] <- rep(0, N)
    } else {
      mu_resamp <- mu
      tau_resamp <- tau
      nu_resamp <- nu
    }
    
    ## Update alpha
    alpha_update <- alpha_update + alpha_incre
    
    ## MCMC kernel
    inits <- list(mu_0 = mu_resamp,
                  tau_0 = tau_resamp,
                  nu_0 = nu_resamp)
    out_jags <- dc_melding_jags_sub3(data = data,
                                     mu_mean = mu_mean_prior,
                                     mu_sd = mu_sd_prior,
                                     tau_alpha = tau_alpha_prior,
                                     tau_beta = tau_beta_prior,
                                     nu_p = nu_p_prior,
                                     alpha_j = alpha_update,
                                     inits = inits,
                                     N = N,
                                     Ntotal = Ntotal,
                                     n_chains = n_chains)
    mu <- out_jags$mu
    mu_array[,i+2] <- mu
    tau <- out_jags$tau
    tau_array[,i+2] <- tau
    nu <- out_jags$nu
    nu_array[,i+2] <- nu
    
    llike <- tstudent_log_particle(data, mu, tau, nu)
    w_log <- w_log + alpha_incre * llike
    if (max(w_log) < min_lim) w_log <- w_log/1e4
    w_log[which(w_log < min_lim)] <- min_lim
    w_log_array[,i+2] <- w_log
    W <- exp(w_log - matrixStats::logSumExp(w_log))
    W_array[,i+2] <- W
  }
  
  time_end <- Sys.time()
  running_time <- difftime(time_end, time_begin)
  
  return(list(mu = mu_array,
              tau = tau_array,
              nu = nu_array,
              W = W_array,
              running_time = running_time))
}







dc_melding_smc_sub5 <- function(data,
                                alpha1_prior,
                                beta1_prior,
                                alpha2_prior,
                                beta2_prior,
                                N,
                                Ntotal,
                                n_chains,
                                alpha_incre) {
  time_begin <- Sys.time()
  
  ESS_min <- N / 2
  min_lim <- log(.Machine$double.xmin)
  
  TT <- length(data)
  nt <- 1 / alpha_incre
  
  # Storage
  sigma1_array <- sigma2_array <- array(rep(NA, (nt+1)*N), dim = c(N, nt+1))
  w_log_array <- W_array <- array(rep(NA, (nt+1)*N), dim = c(N, nt+1))
  
  # Initialization
  sigma1_0 <- rep(5, N)  # rgamma(N, shape = alpha1_prior, rate = beta1_prior)
  sigma1_array[,1] <- sigma1_0
  sigma2_0 <- rgamma(N, shape = alpha2_prior, rate = beta2_prior)
  sigma2_array[,1] <- sigma2_0
  x_0 <- matrix(1, nrow = N, ncol = TT)
  
  w_log_0 <- dgamma(sigma1_0, shape = alpha1_prior, rate = beta1_prior, log = T) + 
    dgamma(sigma2_0, shape = alpha2_prior, rate = beta2_prior, log = T)
  w_log_0 <- w_log_0 + dnorm(x_0[,1], rep(1, N), sigma2_0, log = T)
  for (t in 2:TT) {
    w_log_0 <- w_log_0 + dnorm(x_0[,t], x_0[,t-1], sigma2_0, log = T)
  }
  w_log_0[which(w_log_0 < min_lim)] <- min_lim
  w_log_array[,1] <- w_log_0
  W_0 <- exp(w_log_0 - matrixStats::logSumExp(w_log_0))
  W_array[,1] <- W_0
  
  for (i in 1:(nt-1)) {
    if (i == 1) {
      ### Resample -- optionally
      ESS <- 1 / sum(W_0^2)
      if (ESS < ESS_min) {
        U <- runif(1, 0, 1)
        A <- Sys_resamp(W = W_0, P = N, U = U)
        sigma1_resamp <- sigma1_0[A]
        sigma2_resamp <- sigma2_0[A]
        x_resamp <- x_0[A,]
        w_log_0 <- w_log_array[,1] <- rep(0, N)
      } else {
        sigma1_resamp <- sigma1_0
        sigma2_resamp <- sigma2_0
        x_resamp <- x_0
      }
      
      ## Update alpha
      alpha_update <- alpha_incre
      
      ## MCMC kernel
      inits <- list(sigma1_inits = sigma1_resamp,
                    sigma2_inits = sigma2_resamp,
                    x_inits = x_resamp)
      out_jags <- dc_melding_jags_sub5(data = data,
                                       sigma1_alpha = alpha1_prior,
                                       sigma1_beta = beta1_prior,
                                       sigma2_alpha = alpha2_prior,
                                       sigma2_beta = beta2_prior,
                                       inits = inits,
                                       alpha_j = alpha_update,
                                       N = N,
                                       Ntotal = Ntotal,
                                       n_chains = n_chains)
      sigma1 <- out_jags$sigma1
      sigma1_array[,2] <- sigma1
      sigma2 <- out_jags$sigma2
      sigma2_array[,2] <- sigma2
      x <- out_jags$x
      
      llike <- dnorm(data[1], x[,1], sigma1, log = T)
      for (t in 2:TT) {
        llike <- llike + dnorm(data[t], x[,t], sigma1, log = T)
      }
      w_log <- w_log_0 + alpha_incre * llike
      if (max(w_log) < min_lim) w_log <- w_log/1e4
      w_log[which(w_log < min_lim)] <- min_lim
      w_log_array[,2] <- w_log
      W <- exp(w_log - matrixStats::logSumExp(w_log))
      W_array[,2] <- W
    }
    
    ### Resample -- optionally
    ESS <- 1 / sum(W^2)
    if (ESS < ESS_min) {
      U <- runif(1, 0, 1)
      A <- Sys_resamp(W = W_0, P = N, U = U)
      sigma1_resamp <- sigma1[A]
      sigma2_resamp <- sigma2[A]
      x_resamp <- x[A,]
      w_log <- w_log_array[,i+1] <- rep(0, N)
    } else {
      sigma1_resamp <- sigma1
      sigma2_resamp <- sigma2
      x_resamp <- x
    }
    
    ## Update alpha
    alpha_update <- alpha_update + alpha_incre
    
    ## MCMC kernel
    inits <- list(sigma1_inits = sigma1_resamp,
                  sigma2_inits = sigma2_resamp,
                  x_inits = x_resamp)
    out_jags <- dc_melding_jags_sub5(data = data,
                                     sigma1_alpha = alpha1_prior,
                                     sigma1_beta = beta1_prior,
                                     sigma2_alpha = alpha2_prior,
                                     sigma2_beta = beta2_prior,
                                     inits = inits,
                                     alpha_j = alpha_update,
                                     N = N,
                                     Ntotal = Ntotal,
                                     n_chains = n_chains)
    sigma1 <- out_jags$sigma1
    sigma1_array[,i+2] <- sigma1
    sigma2 <- out_jags$sigma2
    sigma2_array[,i+2] <- sigma2
    x <- out_jags$x
    
    llike <- dnorm(data[1], x[,1], sigma1, log = T)
    for (t in 2:TT) {
      llike <- llike + dnorm(data[t], x[,t], sigma1, log = T)
    }
    w_log <- w_log + alpha_incre * llike
    if (max(w_log) < min_lim) w_log <- w_log/1e4
    w_log[which(w_log < min_lim)] <- min_lim
    w_log_array[,i+2] <- w_log
    W <- exp(w_log - matrixStats::logSumExp(w_log))
    W_array[,i+2] <- W
  }
  
  time_end <- Sys.time()
  running_time <- difftime(time_end, time_begin)
  
  return(list(sigma1 = sigma1_array,
              sigma2 = sigma2_array,
              W = W_array,
              running_time = running_time))
}







dc_melding_smc_sub11 <- function(data,
                                 mu_prior,
                                 sigma_prior,
                                 alpha_prior,
                                 beta_prior,
                                 N,
                                 Ntotal,
                                 n_chains,
                                 alpha_incre) {
  time_begin <- Sys.time()
  
  ESS_min <- N / 2
  min_lim <- log(.Machine$double.xmin)
  
  nt <- 1 / alpha_incre
  
  # Storage
  mu_array <- sigma_array <- array(rep(NA, (nt+1)*N), dim = c(N, nt+1))
  w_log_array <- W_array <- array(rep(NA, (nt+1)*N), dim = c(N, nt+1))
  
  # Initialization
  mu_0 <- rnorm(N, mu_prior, sigma_prior)
  mu_array[,1] <- mu_0
  sigma_0 <- rep(sd(data), N)  # rgamma(N, shape = alpha_prior, rate = beta_prior)
  sigma_array[,1] <- sigma_0
  
  w_log_0 <- dnorm(mu_0, mu_prior, sigma_prior, T) + 
    dgamma(sigma_0, shape = alpha_prior, rate = beta_prior, log = T)
  w_log_0[which(w_log_0 < min_lim)] <- min_lim
  w_log_array[,1] <- w_log_0
  W_0 <- exp(w_log_0 - matrixStats::logSumExp(w_log_0))
  W_array[,1] <- W_0
  
  for (i in 1:(nt-1)) {
    if (i == 1) {
      ### Resampling -- optionally
      ESS <- 1 / sum(W_0^2)
      if (ESS < ESS_min) {
        U <- runif(1, 0, 1)
        A <- Sys_resamp(W = W_0, P = N, U = U)
        mu_resamp <- mu_0[A]
        sigma_resamp <- sigma_0[A]
        w_log_0 <- w_log_array[,1] <- rep(0, N)
      } else {
        mu_resamp <- mu_0
        sigma_resamp <- sigma_0
      }
      
      ## Update alpha
      alpha_update <- alpha_incre
      
      ## MCMC kernel
      inits <- list(mu_0 = mu_resamp,
                    sigma_0 = sigma_resamp)
      out_jags <- dc_melding_jags_sub11(data = data,
                                        mu_mean = mu_prior,
                                        mu_sd = sigma_prior,
                                        sigma_alpha = alpha_prior,
                                        sigma_beta = beta_prior,
                                        inits = inits,
                                        alpha_j = alpha_update,
                                        N = N,
                                        Ntotal = Ntotal,
                                        n_chains = n_chains)
      mu <- out_jags$mu
      mu_array[,2] <- mu
      sigma <- out_jags$sigma
      sigma_array[,2] <- sigma
      
      llike <- normal_uni_log_particle(data, mu, sigma)
      w_log <- w_log_0 + alpha_incre * llike
      if (max(w_log) < min_lim) w_log <- w_log/1e4
      w_log[which(w_log < min_lim)] <- min_lim
      w_log_array[,2] <- w_log
      W <- exp(w_log - matrixStats::logSumExp(w_log))
      W_array[,2] <- W
    }
    
    ### Resampling -- optionally
    ESS <- 1 / sum(W^2)
    if (ESS < ESS_min) {
      U <- runif(1, 0, 1)
      A <- Sys_resamp(W = W, P = N, U = U)
      mu_resamp <- mu[A]
      sigma_resamp <- sigma[A]
      w_log <- w_log_array[,i+1] <- rep(0, N)
    } else {
      mu_resamp <- mu
      sigma_resamp <- sigma
    }
    
    ## Update alpha
    alpha_update <- alpha_update + alpha_incre
    
    ## MCMC kernel
    inits <- list(mu_0 = mu_resamp,
                  sigma_0 = sigma_resamp)
    out_jags <- dc_melding_jags_sub11(data = data,
                                      mu_mean = mu_prior,
                                      mu_sd = sigma_prior,
                                      sigma_alpha = alpha_prior,
                                      sigma_beta = beta_prior,
                                      inits = inits,
                                      alpha_j = alpha_update,
                                      N = N,
                                      Ntotal = Ntotal,
                                      n_chains = n_chains)
    mu <- out_jags$mu
    mu_array[,i+2] <- mu
    sigma <- out_jags$sigma
    sigma_array[,i+2] <- sigma
    
    llike <- normal_uni_log_particle(data, mu, sigma)
    w_log <- w_log_array[,i+1] + alpha_incre * llike
    if (max(w_log) < min_lim) w_log <- w_log/1e4
    w_log[which(w_log < min_lim)] <- min_lim
    w_log_array[,i+2] <- w_log
    W <- exp(w_log - matrixStats::logSumExp(w_log))
    W_array[,i+2] <- W
  }
  
  time_end <- Sys.time()
  running_time <- difftime(time_end, time_begin)
  
  return(list(mu = mu_array,
              sigma = sigma_array,
              W = W_array,
              running_time = running_time))
}







dc_melding_smc_sub9 <- function(data,
                                mu_mean_prior,
                                mu_sd_prior,
                                tau_alpha_prior,
                                tau_beta_prior,
                                nu_p_prior,
                                N,
                                Ntotal,
                                n_chains,
                                alpha_incre) {
  time_begin <- Sys.time()
  
  ESS_min <- N / 2
  min_lim <- log(.Machine$double.xmin)
  
  nt <- 1 / alpha_incre
  
  # Storage
  mu_array <- tau_array <- nu_array <- array(rep(NA, (nt+1)*N), dim = c(N, nt+1))
  w_log_array <- W_array <- array(rep(NA, (nt+1)*N), dim = c(N, nt+1))
  
  # Initialization
  mu_0 <- rep(0, N) #rnorm(N, mu_prior, sigma_prior)
  mu_array[,1] <- mu_0
  tau_0 <- rep(1/sd(data)^2, N)  #rgamma(N, shape = alpha_prior, rate = beta_prior)
  tau_array[,1] <- tau_0
  nu_0 <- rcat(N, nu_p_prior)
  nu_array[,1] <- nu_0
  
  w_log_0 <- normal_uni_log(mu_0, mu_mean_prior, mu_sd_prior) + 
    gamma_log(tau_0, tau_alpha_prior, tau_beta_prior)
  w_log_0 <- rep(0, N)
  if (max(w_log_0) < min_lim) w_log_0 <- w_log_0/1e4
  w_log_0[which(w_log_0 < min_lim)] <- min_lim
  w_log_array[,1] <- w_log_0
  W_0 <- exp(w_log_0 - matrixStats::logSumExp(w_log_0))
  W_array[,1] <- W_0
  
  for (i in 1:(nt-1)) {
    if (i == 1) {
      ### Resampling -- optionally
      ESS <- 1 / sum(W_0^2)
      if (ESS < ESS_min) {
        U <- runif(1, 0, 1)
        A <- Sys_resamp(W = W_0, P = N, U = U)
        mu_resamp <- mu_0[A]
        tau_resamp <- tau_0[A]
        nu_resamp <- nu_0[A]
        w_log_0 <- w_log_array[,1] <- rep(0, N)
      } else {
        mu_resamp <- mu_0
        tau_resamp <- tau_0
        nu_resamp <- nu_0
      }
      
      ## Update alpha
      alpha_update <- alpha_incre
      
      ## MCMC kernel
      inits <- list(mu_0 = mu_resamp,
                    tau_0 = tau_resamp,
                    nu_0 = nu_resamp)
      out_jags <- dc_melding_jags_sub9(data = data,
                                       mu_mean = mu_mean_prior,
                                       mu_sd = mu_sd_prior,
                                       tau_alpha = tau_alpha_prior,
                                       tau_beta = tau_beta_prior,
                                       nu_p = nu_p_prior,
                                       alpha_j = alpha_update,
                                       inits = inits,
                                       N = N,
                                       Ntotal = Ntotal,
                                       n_chains = n_chains)
      mu <- out_jags$mu
      mu_array[,2] <- mu
      tau <- out_jags$tau
      tau_array[,2] <- tau
      nu <- out_jags$nu
      nu_array[,2] <- nu
      
      llike <- tstudent_log_particle(data, mu, tau, nu)
      w_log <- w_log_0 + alpha_incre * llike
      if (max(w_log) < min_lim) w_log <- w_log/1e4
      w_log[which(w_log < min_lim)] <- min_lim
      w_log_array[,2] <- w_log
      W <- exp(w_log - matrixStats::logSumExp(w_log))
      W_array[,2] <- W
    }
    
    ### Resampling -- optionally
    ESS <- 1 / sum(W^2)
    if (ESS < ESS_min) {
      U <- runif(1, 0, 1)
      A <- Sys_resamp(W = W, P = N, U = U)
      mu_resamp <- mu[A]
      tau_resamp <- tau[A]
      nu_resamp <- nu[A]
      w_log <- w_log_array[,i+1] <- rep(0, N)
    } else {
      mu_resamp <- mu
      tau_resamp <- tau
      nu_resamp <- nu
    }
    
    ## Update alpha
    alpha_update <- alpha_update + alpha_incre
    
    ## MCMC kernel
    inits <- list(mu_0 = mu_resamp,
                  tau_0 = tau_resamp,
                  nu_0 = nu_resamp)
    out_jags <- dc_melding_jags_sub9(data = data,
                                     mu_mean = mu_mean_prior,
                                     mu_sd = mu_sd_prior,
                                     tau_alpha = tau_alpha_prior,
                                     tau_beta = tau_beta_prior,
                                     nu_p = nu_p_prior,
                                     alpha_j = alpha_update,
                                     inits = inits,
                                     N = N,
                                     Ntotal = Ntotal,
                                     n_chains = n_chains)
    mu <- out_jags$mu
    mu_array[,i+2] <- mu
    tau <- out_jags$tau
    tau_array[,i+2] <- tau
    nu <- out_jags$nu
    nu_array[,i+2] <- nu
    
    llike <- tstudent_log_particle(data, mu, tau, nu)
    w_log <- w_log + alpha_incre * llike
    if (max(w_log) < min_lim) w_log <- w_log/1e4
    w_log[which(w_log < min_lim)] <- min_lim
    w_log_array[,i+2] <- w_log
    W <- exp(w_log - matrixStats::logSumExp(w_log))
    W_array[,i+2] <- W
  }
  
  time_end <- Sys.time()
  running_time <- difftime(time_end, time_begin)
  
  return(list(mu = mu_array,
              tau = tau_array,
              nu = nu_array,
              W = W_array,
              running_time = running_time))
}








dc_melding_smc_sub7 <- function(data,
                                mu1_prior,
                                sigma1_prior,
                                mu2_prior,
                                sigma2_prior,
                                N,
                                Ntotal,
                                n_chains,
                                alpha_incre) {
  time_begin <- Sys.time()
  
  ESS_min <- N / 2
  min_lim <- log(.Machine$double.xmin)
  
  TT <- length(data)
  nt <- 1 / alpha_incre
  
  # Storage
  mu1_array <- mu2_array <- array(rep(NA, (nt+1)*N), dim = c(N, nt+1))
  w_log_array <- W_array <- array(rep(NA, (nt+1)*N), dim = c(N, nt+1))
  
  # Initialization
  mu1_0 <- rnorm(N, mu1_prior, sigma1_prior)
  mu1_array[,1] <- mu1_0
  mu2_0 <- rnorm(N, mu2_prior, sigma2_prior)
  mu2_array[,1] <- mu2_0
  x_0 <- matrix(rep(exp(mu1_0), TT), nrow = N)
  
  w_log_0 <- dnorm(mu1_0, mu1_prior, sigma1_prior, log = T) +
    dnorm(mu2_0, mu2_prior, sigma2_prior, log = T)
  w_log_0 <- w_log_0 + dlnorm(x_0[,1], mu1_0, rep(.1, N), log = T)
  for (t in 2:TT) {
    w_log_0 <- w_log_0 + dlnorm(x_0[,t], x_0[,t-1] + mu1_0, rep(.1, N), log = T)
  }
  w_log_0[which(w_log_0 < min_lim)] <- min_lim
  w_log_array[,1] <- w_log_0
  W_0 <- exp(w_log_0 - matrixStats::logSumExp(w_log_0))
  W_array[,1] <- W_0
  
  for (i in 1:(nt-1)) {
    if (i == 1) {
      ### Resample -- optionally
      ESS <- 1 / sum(W_0^2)
      if (ESS < ESS_min) {
        U <- runif(1, 0, 1)
        A <- Sys_resamp(W = W_0, P = N, U = U)
        mu1_resamp <- mu1_0[A]
        mu2_resamp <- mu2_0[A]
        x_resamp <- x_0[A,]
        w_log_0 <- w_log_array[,1] <- rep(0, N)
      } else {
        mu1_resamp <- mu1_0
        mu2_resamp <- mu2_0
        x_resamp <- x_0
      }
      
      ## Update alpha
      alpha_update <- alpha_incre
      
      ## MCMC kernel
      inits <- list(mu1_inits = mu1_resamp,
                    mu2_inits = mu2_resamp,
                    x_inits = x_resamp)
      out_jags <- dc_melding_jags_sub7(data = data,
                                       mu1_mu = mu1_prior,
                                       mu1_sigma = sigma1_prior,
                                       mu2_mu = mu2_prior,
                                       mu2_sigma = sigma2_prior,
                                       inits = inits,
                                       alpha_j = alpha_update,
                                       N = N,
                                       Ntotal = Ntotal,
                                       n_chains = n_chains)
      mu1 <- out_jags$mu1
      mu1_array[,2] <- mu1
      mu2 <- out_jags$mu2
      mu2_array[,2] <- mu2
      x <- out_jags$x
      
      llike <- dnorm(data[1], mu2, x[,1], log = T)
      for (t in 2:TT) {
        llike <- llike + dnorm(data[t], mu2, x[,t], log = T)
      }
      w_log <- w_log_0 + alpha_incre * llike
      if (max(w_log) < min_lim) w_log <- w_log/1e4
      w_log[which(w_log < min_lim)] <- min_lim
      w_log_array[,2] <- w_log
      W <- exp(w_log - matrixStats::logSumExp(w_log))
      W_array[,2] <- W
    }
    
    ### Resample -- optionally
    ESS <- 1 / sum(W^2)
    if (ESS < ESS_min) {
      U <- runif(1, 0, 1)
      A <- Sys_resamp(W = W_0, P = N, U = U)
      mu1_resamp <- mu1[A]
      mu2_resamp <- mu2[A]
      x_resamp <- x[A,]
      w_log <- w_log_array[,i+1] <- rep(0, N)
    } else {
      mu1_resamp <- mu1
      mu2_resamp <- mu2
      x_resamp <- x
    }
    
    ## Update alpha
    alpha_update <- alpha_update + alpha_incre
    
    ## MCMC kernel
    inits <- list(mu1_inits = mu1_resamp,
                  mu2_inits = mu2_resamp,
                  x_inits = x_resamp)
    out_jags <- dc_melding_jags_sub7(data = data,
                                     mu1_mu = mu1_prior,
                                     mu1_sigma = sigma1_prior,
                                     mu2_mu = mu2_prior,
                                     mu2_sigma = sigma2_prior,
                                     inits = inits,
                                     alpha_j = alpha_update,
                                     N = N,
                                     Ntotal = Ntotal,
                                     n_chains = n_chains)
    mu1 <- out_jags$mu1
    mu1_array[,i+2] <- mu1
    mu2 <- out_jags$mu2
    mu2_array[,i+2] <- mu2
    x <- out_jags$x
    
    llike <- dnorm(data[1], mu2, x[,1], log = T)
    for (t in 2:TT) {
      llike <- llike + dnorm(data[t], mu2, x[,t], log = T)
    }
    w_log <- w_log + alpha_incre * llike
    if (max(w_log) < min_lim) w_log <- w_log/1e4
    w_log[which(w_log < min_lim)] <- min_lim
    w_log_array[,i+2] <- w_log
    W <- exp(w_log - matrixStats::logSumExp(w_log))
    W_array[,i+2] <- W
  }
  
  time_end <- Sys.time()
  running_time <- difftime(time_end, time_begin)
  
  return(list(mu1 = mu1_array,
              mu2 = mu2_array,
              W = W_array,
              running_time = running_time))
}
















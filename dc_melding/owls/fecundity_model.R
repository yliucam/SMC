fecundity_model <- function(data,
                            lower_prior,
                            upper_prior,
                            N,
                            Ntotal,
                            n_chains,
                            alpha) {
  time_begin <- Sys.time()
  
  ESS_min <- N / 2
  min_lim <- log(.Machine$double.xmin)
  
  ti <- dim(data)[1]
  nt <- 1 / alpha
  
  # Storage
  rho_array <- array(rep(NA, (nt+1)*N), dim = c(N, nt+1))
  w_log_array <- array(rep(NA, (nt+1)*N), dim = c(N, nt+1))
  W_array <- array(rep(NA, (nt+1)*N), dim = c(N, nt+1))
  
  # Initialization
  # rho_0 <- rep(5, N)
  rho_0 <- runif(N, 0, 10)
  rho_array[,1] <- rho_0
  
  w_log_0 <- rep(log(1/10), N)
  w_log_array[,1] <- w_log_0
  W_0 <- exp(w_log_0 - matrixStats::logSumExp(w_log_0))
  W_array[,1] <- W_0
  
  for (i in 1:(nt-1)){
    if (i == 1) {
      ## Resampling -- optionally (no need actually for this initialization)
      ESS <- 1 / sum(W_0^2)
      if (ESS < ESS_min) {
        U <- runif(1, 0, 1)
        A <- Sys_resamp(W = W_0, P = N, U = U)
        rho_resamp <- rho_0[A]
        w_log_0 <- w_log_array[,1] <- rep(0, N)
        W_0 <- W_array[,1] <- rep(1, N)
      } else {
        rho_resamp <- rho_0
      }
      
      ## Update alpha
      alpha_update <- alpha
      
      ## MCMC kernel
      inits <- list(rho_0 = rho_resamp)
      out_jags <- fecundity_jags(data = data,
                                 lower_prior = lower_prior,
                                 upper_prior = upper_prior,
                                 alpha_j = alpha_update,
                                 inits = inits,
                                 N = N,
                                 Ntotal = Ntotal,
                                 n_chains = n_chains)
      rho <- out_jags$rho_jags
      rho_array[,2] <- rho
      
      ## Update weights
      llike <- fecundity_llike(data = data, rho = rho)
      w_log <- w_log_0 + alpha * llike
      if (max(w_log) < min_lim) w_log <- w_log/1e4
      w_log[which(w_log < min_lim)] <- min_lim
      w_log_array[,2] <- w_log
      W <- exp(w_log - matrixStats::logSumExp(w_log))
      W_array[,2] <- W
    }
    
    ## Resampling -- optionally (no need actually for this initialization)
    ESS <- 1 / sum(W^2)
    if (ESS < ESS_min) {
      U <- runif(1, 0, 1)
      A <- Sys_resamp(W = W, P = N, U = U)
      rho_resamp <- rho[A]
      w_log <- w_log_array[,i+1] <- rep(0, N)
      W <- W_array[,i+1] <- rep(1, N)
    } else {
      rho_resamp <- rho
    }
    
    ## Update alpha
    alpha_update <- alpha_update + alpha
    
    ## MCMC kernel
    inits <- list(rho_0 = rho_resamp)
    out_jags <- fecundity_jags(data = data,
                               lower_prior = lower_prior,
                               upper_prior = upper_prior,
                               alpha_j = alpha_update,
                               inits = inits,
                               N = N,
                               Ntotal = Ntotal,
                               n_chains = n_chains)
    rho <- out_jags$rho_jags
    rho_array[,i+2] <- rho
    
    ## Update weights
    llike <- fecundity_llike(data = data, rho = rho)
    w_log <- w_log + alpha * llike
    if (max(w_log) < min_lim) w_log <- w_log/1e4
    w_log[which(w_log < min_lim)] <- min_lim
    w_log_array[,i+2] <- w_log
    W <- exp(w_log - matrixStats::logSumExp(w_log))
    W_array[,i+2] <- W
  }
  
  time_end <- Sys.time()
  time_run <- time_end - time_begin
  
  return(list(rho = rho_array,
              W = W_array,
              tim = time_run))
}




fecundity_jags <- function(data,
                           lower_prior,
                           upper_prior,
                           alpha_j,
                           inits,
                           N,
                           Ntotal,
                           n_chains) {
  ti <- dim(data)[1]
  
  DATA_Nt <- matrix(rep(data[,1], N/n_chains), nrow = N/n_chains, byrow = T)
  DATA_nt <- matrix(rep(data[,2], N/n_chains), nrow = N/n_chains, byrow = T)
  
  rho_0 <- matrix(inits$rho_0, ncol = n_chains)
  
  model <- "model{
    for (i in 1:N) {
      rho[i] ~ dunif(lower, upper)
      
      for (j in 1:ti) {
        rate[i,j] <- Nt[i,j] * rho[i]
        logL[i,j] <- logdensity.pois(nt[i,j], rate[i,j])
        phi[i,j] <- -alpha_j * logL[i,j] + 1000
        zeros[i,j] ~ dpois(phi[i,j])
      }
    }
  }"
  
  ## Parallel computation
  cl <- makeCluster(n_chains)
  registerDoParallel(cl)
  clusterExport(cl, c("DATA_Nt", "DATA_nt", "rho_0", "model"), envir = environment())
  
  mcmc_list <- foreach(i = 1:n_chains, .packages = "rjags") %dopar% {
    foo <- jags.model(textConnection(model),
                      data = list(Nt = DATA_Nt,
                                  nt = DATA_nt,
                                  lower = lower_prior,
                                  upper = upper_prior,
                                  alpha_j = alpha_j,
                                  zeros = matrix(rep(0, ti*N/n_chains), ncol = ti),
                                  N = N/n_chains,
                                  ti = ti),
                      inits = list(rho = rho_0[,i]),
                      quiet = T)
    jags <- coda.samples(foo,
                         variable.names = c("rho"),
                         n.iter = Ntotal,
                         progress.bar = "none")
  }
  
  stopCluster(cl)
  
  jags_out <- as.matrix(mcmc_list)
  jags_out1 <- as.matrix(mcmc_list[[1]]) # This copy is for extracting the column names
  
  # In terms of the partial matching of variable/column names to extract and store posterior samples
  rho_res <- lapply(mcmc_list, function(jags_out) jags_out[,grep("rho", colnames(jags_out1))])
  rho_res <- unlist(rho_res)
  rho <- matrix(rho_res, nrow = Ntotal)[Ntotal,]
  
  return(list(rho_jags = rho))
}

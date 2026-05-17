dc_melding_jags_sub6 <- function(data,
                                 mu,
                                 sigma,
                                 alpha_prior,
                                 beta_prior,
                                 inits,
                                 alpha_j,
                                 N,
                                 Ntotal,
                                 n_chains) {
  TT <- length(data)
  
  DATA <- matrix(rep(data, N/n_chains), nrow = N/n_chains, byrow = T)
  
  rho_0 <- matrix(inits$rho_inits, ncol = n_chains)
  x_0 <- array(inits$x_inits, dim = c(N/n_chains, n_chains, TT))
  x_0 <- aperm(x_0, c(1,3,2))
  
  inits <- list(rho_inits = rho_0,
                x_inits = x_0)
  
  model <- "model{
    for (i in 1:N) {
      rho[i] ~ dbeta(alpha, beta)
      
      tau[i] <- pow(sigma[i], -2)
      x_0[i] ~ dnorm(mu[i], tau[i])
      xmean[i,1] <- mu[i] + rho[i] * (x_0[i] - mu[i])
      x[i,1] ~ dnorm(xmean[i,1], tau[i])
      Ytau[i,1] <- 1/exp(x[i,1])
      Y[i,1] ~ dnorm(0, alpha_j * Ytau[i,1])
      for (t in 2:TT) {
        xmean[i,t] <- mu[i] + rho[i] * (x[i,t-1] - mu[i])
        x[i,t] ~ dnorm(xmean[i,t], tau[i])
        Ytau[i,t] <- 1/exp(x[i,t])
        Y[i,t] ~ dnorm(0, alpha_j * Ytau[i,t])
      }
    }
  }"
  
  
  cl <- makeCluster(n_chains)
  registerDoParallel(cl)
  clusterExport(cl, c("DATA", "inits", "mu", "sigma", "model"),
                envir = environment())
  
  mcmc_list <- foreach(i = 1:n_chains, .packages = "rjags") %dopar% {
    foo <- jags.model(textConnection(model),
                      data = list(Y = DATA,
                                  mu = mu,
                                  sigma = sigma,
                                  alpha = alpha_prior,
                                  beta = beta_prior,
                                  alpha_j = alpha_j,
                                  N = N/n_chains,
                                  TT = TT),
                      inits = list(rho = inits$rho_inits[,i],
                                   x = inits$x_inits[,,i]))
    
    out_jags <- coda.samples(foo,
                             variable.names = c("rho", "x"),
                             n.iter = Ntotal)
  }
  
  stopCluster(cl)
  
  jags_out <- as.matrix(mcmc_list)
  jags_out1 <- as.matrix(mcmc_list[[1]])
  
  var_names <- c("rho", "x")
  var_names_index <- 1:TT
  
  var_names_full <- rep(NA, TT+1)
  var_names_full[1] <- paste0("^", var_names[1], "\\[")
  for (t in 1:TT) {
    var_names_full[t+1] <- paste0("^", var_names[2], "\\[[0-9]+,", var_names_index[t], "\\]")
  }
  
  rho_res <- lapply(mcmc_list, function(jags_out) jags_out[,grep(var_names_full[1], colnames(jags_out1))])
  rho_res <- unlist(rho_res)
  rho_res <- matrix(rho_res, nrow = Ntotal)[Ntotal,]
  
  x_res <- matrix(rep(NA, TT*N), ncol = TT)
  for (t in 1:TT) {
    x_temp <- lapply(mcmc_list, function(jags_out) jags_out[,grep(var_names_full[t+1], colnames(jags_out1))])
    x_temp <- unlist(x_temp)
    x_res[,t] <- matrix(x_temp, nrow = Ntotal)[Ntotal,]
  }
  
  return(list(rho = rho_res,
              x = x_res))
}







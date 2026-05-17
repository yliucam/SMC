# JAGs for t-student distribution t(mu, tau, nu) - mu: location
#                                                - tau: precision > 0
#                                                - nu: degrees of freedom -- natural numbers
dc_melding_jags_sub2 <- function(data,
                                 mu,
                                 tau,
                                 nu_p,
                                 alpha_j,
                                 inits,
                                 N,
                                 Ntotal,
                                 n_chains) { 
  n <- length(data)
  
  DATA <- matrix(rep(data, N/n_chains), nrow = N/n_chains, byrow = T)
  
  model <- "model {
    for (i in 1:N) {
      nu[i] ~ dcat(nu_p[])
      
      for (j in 1:n) {
        logL[i,j] <- loggam((nu[i]+1)/2) - loggam(nu[i]/2) + .5 * log(tau[i]/nu[i]) -
          .5 * (nu[i]+1) * log(1 + tau[i] * pow((Y[i,j]-mu[i]), 2) / nu[i])
        phi[i,j] <- -alpha_j * logL[i,j] + 1000
        zeros[i,j] ~ dpois(phi[i,j])
      }
    }
  }"
  
  nu_0 <- matrix(inits$nu_inits, ncol = n_chains)
  
  # Parallel computation
  cl <- makeCluster(n_chains)
  registerDoParallel(cl)
  clusterExport(cl, c("DATA", "nu_0", "mu", "tau", "nu_p", "model"), envir = environment())
  
  mcmc_list <- foreach(i = 1:n_chains, .packages = "rjags") %dopar% {
    foo <- jags.model(textConnection(model),
                      data = list(Y = DATA,
                                  mu = mu,
                                  tau = tau,
                                  nu_p = nu_p,
                                  alpha_j = alpha_j,
                                  zeros = matrix(rep(0, n*N/n_chains), ncol=n),
                                  n = n,
                                  N = N/n_chains),
                      inits = list(nu = nu_0[,i]),
                      n.chains = 1,
                      quiet = T)
    
    # update(foo, burn_in)
    jags <- coda.samples(model = foo, 
                         variable.names = c("nu"),
                         n.iter = Ntotal,
                         progress.bar = "none")
  }
  
  stopCluster(cl)
  
  jags_out <- as.matrix(mcmc_list)
  jags_out1 <- as.matrix(mcmc_list[[1]]) # This copy is for extracting the column names
  
  # In terms of the partial matching of variable/column names to extract and store posterior samples
  nu_res <- lapply(mcmc_list, function(jags_out) jags_out[,grep("nu", colnames(jags_out1))])
  nu_res <- unlist(nu_res)
  nu_res <- matrix(nu_res, nrow = Ntotal)[Ntotal,]
  
  return(list(nu = nu_res))
}
dc_melding_jags_sub1 <- function(data,
                                 mu_mean,
                                 mu_sd,
                                 sigma_alpha,
                                 sigma_beta,
                                 inits,
                                 alpha_j,
                                 N,
                                 Ntotal,
                                 n_chains) {
  n <- length(data)
  
  DATA <- matrix(rep(data, N/n_chains), nrow = N/n_chains, byrow = T)
  
  mu_0 <- matrix(inits$mu_0, ncol = n_chains)
  sigma_0 <- matrix(inits$sigma_0, ncol = n_chains)
  
  model <- "model{
    for (i in 1:N) {
      mu[i] ~ dnorm(mu_mean, 1/mu_var)
      sigma[i] ~ dgamma(sigma_alpha, sigma_beta) #dunif(sigma_lower, sigma_upper)
      
      tau[i] <- pow(sigma[i], -2)
      
      for (j in 1:n) {
        Y[i,j] ~ dnorm(mu[i], alpha_j*tau[i])
      }
    }
  }"
  
  
  cl <- makeCluster(n_chains)
  registerDoParallel(cl)
  clusterExport(cl, c("DATA", "mu_0", "sigma_0", "sigma_alpha", "sigma_beta",
                      "model"), 
                envir = environment())
  
  mcmc_list <- foreach(i = 1:n_chains, .packages = "rjags") %dopar% {
    foo <- jags.model(textConnection(model),
                      data = list(Y = DATA,
                                  mu_mean = mu_mean,
                                  mu_var = mu_sd^2,
                                  sigma_alpha = sigma_alpha,
                                  sigma_beta = sigma_beta,
                                  alpha_j = alpha_j,
                                  n = n,
                                  N = N/n_chains),
                      inits = list(mu = mu_0[,i],
                                   sigma = sigma_0[,i]),
                      n.chains = 1,
                      quiet = T)
    
    jags <- coda.samples(model = foo, 
                         variable.names = c("mu", "sigma"),
                         n.iter = Ntotal,
                         progress.bar = "none")
  }
  
  stopCluster(cl)
  
  jags_out <- as.matrix(mcmc_list)
  jags_out1 <- as.matrix(mcmc_list[[1]]) # This copy is for extracting the column names
  
  var_names <- c("mu", "sigma")
  # Transform the variables names into the format that can be used by grep()
  # for partial matching
  for (i in 1:length(var_names)) {
    var_names[i] <- paste0("^", var_names[i], "\\[")
  }
  
  # In terms of the partial matching of variable/column names to extract and store posterior samples
  mu_res <- lapply(mcmc_list, function(jags_out) jags_out[,grep(var_names[1], colnames(jags_out1))])
  mu_res <- unlist(mu_res)
  mu <- matrix(mu_res, nrow = Ntotal)[Ntotal,]
  
  sigma_res <- lapply(mcmc_list, function(jags_out) jags_out[,grep(var_names[2], colnames(jags_out1))])
  sigma_res <- unlist(sigma_res)
  sigma <- matrix(sigma_res, nrow = Ntotal)[Ntotal,]
  
  return(list(mu = mu,
              sigma = sigma))
}









dc_melding_jags_sub3 <- function(data,
                                 mu_mean,
                                 mu_sd,
                                 tau_alpha,
                                 tau_beta,
                                 nu_p,
                                 alpha_j,
                                 inits,
                                 N,
                                 Ntotal,
                                 n_chains) {
  n <- length(data)
  
  DATA <- matrix(rep(data, N/n_chains), nrow = N/n_chains, byrow = T)
  
  # Prepare initial values for parallel computation 
  mu_0 <- matrix(inits$mu_0, ncol = n_chains)
  tau_0 <- matrix(inits$tau_0, ncol = n_chains)
  nu_0 <- matrix(inits$nu_0, ncol = n_chains)
  
  model <- "model{
    for (i in 1:N) {
      mu[i] ~ dnorm(mu_mean, 1/mu_var)
      tau[i] ~ dgamma(tau_alpha, tau_beta)
      nu[i] ~ dcat(nu_p[])
      
      for (j in 1:n) {
        logL[i,j] <- loggam((nu[i]+1)/2) - loggam(nu[i]/2) + .5 * log(tau[i]/nu[i]) -
            .5 * (nu[i]+1) * log(1 + tau[i] * pow((Y[i,j]-mu[i]), 2) / nu[i])
        phi[i,j] <- -alpha_j * logL[i,j] + 1000
        zeros[i,j] ~ dpois(phi[i,j])
      }
    }
  }"
  
  
  cl <- makeCluster(n_chains)
  registerDoParallel(cl)
  clusterExport(cl, c("DATA", "mu_0", "tau_0", "nu_0", "mu_mean", "mu_sd",
                      "tau_alpha", "tau_beta", "nu_p", "model"), envir = environment())
  
  mcmc_list <- foreach(i = 1:n_chains, .packages = "rjags") %dopar% {
    foo <- jags.model(textConnection(model),
                      data = list(Y = DATA,
                                  mu_mean = mu_mean,
                                  mu_var = mu_sd^2,
                                  tau_alpha = tau_alpha,
                                  tau_beta = tau_beta,
                                  nu_p = nu_p,
                                  alpha_j = alpha_j,
                                  n = n,
                                  N = N/n_chains,
                                  zeros = matrix(rep(0, n*N/n_chains), ncol = n)),
                      inits = list(mu = mu_0[,i],
                                   tau = tau_0[,i],
                                   nu = nu_0[,i]),
                      n.chains = 1,
                      quiet = T)
    jags <- coda.samples(model = foo, 
                         variable.names = c("mu", "tau", "nu"),
                         n.iter = Ntotal,
                         progress.bar = "none")
  }
  
  stopCluster(cl)
  
  jags_out <- as.matrix(mcmc_list)
  jags_out1 <- as.matrix(mcmc_list[[1]]) # This copy is for extracting the column names
  
  var_names <- c("mu", "tau", "nu")
  # Transform the variables names into the format that can be used by grep()
  # for partial matching
  for (i in 1:length(var_names)) {
    var_names[i] <- paste0("^", var_names[i], "\\[")
  }
  
  # In terms of the partial matching of variable/column names to extract and store posterior samples
  mu_res <- lapply(mcmc_list, function(jags_out) jags_out[,grep(var_names[1], colnames(jags_out1))])
  mu_res <- unlist(mu_res)
  mu <- matrix(mu_res, nrow = Ntotal)[Ntotal,]
  
  tau_res <- lapply(mcmc_list, function(jags_out) jags_out[,grep(var_names[2], colnames(jags_out1))])
  tau_res <- unlist(tau_res)
  tau <- matrix(tau_res, nrow = Ntotal)[Ntotal,]
  
  nu_res <- lapply(mcmc_list, function(jags_out) jags_out[,grep(var_names[3], colnames(jags_out1))])
  nu_res <- unlist(nu_res)
  nu <- matrix(nu_res, nrow = Ntotal)[Ntotal,]
  
  return(list(mu = mu, 
              tau = tau,
              nu = nu))
}







dc_melding_jags_sub5 <- function(data,
                                 sigma1_alpha,
                                 sigma1_beta,
                                 sigma2_alpha,
                                 sigma2_beta,
                                 inits,
                                 alpha_j,
                                 N,
                                 Ntotal,
                                 n_chains) {
  TT <- length(data)
  
  DATA <- matrix(rep(data, N/n_chains), nrow = N/n_chains, byrow = T)
  
  sigma1_0 <- matrix(inits$sigma1_inits, ncol = n_chains)
  sigma2_0 <- matrix(inits$sigma2_inits, ncol = n_chains)
  x_0 <- array(inits$x_inits, dim = c(N/n_chains, n_chains, TT))
  x_0 <- aperm(x_0, c(1,3,2))
  
  inits <- list(sigma1_inits = sigma1_0,
                sigma2_inits = sigma2_0,
                x_inits = x_0)
  
  model <- "model {
    for (i in 1:N) {
      sigma1[i] ~ dgamma(alpha1, beta1)
      sigma2[i] ~ dgamma(alpha2, beta2)
      
      tau1[i] <- pow(sigma1[i], -2)
      tau2[i] <- pow(sigma2[i], -2)
      
      x_0[i] ~ dnorm(1, tau2[i])
      x[i,1] ~ dnorm(x_0[i], tau2[i])
      
      Y[i,1] ~ dnorm(x[i,1], alpha_j * tau1[i])
      
      for (t in 2:TT) {
        x[i,t] ~ dnorm(x[i,t-1], tau2[i])
        Y[i,t] ~ dnorm(x[i,t], alpha_j * tau1[i])
      }
    }
  }"
  
  
  cl <- makeCluster(n_chains)
  registerDoParallel(cl)
  clusterExport(cl, c("DATA", "inits", "model"), envir = environment())
  
  mcmc_list <- foreach(i = 1:n_chains, .packages = "rjags") %dopar% {
    foo <- jags.model(textConnection(model),
                      data = list(Y = DATA,
                                  alpha1 = sigma1_alpha,
                                  beta1 = sigma1_beta,
                                  alpha2 = sigma2_alpha,
                                  beta2 = sigma2_beta,
                                  alpha_j = alpha_j,
                                  N = N/n_chains,
                                  TT = TT),
                      inits = list(sigma1 = inits$sigma1_inits[,i],
                                   sigma2 = inits$sigma2_inits[,i],
                                   x = inits$x_inits[,,i]))
    
    out_jags <- coda.samples(foo,
                             variable.names = c("sigma1", "sigma2", "x"),
                             n.iter = Ntotal)
  }
  
  stopCluster(cl)
  
  jags_out <- as.matrix(mcmc_list)
  jags_out1 <- as.matrix(mcmc_list[[1]])
  
  var_names <- c("sigma1", "sigma2", "x")
  var_names_index <- 1:TT
  
  var_names_full <- rep(NA, TT+2)
  var_names_full[1] <- paste0("^", var_names[1], "\\[")
  var_names_full[2] <- paste0("^", var_names[2], "\\[")
  for (t in 1:TT) {
    var_names_full[t+2] <- paste0("^", var_names[3], "\\[[0-9]+,", var_names_index[t], "\\]")
  }
  
  sigma1_res <- lapply(mcmc_list, function(jags_out) jags_out[,grep(var_names_full[1], colnames(jags_out1))])
  sigma1_res <- unlist(sigma1_res)
  sigma1_res <- matrix(sigma1_res, nrow = Ntotal)[Ntotal,]
  
  sigma2_res <- lapply(mcmc_list, function(jags_out) jags_out[,grep(var_names_full[2], colnames(jags_out1))])
  sigma2_res <- unlist(sigma2_res)
  sigma2_res <- matrix(sigma2_res, nrow = Ntotal)[Ntotal,]
  
  x_res <- matrix(rep(NA, TT*N), ncol = TT)
  for (t in 1:TT) {
    x_temp <- lapply(mcmc_list, function(jags_out) jags_out[,grep(var_names_full[t+2], colnames(jags_out1))])
    x_temp <- unlist(x_temp)
    x_res[,t] <- matrix(x_temp, nrow = Ntotal)[Ntotal,]
  }
  
  return(list(sigma1 = sigma1_res,
              sigma2 = sigma2_res,
              x = x_res))
}








dc_melding_jags_sub11 <- function(data,
                                  mu_mean,
                                  mu_sd,
                                  sigma_alpha,
                                  sigma_beta,
                                  inits,
                                  alpha_j,
                                  N,
                                  Ntotal,
                                  n_chains) {
  n <- length(data)
  
  DATA <- matrix(rep(data, N/n_chains), nrow = N/n_chains, byrow = T)
  
  mu_0 <- matrix(inits$mu_0, ncol = n_chains)
  sigma_0 <- matrix(inits$sigma_0, ncol = n_chains)
  
  model <- "model{
    for (i in 1:N) {
      mu[i] ~ dnorm(mu_mean, 1/mu_var)
      sigma[i] ~ dgamma(sigma_alpha, sigma_beta)
      
      tau[i] <- pow(sigma[i], -2)
      
      for (j in 1:n) {
        Y[i,j] ~ dnorm(mu[i], alpha_j*tau[i])
      }
    }
  }"
  
  
  cl <- makeCluster(n_chains)
  registerDoParallel(cl)
  clusterExport(cl, c("DATA", "mu_0", "sigma_0", "sigma_alpha", "sigma_beta",
                      "model"), 
                envir = environment())
  
  mcmc_list <- foreach(i = 1:n_chains, .packages = "rjags") %dopar% {
    foo <- jags.model(textConnection(model),
                      data = list(Y = DATA,
                                  mu_mean = mu_mean,
                                  mu_var = mu_sd^2,
                                  sigma_alpha = sigma_alpha,
                                  sigma_beta = sigma_beta,
                                  alpha_j = alpha_j,
                                  n = n,
                                  N = N/n_chains),
                      inits = list(mu = mu_0[,i],
                                   sigma = sigma_0[,i]),
                      n.chains = 1,
                      quiet = T)
    
    jags <- coda.samples(model = foo, 
                         variable.names = c("mu", "sigma"),
                         n.iter = Ntotal,
                         progress.bar = "none")
  }
  
  stopCluster(cl)
  
  jags_out <- as.matrix(mcmc_list)
  jags_out1 <- as.matrix(mcmc_list[[1]]) # This copy is for extracting the column names
  
  var_names <- c("mu", "sigma")
  # Transform the variables names into the format that can be used by grep()
  # for partial matching
  for (i in 1:length(var_names)) {
    var_names[i] <- paste0("^", var_names[i], "\\[")
  }
  
  # In terms of the partial matching of variable/column names to extract and store posterior samples
  mu_res <- lapply(mcmc_list, function(jags_out) jags_out[,grep(var_names[1], colnames(jags_out1))])
  mu_res <- unlist(mu_res)
  mu <- matrix(mu_res, nrow = Ntotal)[Ntotal,]
  
  sigma_res <- lapply(mcmc_list, function(jags_out) jags_out[,grep(var_names[2], colnames(jags_out1))])
  sigma_res <- unlist(sigma_res)
  sigma <- matrix(sigma_res, nrow = Ntotal)[Ntotal,]
  
  return(list(mu = mu,
              sigma = sigma))
}









dc_melding_jags_sub9 <- function(data,
                                 mu_mean,
                                 mu_sd,
                                 tau_alpha,
                                 tau_beta,
                                 nu_p,
                                 alpha_j,
                                 inits,
                                 N,
                                 Ntotal,
                                 n_chains) {
  n <- length(data)
  
  DATA <- matrix(rep(data, N/n_chains), nrow = N/n_chains, byrow = T)
  
  # Prepare initial values for parallel computation 
  mu_0 <- matrix(inits$mu_0, ncol = n_chains)
  tau_0 <- matrix(inits$tau_0, ncol = n_chains)
  nu_0 <- matrix(inits$nu_0, ncol = n_chains)
  
  model <- "model{
    for (i in 1:N) {
      mu[i] ~ dnorm(mu_mean, 1/mu_var)
      tau[i] ~ dgamma(tau_alpha, tau_beta)
      nu[i] ~ dcat(nu_p[])
      
      for (j in 1:n) {
        logL[i,j] <- loggam((nu[i]+1)/2) - loggam(nu[i]/2) + .5 * log(tau[i]/nu[i]) -
            .5 * (nu[i]+1) * log(1 + tau[i] * pow((Y[i,j]-mu[i]), 2) / nu[i])
          phi[i,j] <- -alpha_j * logL[i,j] + 1000
          zeros[i,j] ~ dpois(phi[i,j])
      }
    }
  }"
  
  
  cl <- makeCluster(n_chains)
  registerDoParallel(cl)
  clusterExport(cl, c("DATA", "mu_0", "tau_0", "nu_0", "mu_mean", "mu_sd",
                      "tau_alpha", "tau_beta", "nu_p", "model"), envir = environment())
  
  mcmc_list <- foreach(i = 1:n_chains, .packages = "rjags") %dopar% {
    foo <- jags.model(textConnection(model),
                      data = list(Y = DATA,
                                  mu_mean = mu_mean,
                                  mu_var = mu_sd^2,
                                  tau_alpha = tau_alpha,
                                  tau_beta = tau_beta,
                                  nu_p = nu_p,
                                  alpha_j = alpha_j,
                                  n = n,
                                  N = N/n_chains,
                                  zeros = matrix(rep(0, n*N/n_chains), ncol = n)),
                      inits = list(mu = mu_0[,i],
                                   tau = tau_0[,i],
                                   nu = nu_0[,i]),
                      n.chains = 1,
                      quiet = T)
    jags <- coda.samples(model = foo, 
                         variable.names = c("mu", "tau", "nu"),
                         n.iter = Ntotal,
                         progress.bar = "none")
  }
  
  stopCluster(cl)
  
  jags_out <- as.matrix(mcmc_list)
  jags_out1 <- as.matrix(mcmc_list[[1]]) # This copy is for extracting the column names
  
  var_names <- c("mu", "tau", "nu")
  # Transform the variables names into the format that can be used by grep()
  # for partial matching
  for (i in 1:length(var_names)) {
    var_names[i] <- paste0("^", var_names[i], "\\[")
  }
  
  # In terms of the partial matching of variable/column names to extract and store posterior samples
  mu_res <- lapply(mcmc_list, function(jags_out) jags_out[,grep(var_names[1], colnames(jags_out1))])
  mu_res <- unlist(mu_res)
  mu <- matrix(mu_res, nrow = Ntotal)[Ntotal,]
  
  tau_res <- lapply(mcmc_list, function(jags_out) jags_out[,grep(var_names[2], colnames(jags_out1))])
  tau_res <- unlist(tau_res)
  tau <- matrix(tau_res, nrow = Ntotal)[Ntotal,]
  
  nu_res <- lapply(mcmc_list, function(jags_out) jags_out[,grep(var_names[3], colnames(jags_out1))])
  nu_res <- unlist(nu_res)
  nu <- matrix(nu_res, nrow = Ntotal)[Ntotal,]
  
  return(list(mu = mu, 
              tau = tau,
              nu = nu))
}








dc_melding_jags_sub7 <- function(data,
                                 mu1_mu,
                                 mu1_sigma,
                                 mu2_mu,
                                 mu2_sigma,
                                 inits,
                                 alpha_j,
                                 N,
                                 Ntotal,
                                 n_chains) {
  TT <- length(data)
  
  DATA <- matrix(rep(data, N/n_chains), nrow = N/n_chains, byrow = T)
  
  mu1_0 <- matrix(inits$mu1_inits, ncol = n_chains)
  mu2_0 <- matrix(inits$mu2_inits, ncol = n_chains)
  x_0 <- array(inits$x_inits, dim = c(N/n_chains, n_chains, TT))
  x_0 <- aperm(x_0, c(1,3,2))
  
  inits <- list(mu1_inits = mu1_0,
                mu2_inits = mu2_0,
                x_inits = x_0)
  
  model <- "model {
    for (i in 1:N) {
      mu1[i] ~ dnorm(mean1, tau1)
      mu2[i] ~ dnorm(mean2, tau2)
      
      x_0[i] ~ dlnorm(mu1[i], 1/.01)
      x[i,1] ~ dlnorm(mu1[i] + x_0[i], 1/.01)
      
      tau_x[i,1] <- pow(x[i,1], -2) + 1e-10 # Adding 1e-10 for the stability
      
      Y[i,1] ~ dnorm(mu2[i], alpha_j * x[i,1])
      
      for (t in 2:TT) {
        x[i,t] ~ dlnorm(mu1[i] + x[i,t-1], 1/.01)
        tau_x[i,t] <- pow(x[i,t], -2) + 1e-10 # Adding 1e-10 for the stability
        Y[i,t] ~ dnorm(mu2[i], alpha_j * tau_x[i,t])
      }
    }
  }"
  
  
  cl <- makeCluster(n_chains)
  registerDoParallel(cl)
  clusterExport(cl, c("DATA", "inits", "model"), envir = environment())
  
  mcmc_list <- foreach(i = 1:n_chains, .packages = "rjags") %dopar% {
    foo <- jags.model(textConnection(model),
                      data = list(Y = DATA,
                                  mean1 = mu1_mu,
                                  tau1 = 1/mu1_sigma^2,
                                  mean2 = mu2_mu,
                                  tau2 = 1/mu2_sigma^2,
                                  alpha_j = alpha_j,
                                  N = N/n_chains,
                                  TT = TT),
                      inits = list(mu1 = inits$mu1_inits[,i],
                                   mu2 = inits$mu2_inits[,i],
                                   x = inits$x_inits[,,i]))
    
    out_jags <- coda.samples(foo,
                             variable.names = c("mu1", "mu2", "x"),
                             n.iter = Ntotal)
  }
  
  stopCluster(cl)
  
  jags_out <- as.matrix(mcmc_list)
  jags_out1 <- as.matrix(mcmc_list[[1]])
  
  var_names <- c("mu1", "mu2", "x")
  var_names_index <- 1:TT
  
  var_names_full <- rep(NA, TT+2)
  var_names_full[1] <- paste0("^", var_names[1], "\\[")
  var_names_full[2] <- paste0("^", var_names[2], "\\[")
  for (t in 1:TT) {
    var_names_full[t+2] <- paste0("^", var_names[3], "\\[[0-9]+,", var_names_index[t], "\\]")
  }
  
  mu1_res <- lapply(mcmc_list, function(jags_out) jags_out[,grep(var_names_full[1], colnames(jags_out1))])
  mu1_res <- unlist(mu1_res)
  mu1_res <- matrix(mu1_res, nrow = Ntotal)[Ntotal,]
  
  mu2_res <- lapply(mcmc_list, function(jags_out) jags_out[,grep(var_names_full[2], colnames(jags_out1))])
  mu2_res <- unlist(mu2_res)
  mu2_res <- matrix(mu2_res, nrow = Ntotal)[Ntotal,]
  
  x_res <- matrix(rep(NA, TT*N), ncol = TT)
  for (t in 1:TT) {
    x_temp <- lapply(mcmc_list, function(jags_out) jags_out[,grep(var_names_full[t+2], colnames(jags_out1))])
    x_temp <- unlist(x_temp)
    x_res[,t] <- matrix(x_temp, nrow = Ntotal)[Ntotal,]
  }
  
  return(list(mu1 = mu1_res,
              mu2 = mu2_res,
              x = x_res))
}



























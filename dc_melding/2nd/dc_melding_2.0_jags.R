# MCMC for submodels
dc_melding_sub_jags <- function(data,
                                mu_phi,
                                sigma_phi,
                                N,
                                Ntotal) {
  n <- dim(data)[1]
  d <- dim(data)[2]
  
  library(rjags)
  
  model <- "model {
    
    #phi12 ~ dnorm(mu1, tau_phi1)
    #phi23 ~ dnorm(mu2, tau_phi2)
    #tau_phi1 <- pow(sigma1, -2)
    #tau_phi2 <- pow(sigma2, -2)
    phi12 ~ dunif(-1000, 1000)
    phi23 ~ dunif(-1000, 1000)
    
    psi1 ~ dgamma(1, 1)
    psi3 ~ dgamma(1, 1)
    tau1 <- pow(psi1, -2)
    tau3 <- pow(psi3, -2)
    
    for (i in 1:n) {
      y1[i] ~ dnorm(phi12, tau1)
      y3[i] ~ dnorm(phi23, tau3)
    }
  }
  "
  
  parameters <- c("phi12", "phi23",
                  "psi1", "psi3")
  
  burn_in <- 0
  
  steps <- Ntotal
  
  thin <- 1
  
  foo <- jags.model(textConnection(model),
                    data = list(y1=data[,1], y3=data[,2], n=n,
                                mu1=mu_phi[1], mu2=mu_phi[2],
                                sigma1=sigma_phi[1], sigma2=sigma_phi[2]),
                    quiet = T)
  
  update(foo)
  
  mu <- array(rep(NA, N*Ntotal*d), dim = c(N, Ntotal, d))
  sigma <- array(rep(NA, N*Ntotal*d), dim = c(N, Ntotal, d))
  
  for (i in 1:N) {
    jags <- coda.samples(model = foo, variable.names = parameters, n.iter = steps, thin = thin)
    out_bugs <- as.matrix(jags)
    for (j in 1:d) {
      mu[i,,j] <- out_bugs[,j]
      sigma[i,,j] <- out_bugs[,j+d]
    }
  }
  
  min_lim <- log(.Machine$double.xmin)
  w_log <- array(rep(NA, N*Ntotal*d), dim = c(N, Ntotal, d))
  W <- array(rep(NA, N*Ntotal*d), dim = c(N, Ntotal, d))
  for (i in 1:Ntotal) {
    for (j in 1:d) {
      w_log[,i,j] <- apply(cbind(data[,j], mu[,i,j], sigma[,i,j]), 1, function(x) dnorm(x[1], x[2], x[3], log = T))
      w_log[which(w_log < min_lim),i,j] <- min_lim
      W[,i,j] <- exp(w_log[,i,j] - matrixStats::logSumExp(w_log[,i,j]))
    }
  }
  
  return(list(mu=mu, sigma=sigma, w_log=w_log, W=W))
}





# MCMC for the SMC sampler
dc_melding_SMC_jags <- function(data,
                                phi,
                                psi_alpha,
                                psi_beta,
                                psi2_init,
                                alpha,
                                N,
                                Ntotal) {
  library(rjags)
  
  n <- length(data)
  
  model <- "model {
    for (i in 1:n) {
      z[i] ~ dpois(phi[i])
      loglik[i] <- logdensity.norm(y[i], mu, tau)
      phi[i] <- -alpha * loglik[i] + 100
    }
    
    tau <- pow(sigma, -2)
    sigma ~ dgamma(psi_alpha, psi_beta)
  }
  "
  
  parameters <- c("sigma")
  
  burn_in <- 0
  
  steps <- Ntotal
  
  thin <- 1
  
  psi2 <- array(rep(NA, N*Ntotal), dim = c(N, Ntotal))
  
  for (i in 1:N) {
    foo <- jags.model(textConnection(model),
                      data = list(y=data, z=rep(0, n), n=n, mu=phi[i], alpha=alpha,
                                  psi_alpha=psi_alpha, psi_beta=psi_beta),
                      inits = list(sigma=psi2_init[i]),
                      quiet = T)
    update(foo)
    jags <- coda.samples(model = foo, variable.names = parameters, n.iter = steps, thin = thin)
    out_bugs <- as.matrix(jags)
    psi2[i,] <- out_bugs
  }
  
  return(psi2)
}










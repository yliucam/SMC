pmcmc_sub6_kernel_init <- function(data,
                                   inits,
                                   Nparam,
                                   Nstate,
                                   Ntotal,
                                   ncore) {
  Code_init <- nimbleCode({
    mu ~ dnorm(0, sd = 2) # For efficient parallel computation using Nimbles's PMCMC (constants are not changable for compiled models), we treat mu as a parameter here, but we will give it an extremely small random-walk scale.
    sigma ~ dgamma(2, 2) # Same reason as mu
    rho ~ dbeta(9, 1)
    
    u0 ~ dnorm(0, 1)
    x0 <- mu + sigma * u0
    xmean <- mu + rho * (x0 - mu)
    x <- xmean + sigma * u0
    y ~ dnorm(0, sd = exp(x/2))
  })
  
  
  mu_0 <- inits$mu_0
  sigma_0 <- inits$sigma_0
  rho_0 <- inits$rho_0
  theta_0 <- inits$theta_0
  
  inits_collection <- cbind(mu_0, sigma_0, rho_0)
  
  cl <- parallel::makeCluster(ncore)
  clusterExport(cl, c("Code_init", "data", "theta_0", 
                      "inits_collection", "Nparam", "Ntotal"),
                envir = environment())
  clusterEvalQ(cl, library(nimbleSMC))
  
  clusterEvalQ(cl, {
    Model <- nimbleModel(Code_init,
                         data = list(y = NA))
    MCMCconf <- configureMCMC(Model, nodes=NULL)
    MCMCconf$addSampler(c("mu", "sigma", "rho"), 
                        type="RW_block",
                        propCov = diag(c(1e-10, 1e-10, 1)^2))
    MCMCsampler <- buildMCMC(MCMCconf)
    cModel <- compileNimble(Model)
    cMCMC <- compileNimble(MCMCsampler, project=Model)
  })
  
  parMCMC <- function(p, theta_0, inits_collection, y, Ntotal) {
    cModel$x <- theta_0[p]
    cModel$setInits(list(mu = inits_collection[p,1],
                         sigma = inits_collection[p,2],
                         rho = inits_collection[p,3]))
    cModel$y <- y
    cMCMC$run(Ntotal)
    cModel$rho
  }
  
  clusterExport(cl, c("parMCMC"), envir = environment())
  
  res_MCMC <- parLapply(cl, 1:Nparam, 
                        fun = parMCMC,
                        theta_0 = theta_0,
                        inits_collection=inits_collection, 
                        y=data,
                        Ntotal = Ntotal)
  
  stopCluster(cl)
  
  rho <- sapply(res_MCMC, function(x) x)
  
  return(list(rho = rho))
}




pmcmc_sub2_kernel <- function(data, 
                              inits,
                              Nparam,
                              Nstate,
                              Ntotal,
                              ncore,
                              TT,
                              t) {
  Code <- nimbleCode({
    mu ~ dnorm(0, sd = 2) # For efficient parallel computation using Nimbles's PMCMC (constants are not changable for compiled models), we treat mu as a parameter here, but we will give it an extremely small random-walk scale.
    sigma ~ dgamma(2, 2) # Same reason as mu
    rho ~ dbeta(9, 1)
    
    x0 ~ dnorm(mu, sd = sigma)
    xmean[1] <- mu + rho * (x0 - mu)
    x[1] ~ dnorm(xmean[1], sd = sigma)
    for (t in 2:TT) {
      xmean[t] <- mu + rho * (x[t-1] - mu)
      x[t] ~ dnorm(xmean[t], sd = sigma)
    }
    
    
    for (t in 1:TT) {
      yvar[t] <- exp(x[t])
      y[t] ~ dnorm(0, var = yvar[t])
    }
  })
  
  mu_0 <- inits$mu_0
  sigma_0 <- inits$sigma_0
  rho_0 <- inits$rho_0
  # phi_0 <- inits$phi_0
  theta_0 <- inits$theta_0
  
  inits_collection <- cbind(mu_0, sigma_0, rho_0)
  
  cl <- parallel::makeCluster(ncore)
  clusterExport(cl, c("Code", "data", "theta_0", "t", "TT",
                      "inits_collection", "Nparam", "Nstate", "Ntotal"),
                envir = environment())
  clusterEvalQ(cl, library(nimbleSMC))
  
  clusterEvalQ(cl, {
    Model <- nimbleModel(Code,
                         data = list(y = data),
                         constants = list(TT = TT))
    PMCMCconf <- configureMCMC(Model, nodes=NULL)
    PMCMCconf$addSampler(c("mu","sigma","rho"), 
                         type="RW_PF_block",
                         control = list(latents = "x",
                                        pfControl = list(Np = Nstate),
                                        propCov = diag(c(1e-10, 1e-10, 1)^2)))
    PMCMCsampler <- buildMCMC(PMCMCconf)
    cPMCMC <- compileNimble(Model, PMCMCsampler)
  })
  
  parPMCMC <- function(p, theta_0, inits_collection, y, Ntotal, t) {
    if (t == 1) {
      cPMCMC$Model$x[1] <- theta_0[p]
      cPMCMC$Model$x[2] <- NA
    } else {
      cPMCMC$Model$x[1:t] <- theta_0[,p]
      cPMCMC$Model$x[t+1] <- NA
    }
    cPMCMC$Model$setInits(list(mu = inits_collection[p,1],
                               sigma = inits_collection[p,2],
                               rho = inits_collection[p,3]))
    cPMCMC$PMCMCsampler$run(Ntotal)
    cPMCMC$Model$rho
  }
  
  clusterExport(cl, c("parPMCMC"), envir = environment())
  
  res_PMCMC <- parLapply(cl, 1:Nparam, 
                         fun = parPMCMC,
                         theta_0 = theta_0,
                         inits_collection = inits_collection, 
                         y = data,
                         Ntotal = Ntotal,
                         t = t)
  
  stopCluster(cl)
  
  rho <- sapply(res_PMCMC, function(x) x)
  
  return(list(rho = rho))
}
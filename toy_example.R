library(nimbleSMC)
library(rjags)

rho <- .9
sigma_x <- 1
sigma_y <- .2

x <- rep(NA, 100)
y <- rep(NA, 100)
x[1] <- rnorm(1, 0, sigma_x)
y[1] <- x[1] + rnorm(1, 0, sigma_y)
for (i in 2:100) {
  x[i] <- rho * x[i-1] + rnorm(1, 0, sigma_x)
  y[i] <- x[i] + rnorm(1, 0, sigma_y)
}




#---------------regular_MCMC------------------#

mcmc_code <- "model
{
  rho ~ dunif(-1, 1)
  sig2_x_inv ~ dgamma(2, 2)
  sig2_y_inv ~ dgamma(2, 2)
  x[1] ~ dnorm(0, sig2_x_inv)
  v[1] ~ dnorm(0, sig2_y_inv)
  y[1] ~ dnorm(x[1], sig2_y_inv)
  
  for (i in 2:t) {
    u[i] ~ dnorm(0, sig2_x_inv)
    x[i] <- rho * x[i-1] + u[i]  
    v[i] ~ dnorm(0, sig2_y_inv)
    y[i] ~ dnorm(x[i], sig2_y_inv)
  }
}
"

mcmc_model <- jags.model(textConnection(mcmc_code), data = list(y=y, t=100),
                         inits = list(sig2_x_inv=1/1000, sig2_y_inv=1/1000, rho=0))

burn_in <- 5000

update(mcmc_model, burn_in)

out_mcmc <- coda.samples(model = mcmc_model, variable.names = c('x[100]', 'sig2_x_inv', 'sig2_y_inv', 'rho'),
                         n.iter = 10000, thin = 2)
out_mcmc <- as.mcmc(out_mcmc)
traceplot(out_mcmc[,'rho'])
traceplot(out_mcmc[,'sig2_x_inv'])
traceplot(out_mcmc[,'sig2_y_inv'])
traceplot(out_mcmc[,'x[100]'])
hist(out_mcmc[,'rho'])
hist(sqrt(1/out_mcmc[,'sig2_x_inv']**2))
hist(sqrt(1/out_mcmc[,'sig2_y_inv']**2))
hist(out_mcmc[,'x[100]'])





#------------------PMCMC_nimble---------------#

nimble_code <- nimbleCode({
  rho ~ dunif(-1, 1)
  sig2_x_inv ~ dgamma(2, 2)
  sig2_y_inv ~ dgamma(2, 2)
  x[1] ~ dnorm(0, sig2_x_inv)
  v[1] ~ dnorm(0, sig2_y_inv)
  y[1] ~ dnorm(x[1], sig2_y_inv)
  
  for (i in 2:t) {
    x[i] ~ dnorm(rho * x[i-1], sd = sqrt(1/sig2_x_inv))
    y[i] ~ dnorm(x[i], sd = sqrt(1/sig2_y_inv))
  }
})

nimble_model <- nimbleModel(nimble_code,
                            data = list(y=y),
                            constants = list(t=100),
                            inits = list(rho=0, sig2_x_inv=1/1000, sig2_y_inv=1/1000),
                            check=F)

nimble_mcmcconf <- configureMCMC(nimble_model, nodes=NULL)
nimble_mcmcconf$addSampler(target = c('rho', 'sig2_x_inv', 'sig2_y_inv'),
                           type = 'RW_PF_block', control = list(latents='x'))

nimble_MCMC <- buildMCMC(nimble_mcmcconf)
nimble_mcmc <- compileNimble(nimble_model, nimble_MCMC, resetFunctions = T)

nimble_mcmc$nimble_MCMC$run(10000)
















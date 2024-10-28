library(rjags)

hier_bugs <- "model
{
  alpha ~ dunif(.001, 10000)
  
  alpha1 ~ dgamma(alpha, 1)
  alpha2 ~ dgamma(alpha, .5)
  
  for (i in 1:50) {
    p1[i] ~ dbeta(alpha1, 2)
    p2[i] ~ dbeta(alpha2, 10)
    
    for (j in 1:50) {
      y1[j,i] ~ dbin(p1[i], 100)
      y2[j,i] ~ dbin(p2[i], 200)
    }
  }
  
}
" 

data <- list(y1=y1, y2=y2)


parameters <- c('alpha', 'alpha1', 'alpha2')

burn_in <- 30000

steps <- 100000

thin <- 5

foo <- jags.model(textConnection(hier_bugs), data = data)

update(foo, burn_in)

jags <- coda.samples(model=foo, variable.names = parameters, n.iter = steps, thin = thin)



library(posterior, quietly = TRUE)
library(tidybayes)

jags_draws <- as_draws(as.matrix(jags))
summary(jags_draws)


# Double check prior in JAGS model matches the simulated prior

no_data <- list()
prior <- jags.model(textConnection(hier_bugs), data = no_data)
prior_out <- coda.samples(model = prior,
                          variable.names = parameters,
                          n.iter = steps,
                          thin = thin)
prior_out_draws <- as_draws(prior_out)

summary(extract_variable(prior_out_draws, "alpha")) # JAGS prior
summary(rlnorm(10000, 3, 2)) # Prior in R

summary(extract_variable(prior_out_draws, "alpha1")) # JAGS prior
summary(sapply(rlnorm(10000, 3, 2), function(alpha) rgamma(1, alpha, 1))) # Prior in R

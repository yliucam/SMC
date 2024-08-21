library(rjags)

hier_bugs <- "model
{
  alpha0 ~ dlnorm(3, 1/4)

  alpha[1] ~ dgamma(alpha0, 1)
  alpha[2] ~ dgamma(alpha0, .5)

  for (j in 1:4){
    p[1, j] ~ dbeta(alpha[1], 2)
    p[2, j] ~ dbeta(alpha[2], 10)
  }
  
  for (i in 1:20) {
    y11[i] ~ dbin(p[1,1], 100)
    y12[i] ~ dbin(p[1,2], 100)
    y13[i] ~ dbin(p[1,3], 100)
    y14[i] ~ dbin(p[1,4], 100)
    y21[i] ~ dbin(p[2,1], 200)
    y22[i] ~ dbin(p[2,2], 200)
    y23[i] ~ dbin(p[2,3], 200)
    y24[i] ~ dbin(p[2,4], 200)
  }
  
}
" 

data <- list(y11=y11, y12=y12, y13=y13, y14=y14, y21=y21, y22=y22, y23=y23, y24=y24)

parameters <- c('alpha0', 'alpha','p')


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

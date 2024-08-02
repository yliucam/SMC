library(rjags)

hier_bugs <- "model
{
  alpha ~ dlnorm(3, 1/4)
  
  alpha1 ~ dgamma(alpha, 1)
  alpha2 ~ dgamma(alpha, .5)
  
  p11 ~ dbeta(alpha1, 2)
  p12 ~ dbeta(alpha1, 2)
  p13 ~ dbeta(alpha1, 2)
  p14 ~ dbeta(alpha1, 2)
  
  p21 ~ dbeta(alpha2, 10)
  p22 ~ dbeta(alpha2, 10)
  p23 ~ dbeta(alpha2, 10)
  p24 ~ dbeta(alpha2, 10)
  
  for (i in 1:20) {
    y11[i] ~ dbin(p11, 100)
    y12[i] ~ dbin(p12, 100)
    y13[i] ~ dbin(p13, 100)
    y14[i] ~ dbin(p14, 100)
    y21[i] ~ dbin(p21, 200)
    y22[i] ~ dbin(p22, 200)
    y23[i] ~ dbin(p23, 200)
    y24[i] ~ dbin(p24, 200)
  }
  
}
" 

data <- list(y11=y11, y12=y12, y13=y13, y14=y14, y21=y21, y22=y22, y23=y23, y24=y24)

parameters <- c('alpha', 'alpha1', 'alpha2', 'p11', 'p12', 'p13', 'p14', 'p21', 'p22', 'p23', 'p24')

burn_in <- 30000

steps <- 100000

thin <- 5

foo <- jags.model(textConnection(hier_bugs), data = data)

update(foo, burn_in)

out <- coda.samples(model=foo, variable.names = parameters, n.iter = steps, thin = thin)

outmatrix <- as.matrix(out)


# Double check prior in JAGS model matches the simulated prior

library(posterior, quietly = TRUE)
library(tidybayes)

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

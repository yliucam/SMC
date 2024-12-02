source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/1st/dc_melding_1.1_mcmc.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/1st/dc_melding_1.1.R")

set.seed(12345)
#psi3 <- rnorm(1, 1, 100)
psi3 <- rgamma(1, 8, 8)
phi12 <- rnorm(1, psi3, 1)
phi23 <- rnorm(1, psi3, 2)

psi1 <- rgamma(1, 2, 4)
psi2 <- rgamma(1, 9, 6)

y1 <- rnorm(50, phi12, psi1)
y2 <- rnorm(50, phi23, psi2)

debug(dc_melding_1.1)
out <- dc_melding_1.1(data=cbind(y1, y2), 
                      N=2000, 
                      n_leaf=2, 
                      m=10, 
                      alpha=.05, 
                      Ntotal=50)



library(rjags)
model <- " model {
  psi3 ~ dgamma(1, 1)
  
  phi12 ~ dnorm(psi3, 1)
  phi23 ~ dnorm(psi3, 1/4)
  
  psi1 ~ dgamma(1, 1)
  psi2 ~ dgamma(1, 1)
  
  for (i in 1:5) {
    y1[i] ~ dnorm(phi12, 1/psi1^2)
    y2[i] ~ dnorm(phi23, 1/psi2^2)
  }
}
"

parameters <- c("phi12", "phi23", "psi3")

burn_in <- 30000

steps <- 100000

thin <- 5

foo <- jags.model(textConnection(model), data = list(y1=y1, y2=y2))

update(foo, burn_in)

jags <- coda.samples(model=foo, variable.names = parameters, n.iter = steps, thin = thin)

out_bugs <- as.matrix(jags)

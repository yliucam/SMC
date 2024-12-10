library(Rcpp)

sourceCpp("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/1st/dc_melding_1.2_util.cpp")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/1st/dc_melding_1.2.1_mcmc.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/1st/dc_melding_1.2.1.R")

set.seed(12345)
#psi3 <- rnorm(1, 1, 100)
psi0 <- rgamma(1, 10, 50)
phi1 <- rnorm(1, psi0, 1)
phi2 <- rnorm(1, psi0, 10)
phi3 <- rnorm(1, psi0, .1)
phi4 <- rnorm(1, psi0, 8)
phi5 <- rnorm(1, psi0, 2)
phi6 <- rnorm(1, psi0, 100)
phi7 <- rnorm(1, psi0, 50)
phi8 <- rnorm(1, psi0, 20)

psi1 <- rgamma(1, 2, 4)
psi2 <- rgamma(1, 9, 6)
psi3 <- rgamma(1, .1, 5)
psi4 <- rgamma(1, 10, 10)
psi5 <- rgamma(1, .05, 10)
psi6 <- rgamma(1, 5, 10)
psi7 <- rgamma(1, 5, .1)
psi8 <- rgamma(1, .1, .1)

y1 <- rnorm(100, phi1, psi1)
y2 <- rnorm(100, phi2, psi2)
y3 <- rnorm(100, phi3, psi3)
y4 <- rnorm(100, phi4, psi4)
y5 <- rnorm(100, phi5, psi5)
y6 <- rnorm(100, phi6, psi6)
y7 <- rnorm(100, phi7, psi7)
y8 <- rnorm(100, phi8, psi8)

y0 <- rexp(2, psi0)

debug(dc_melding_1.2.1)
out <- dc_melding_1.2.1(data=list(y=cbind(y1, y2, y3, y4, y5, y6, y7, y8), data_0=y0),
                      N=1000, 
                      n_leaf=8, 
                      m=10, 
                      alpha=.05,
                      Ntotal_init=50,
                      Ntotal=50,
                      phi_sigma=c(1, 10, .1, 8, 2, 100, 50, 20))



library(rjags)
model <- " model {
  psi0 ~ dgamma(1, 1)
  
  phi1 ~ dnorm(psi0, 1)
  phi2 ~ dnorm(psi0, 1/100)
  phi3 ~ dnorm(psi0, 1/.01)
  phi4 ~ dnorm(psi0, 1/64)
  phi5 ~ dnorm(psi0, 1/4)
  phi6 ~ dnorm(psi0, 1/10000)
  phi7 ~ dnorm(psi0, 1/250)
  phi8 ~ dnorm(psi0, 1/400)
  
  psi1 ~ dgamma(1, 1)
  psi2 ~ dgamma(1, 1)
  psi3 ~ dgamma(1, 1)
  psi4 ~ dgamma(1, 1)
  psi5 ~ dgamma(1, 1)
  psi6 ~ dgamma(1, 1)
  psi7 ~ dgamma(1, 1)
  psi8 ~ dgamma(1, 1)
  
  for (i in 1:100) {
    y1[i] ~ dnorm(phi1, 1/psi1^2)
    y2[i] ~ dnorm(phi2, 1/psi2^2)
    y3[i] ~ dnorm(phi3, 1/psi3^2)
    y4[i] ~ dnorm(phi4, 1/psi4^2)
    y5[i] ~ dnorm(phi5, 1/psi5^2)
    y6[i] ~ dnorm(phi6, 1/psi6^2)
    y7[i] ~ dnorm(phi7, 1/psi7^2)
    y8[i] ~ dnorm(phi8, 1/psi8^2)
  }
  
  for (i in 1:2) {
    y0[i] ~ dexp(psi0)
  }
}
"

parameters <- c("phi1", "phi2", "phi3", "phi4", "phi5", "phi6", "phi7", "phi8", "psi0")

burn_in <- 30000

steps <- 100000

thin <- 5

foo <- jags.model(textConnection(model), data = list(y1=y1, y2=y2, y3=y3, y4=y4, y5=y5, y6=y6, y7=y7, y8=y8, y0=y0))

update(foo, burn_in)

jags <- coda.samples(model=foo, variable.names = parameters, n.iter = steps, thin = thin)

out_bugs <- as.matrix(jags)

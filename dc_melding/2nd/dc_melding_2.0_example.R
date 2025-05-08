library(Rcpp)

sourceCpp("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/2nd/dc_melding_2.0_util.cpp")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/2nd/dc_melding_2.0_jags.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/2nd/dc_melding_2.0_util.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/2nd/dc_melding_2.0.R")


set.seed(12345)
psi1 <- rgamma(1, 10, 5)
psi3 <- rgamma(1, .1, 1)

psi2 <- rgamma(1, 5, 10)

phi12 <- rnorm(1, 10, 5)
phi23 <- rnorm(1, .1, 10)


y1 <- rnorm(20, phi12, psi1)
y3 <- rnorm(20, phi23, psi3)

y2 <- rnorm(20, phi12+phi23, psi2)

debug(dc_melding_2.0)
out <- dc_melding_2.0(data = list(y1=y1, y2=y2, y3=y3),
                      N = 500,
                      n_sub = 2,
                      m = 100,
                      alpha = .2,
                      lambda = c(1, 1, 1),
                      mu_phi_uni = c(0, 0),
                      sigma_phi_uni = c(1e3, 1e3),
                      mu_phi_mult = c(0, 0),
                      Sigma_phi_mult = cbind(c(1, .8), c(.8, 1)),
                      psi_alpha = 1,
                      psi_beta = 1,
                      Ntotal_sub = 1000,
                      Ntotal = 1000)





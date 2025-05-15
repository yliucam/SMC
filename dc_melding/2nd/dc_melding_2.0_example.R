library(Rcpp)

sourceCpp("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/2nd/dc_melding_2.0_util.cpp")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/2nd/dc_melding_2.0_jags.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/2nd/dc_melding_2.0_smc_mcmc.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/2nd/dc_melding_2.0_util.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/2nd/dc_melding_2.0.R")


set.seed(12345)
psi1 <- rgamma(1, 10, 5)
psi3 <- rgamma(1, 5, 1)

psi2 <- rgamma(1, 5, 10)

phi12 <- rnorm(1, 10, 5)
phi23 <- rnorm(1, .1, 10)


y1 <- rnorm(100, phi12, psi1)
y3 <- rnorm(100, phi23, psi3)

y2 <- rnorm(100, phi12+phi23, psi2)

debug(dc_melding_2.0)
out <- dc_melding_2.0(data = list(y1=y1, y2=y2, y3=y3),
                      N = 5000,
                      n_sub = 2,
                      m = 3000,
                      alpha = .01,
                      lambda = c(1, 1, 1),
                      mu_phi_uni = c(0, 0),
                      sigma_phi_uni = c(1e3, 1e3),
                      mu_phi_mult = c(0, 0),
                      Sigma_phi_mult = cbind(c(1e6, .8*1e6), c(.8*1e6, 1e6)),
                      psi_alpha = 1,
                      psi_beta = 1,
                      Ntotal_sub = 100,
                      Ntotal = 100)

save(out, file = "C:/Users/Yixuan/Documents/codes/SMC/dc_melding/2nd/result.R")





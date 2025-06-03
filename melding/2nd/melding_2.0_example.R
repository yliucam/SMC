library(Rcpp)

sourceCpp("C:/Users/Yixuan/Documents/codes/SMC/melding/2nd/melding_2.0_util.cpp")
source("C:/Users/Yixuan/Documents/codes/SMC/melding/2nd/melding_2.0_jags.R")
source("C:/Users/Yixuan/Documents/codes/SMC/melding/2nd/melding_2.0_util.R")
source("C:/Users/Yixuan/Documents/codes/SMC/melding/2nd/melding_2.0_mcmc.R")
source("C:/Users/Yixuan/Documents/codes/SMC/melding/2nd/melding_2.0.R")


set.seed(12345)
psi1 <- rgamma(1, 10, 5)
psi3 <- rgamma(1, 5, 1)

psi2 <- rgamma(1, 1, 10)

phi12 <- rnorm(1, 10, 5)
phi23 <- rnorm(1, .1, 10)


y1 <- rnorm(100, phi12, psi1)
y3 <- rnorm(100, phi23, psi3)

y2 <- rnorm(100, phi12+phi23, psi2)

debug(melding_2.0)
undebug(melding_2.0)
out_melding <- melding_2.0(data = list(y1=y1, y2=y2, y3=y3),
                           lambda = c(1, 1, 1),
                           Ntotal = 10000,
                           burnin = 5000)

save(out_melding, file = "C:/Users/Yixuan/Documents/codes/SMC/melding/2nd/result_melding.RData")

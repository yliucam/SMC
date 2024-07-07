set.seed(240704)

alpha <- rlnorm(1, 3, 2)

alpha1 <- rgamma(1, alpha, 1)
alpha2 <- rgamma(1, alpha, .5)

p11 <- rbeta(1, alpha1, 2)
p12 <- rbeta(1, alpha1, 10)

p21 <- rbeta(1, alpha2, 2)
p22 <- rbeta(1, alpha2, 10)


y11 <- rbinom(20, 100, p11)
y12 <- rbinom(20, 100, p12)

y21 <- rbinom(20, 200, p21)
y22 <- rbinom(20, 200, p22)


par(mfrow=c(2,2))
hist(y11)
hist(y12)
hist(y21)
hist(y22)

data <- cbind(y11, y12, y21, y22)


library(Rcpp)
library(MASS)
sourceCpp("dc_smc_util.cpp")
source("dc_smc_util.R")
source("dc_smc_leaf_binom.R")

source("dc_smc_algB2_new.R")


gamma_prior <- list(alpha = c(1, 1), beta = c(1, .5))
LN_prior <- list(mu = 0, sigma = 1)

out <- dc_smc_algB2_new(data = data, 
                        n_trial = c(100, 100, 200, 200),
                        N = 5000,
                        nt = 30,
                        gamma_prior = gamma_prior,
                        LN_prior = LN_prior,
                        m = 2,
                        alpha = .2,
                        Ntotal_sub = 5,
                        Ntotal_root = 5)






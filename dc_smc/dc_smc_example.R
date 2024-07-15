set.seed(24001)

alpha <- rlnorm(1, 3, 2)

alpha1 <- rgamma(1, alpha, 1)
alpha2 <- rgamma(1, alpha, .5)

p11 <- rbeta(1, alpha1, 2)
p12 <- rbeta(1, alpha1, 2)
p13 <- rbeta(1, alpha1, 2)
p14 <- rbeta(1, alpha1, 2)

p21 <- rbeta(1, alpha2, 10)
p22 <- rbeta(1, alpha2, 10)
p23 <- rbeta(1, alpha2, 10)
p24 <- rbeta(1, alpha2, 10)

y11 <- rbinom(20, 100, p11)
y12 <- rbinom(20, 100, p12)
y13 <- rbinom(20, 100, p13)
y14 <- rbinom(20, 100, p14)

y21 <- rbinom(20, 200, p21)
y22 <- rbinom(20, 200, p22)
y23 <- rbinom(20, 200, p23)
y24 <- rbinom(20, 200, p24)


par(mfrow=c(2,2))
hist(y11)
hist(y12)
hist(y21)
hist(y22)

data1 <- cbind(y11, y12, y13, y14, y21, y22, y23, y24)


library(Rcpp)
library(MASS)
sourceCpp("dc_smc_util.cpp")
source("dc_smc_util.R")
source("dc_smc_leaf_binom.R")

source("dc_smc_algB2_new.R")


beta_prior <- list(alpha=rep(1, 8), beta=c(rep(2, 4), rep(10, 4)))
n_trial <- c(rep(100, 4), rep(200, 4))
gamma_prior <- list(alpha = c(1, 1), beta = c(1, .5))
LN_prior <- list(mu = 3, sigma = 2)

out_smc <- dc_smc_algB2_new(data = data1, 
                        n_trial = n_trial,
                        N = 1000,
                        nodes_n = c(1, 2, 8),
                        beta_prior = beta_prior,
                        gamma_prior = gamma_prior,
                        LN_prior = LN_prior,
                        m = 2,
                        alpha = .05,
                        Ntotal_sub = 5,
                        Ntotal_root = 5)

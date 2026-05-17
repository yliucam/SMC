library(Rcpp)
library(rjags)
library(doParallel)
library(foreach)
library(msm)
library(matrixStats)

# sourceCpp("dc_melding_owls_util.cpp")
# source("dc_melding_owls_util.R")
# source("merging_loglik.R")
# source("dc_melding_owls_jags.R")
# source("dc_melding_owls.R")
sourceCpp("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/dc_melding_owls_util.cpp")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/dc_melding_owls_util.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/merging_loglik.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/dc_melding_owls_jags.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/dc_melding_owls.R")

count <- read.table("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/data/count.dat")
load("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/results/recap_result.RData")
load("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/results/fecundity_result.RData")

data <- count$V1
alpha0 <- out_recap$alpha[50,1,]
alpha2 <- out_recap$alpha[50,3,]
rho <- out_fec$rho[,50]

N <- 8000
m <- 100

debug(dc_melding_owls_common)
out <- dc_melding_owls_common(data2 = data,
                              N = N,
                              n_sub_var = 3,
                              m = m,
                              alpha = 1,
                              alpha0 = alpha0,
                              alpha2 = alpha2,
                              rho = rho,
                              lambda = c(1, 1, 1),
                              Ntotal = 150,
                              burn_in = 100,
                              n_chains = 1)

save(out, file = "out_owls_dc_melding.RData")


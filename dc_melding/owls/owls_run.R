library(Rcpp)
library(rjags)
library(truncnorm)
library(doParallel)
library(foreach)
library(msm)
library(matrixStats)

sourceCpp("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/dc_melding_owls_util.cpp")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/dc_melding_owls_util.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/dc_melding_owls_jags.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/dc_melding_owls.R")

count <- read.table("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/data/count.dat")
load("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/results/recap_result_N_16000.RData")
load("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/results/fecundity_result_new.RData")

data <- count$V1
alpha0 <- out_recap$alpha[50,1,]
alpha2 <- out_recap$alpha[50,3,]
rho <- out_fec$rho[,11]

debug(dc_melding_owls_common)
out <- dc_melding_owls_common(data2 = data,
                              N = 100,
                              n_sub_var = 3,
                              alpha_incre = .2,
                              alpha0 = alpha0[1:100],
                              alpha2 = alpha2[1:100],
                              rho = rho[1:100],
                              lambda = c(1, 1, 1),
                              Ntotal = 20,
                              n_chains = 10)

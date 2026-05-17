library(Rcpp)
library(rjags)
library(doParallel)
library(foreach)
library(msm)
library(matrixStats)

sourceCpp("dc_melding_owls_util.cpp")
source("dc_melding_owls_util.R")
source("fecundity_model.R")
source("Sys_resamp.R")


fecundity <- read.table("fecundity.dat")

out_fec <- fecundity_model(data = as.matrix(fecundity[1:25,]),
                           lower_prior = 1e-20,
                           upper_prior = 10,
                           N = 12000,
                           Ntotal = 10,
                           n_chains = 50,
                           alpha = .1)

save(out_fec, file = "C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/fecundity_result.RData")
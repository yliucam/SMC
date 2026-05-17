library(Rcpp)
library(rjags)
library(truncnorm)
library(doParallel)
library(foreach)
library(msm)
library(matrixStats)

sourceCpp("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/dc_melding_owls_util.cpp")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/Sys_resamp.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/dc_melding_owls_util.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/dc_melding_owls_jags.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/count_model.R")

count <- read.table("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/data/count.dat")

data <- count$V1

debug(count_model)
out_count <- count_model(data2 = data,
                         N = 500,
                         alpha_incre = .2,
                         Ntotal = 20,
                         n_chains = 10)

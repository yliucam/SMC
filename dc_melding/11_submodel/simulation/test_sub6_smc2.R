library(doMPI)
library(Rcpp)
library(nimbleSMC)
library(rjags)
library(LaplacesDemon)
library(matrixStats)
library(doParallel)
library(foreach)


NUM_REPS <- 50 #500
master.seed <- 20251214
set.seed(master.seed)
seeds <- sample(1:1e7, NUM_REPS)

cl_mpi <- startMPIcluster()
registerDoMPI(cl_mpi)

repN_s <- 1:NUM_REPS
n_s <- 50
sub_s <- 6
PARAMS <- expand.grid(sub = sub_s, rep_N = repN_s)


foreach(jj = 1:nrow(PARAMS)) %dopar% {
  try({
    sourceCpp("dc_melding_11sub_smc2_util.cpp")
    source("Sys_resamp.R")
    source("dc_melding_11sub_smc2_util.R")
    source("dc_melding_11sub_nimble_stage_four.R")
    source("dc_melding_11sub_smc2_stage_four.R")
    
    PATH_PREFIX <-  "/nesi/project/uoa03349/postdoc/dc_melding/11_submodel/results/sub6_smc2/"
    PARAM <- PARAMS[jj,]
    sub_n <- PARAM$sub
    repN <- PARAM$rep_N
    FILENAME <- paste0(PATH_PREFIX, "RES_n", n_s, "_submodel", sub_n, "_smc2_repN", repN, ".RData")
    
    seed <- seeds[repN]
    set.seed(seed)
    phi12 <- rnorm(1, 10, 1)
    phi23 <- rgamma(1, 5, 1)
    phi34  <- rnorm(1, 1, 10)
    phi45 <- rgamma(1, 5, 3)
    
    phi10_11 <- rgamma(1, 4, 2)
    phi9_10 <- rnorm(1, 7, 5)
    phi89 <- rgamma(1, 12, 3)
    phi78 <- rnorm(1, 3, 1)
    
    phi56 <- .178
    phi67 <- -1.024
    
    psi1 <- rgamma(1, 1, 2)
    psi2 <- 5
    psi3 <- 12
    psi4 <- .87
    
    psi11 <- rnorm(1, 4, .8)
    psi10 <- 21
    psi9 <- 2
    psi8 <- .93
    
    psi6 <- .9702
    
    
    y1 <- rnorm(50, phi12, psi1)
    y2 <- rstp(50, mu = phi12, tau = phi23, nu = psi2)
    y3 <- rstp(50, mu = phi34, tau = phi23, nu = psi3)
    
    u4 <- rnorm(10, 0, phi45)
    v4 <- rnorm(10, 0, 1)
    x4 <- integer(10)
    x4_0 <- rnorm(1, phi34, phi45)
    x4[1] <- phi34 + psi4 * (x4_0 - phi34) + u4[1]
    y4 <- integer(10)
    y4[1] <- x4[1] + v4[1]
    for (t in 2:10) {
      x4[t] <- phi34 + psi4 * (x4[t-1] - phi34) + u4[t]
      y4[t] <- x4[t] + v4[t]
    }
    
    x5 <- integer(10)
    x5_0 <- rnorm(1, 1, phi56)
    x5[1] <- rnorm(1, x5_0, phi56)
    y5 <- integer(10)
    y5[1] <- rnorm(1, x5[1], phi45)
    for (t in 2:10) {
      x5[t] <- rnorm(1, x5[t-1], phi56)
      y5[t] <- rnorm(1, x5[t], phi45)
    }
    
    
    y11 <- rnorm(50, psi11, phi10_11)
    y10 <- rstp(50, mu = phi9_10, tau = phi10_11, nu = psi10)
    y9 <- rstp(50, mu = phi9_10, tau = phi89, nu = psi9)
    
    u8 <- rnorm(10, 0, phi89)
    v8 <- rnorm(10, 0, 1)
    x8 <- integer(10)
    x8_0 <- rnorm(1, phi78, phi89)
    x8[1] <- phi78 + psi8 * (x8_0 - phi78) + u8[1]
    y8 <- integer(10)
    y8[1] <- x8[1] + v8[1]
    for (t in 2:10) {
      x8[t] <- phi78 + psi8 * (x8[t-1] - phi78) + u8[t]
      y8[t] <- x8[t] + v8[t]
    }
    
    x7 <- integer(10)
    x7_0 <- exp(rnorm(1, phi67, .1))
    x7[1] <- exp(rnorm(1, x7_0 + phi67, .1))
    y7 <- integer(10)
    y7[1] <- rnorm(1, phi78, x7[1])
    for (t in 2:10) {
      x7[t] <- exp(rnorm(1, x7[t-1] + phi67, .1))
      y7[t] <- rnorm(1, phi78, x7[t])
    }
    
    
    u6 <- rnorm(10, 0, phi56)
    v6 <- rnorm(10, 0, 1)
    x6 <- integer(10)
    x6_0 <- rnorm(1, phi67, phi56)
    x6[1] <- phi67 + psi6 * (x6_0 - phi67) + u6[1]
    y6 <- integer(10)
    y6[1] <- exp(x6[1]/2) * v6[1]
    for (t in 2:10) {
      x6[t] <- phi67 + psi6 * (x6[t-1] - phi67) + u6[t]
      y6[t] <- exp(x6[t]/2) * v6[t]
    }
    
    data <- list(phi12 = phi12, phi23 = phi23, phi34 = phi34, phi45 = phi45,
                 phi56 = phi56, phi67 = phi67, phi78 = phi78, phi89 = phi89,
                 phi9_10 = phi9_10, phi10_11 = phi10_11,
                 psi1 = psi1, psi2 = psi2, psi3 = psi3, psi4 = psi4, 
                 psi6 = psi6, psi8 = psi8, psi9 = psi9, psi10 = psi10, psi11 = psi11,
                 y1 = y1, y2 = y2, y3 = y3, y4 = y4, y5 = y5, y6 = y6, 
                 y7 = y7, y8 = y8, y9 = y9, y10 = y10, y11 = y11)
    
    FILENAME_SUB4 <- paste0("RES_n50_submodel4_repN", repN_s[jj], ".RData")
    FILENAME_SUB8 <- paste0("RES_n50_submodel8_repN", repN_s[jj], ".RData")
    FILENAME_SUB5 <- paste0("RES_n50_submodel5_repN", repN_s[jj], ".RData")
    FILENAME_SUB7 <- paste0("RES_n50_submodel7_repN", repN_s[jj], ".RData")
    FILEADDRESS_SUB4 <- paste0("/nesi/project/uoa03349/postdoc/dc_melding/11_submodel/results/sub4/", FILENAME_SUB4)
    FILEADDRESS_SUB8 <- paste0("/nesi/project/uoa03349/postdoc/dc_melding/11_submodel/results/sub8/", FILENAME_SUB8)
    FILEADDRESS_SUB5 <- paste0("/nesi/project/uoa03349/postdoc/dc_melding/11_submodel/results/sub5/", FILENAME_SUB5)
    FILEADDRESS_SUB7 <- paste0("/nesi/project/uoa03349/postdoc/dc_melding/11_submodel/results/sub7/", FILENAME_SUB7)
    
    load(FILEADDRESS_SUB4)
    A4 <- out$phi_index[,2,6]
    load(FILEADDRESS_SUB8)
    A8 <- out$phi_index[,1,6]
    load(FILEADDRESS_SUB5)
    out5 <- out$sigma2[A4,6]
    load(FILEADDRESS_SUB7)
    out7 <- out$mu1[A8,6]
    
    out <- smc2_sub6(data = y6,
                     phi56 = out5,
                     phi67 = out7,
                     Nstate = 50,
                     Nparam = 10000,
                     n_sub = 2,
                     lambda = c(1/2, 1/2, 1/2),
                     Ntotal = 5,
                     TT = 10,
                     ncore = 10)
    
    save(list = c("out", "data", "seed"), file = FILENAME)
  })
}

closeCluster(cl_mpi)
mpi.quit()




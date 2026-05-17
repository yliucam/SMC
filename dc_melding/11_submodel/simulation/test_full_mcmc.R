library(doMPI)
library(truncnorm)
library(LaplacesDemon)
library(rjags)

NUM_REPS <- 500 #500
master.seed <- 20251214
set.seed(master.seed)
seeds <- sample(1:1e7, NUM_REPS)

cl_mpi <- startMPIcluster()
registerDoMPI(cl_mpi)

repN_s <- 1:NUM_REPS
n_s <- 50

foreach(jj = 1:length(repN_s)) %dopar% {
  try({
    PATH_PREFIX <-  "/nesi/project/uoa03349/postdoc/dc_melding/11_submodel/results/full_mcmc/"
    FILENAME <- paste0(PATH_PREFIX, "RES_n", n_s, "_full_mcmc_repN", repN_s[jj], ".RData")
    
    seed <- seeds[repN_s[jj]]
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
    
    model <- "model{
      phi12 ~ dnorm(0, 1/4)
      phi23 ~ dgamma(2, 2)
      phi34 ~ dnorm(0, 1/4)
      phi45 ~ dgamma(2, 2)
      
      psi1 ~ dgamma(2, 2)
      psi2 ~ dcat(nu_p[])
      psi3 ~ dcat(nu_p[])
      psi4 ~ dbeta(9, 1)
      
      phi10_11 ~ dgamma(2, 2)
      phi9_10 ~ dnorm(0, 1/4)
      phi89 ~ dgamma(12, 2)
      phi78 ~ dnorm(0, 1/4)
      
      psi11 ~ dnorm(0, 1/4)
      psi10 ~ dcat(nu_p[])
      psi9 ~ dcat(nu_p[])
      psi8 ~ dbeta(9, 1)
      
      phi56 ~ dgamma(2, 2)
      phi67 ~ dnorm(0, 1/4)
      
      psi6 ~ dbeta(9, 1)
      
      # model 1, 2 and 3
      for (i in 1:n) {
        y1[i] ~ dnorm(phi12, 1/psi1^2)
        y2[i] ~ dt(phi12, phi23, psi2)
        y3[i] ~ dt(phi34, phi23, psi3)
      }
      
      # model 4
      tau45 <- pow(phi45, -2)
      x4_0 ~ dnorm(phi34, tau45)
      x4mean[1] <- phi34 + psi4 * (x4_0 - phi34)
      x4[1] ~ dnorm(x4mean[1], tau45)
      y4[1] ~ dnorm(x4[1], 1)
      for (t in 2:TT) {
        x4mean[t] <- phi34 + psi4 * (x4[t-1] - phi34)
        x4[t] ~ dnorm(x4mean[t], tau45)
        y4[t] ~ dnorm(x4[t], 1)
      }
      
      # model 5
      tau56 <- pow(phi56, -2)
      x5_0 ~ dnorm(1, tau56)
      x5[1] ~ dnorm(x5_0, tau56)
      y5[1] ~ dnorm(x5[1], tau45)
      for (t in 2:TT) {
        x5[t] ~ dnorm(x5[t-1], tau56)
        y5[t] ~ dnorm(x5[t], tau45)
      }
      
      # model 9, 10 and 11
      for (i in 1:n) {
        y11[i] ~ dnorm(psi11, 1/phi10_11^2)
        y10[i] ~ dt(phi9_10, phi10_11, psi10)
        y9[i] ~ dt(phi9_10, phi89, psi9)
      }
      
      # model 8
      tau89 <- pow(phi89, -2)
      x8_0 ~ dnorm(phi78, tau89)
      x8mean[1] <- phi78 + psi8 * (x8_0 - phi78)
      x8[1] ~ dnorm(x8mean[1], tau89)
      y8[1] ~ dnorm(x8[1], 1)
      for (t in 2:TT) {
        x8mean[t] <- phi78 + psi8 * (x8[t-1] - phi78)
        x8[t] ~ dnorm(x8mean[t], tau89)
        y8[t] ~ dnorm(x8[t], 1)
      }
      
      # model 7
      x7_0 ~ dlnorm(phi67, 1/.01)
      x7[1] ~ dlnorm(x7_0 + phi67, 1/.01)
      tau7[1] <- pow(x7[1], -2) + 1e-10
      y7[1] ~ dnorm(phi78, tau7[1])
      for (t in 2:TT) {
        x7[t] ~ dlnorm(x7[t-1] + phi67, 1/.01)
        tau7[t] <- pow(x7[t], -2) + 1e-10
        y7[t] ~ dnorm(phi78, tau7[t])
      }
      
      # model 6
      x6_0 ~ dnorm(phi67, tau56)
      x6mean[1] <- phi67 + psi6 * (x6_0 - phi67)
      x6[1] ~ dnorm(x6mean[1], tau56)
      y6tau[1] <- 1/exp(x6[1])
      y6[1] ~ dnorm(0, y6tau[1])
      for (t in 2:TT) {
        x6mean[t] <- phi67 + psi6 * (x6[t-1] - phi67)
        x6[t] ~ dnorm(x6mean[t], tau56)
        y6tau[t] <- 1/exp(x6[t])
        y6[t] ~ dnorm(0, y6tau[t])
      }
    }"
    
    data_list <- list(y1 = data$y1, y2 = data$y2, y3 = data$y3, y4 = data$y4, 
                      y5 = data$y5, y6 = data$y6, y7 = data$y7, y8 = data$y8, 
                      y9 = data$y9, y10 = data$y10, y11 = data$y11,
                      n = 50, TT = 10, nu_p = rep(1/30, 30))
    
    parameter <- c("phi12", "phi23", "phi34", "phi45", "phi56", "phi67",
                   "phi78", "phi89", "phi9_10", "phi10_11",
                   "psi1", "psi2", "psi3", "psi4", "psi6",
                   "psi8", "psi9", "psi10", "psi11")
    
    
    foo <- jags.model(textConnection(model), 
                      data = data_list,
                      inits = list(x7 = rep(1e-5, 10)))
    
    update(foo, 5000)
    
    out_jags <- coda.samples(foo,
                             variable.names = parameter,
                             n.iter = 5000)
    out_jags <- as.matrix(out_jags)
    
    save(list = c("out_jags", "data", "seed"), file = FILENAME)
  })
}

closeCluster(cl_mpi)
mpi.quit()













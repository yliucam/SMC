library(Rcpp)
library(rjags)
library(truncnorm)
library(LaplacesDemon)
library(doParallel)
library(foreach)


set.seed(1234)

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





#---------------------------sub1-----------------------------------#
sourceCpp("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_util.cpp")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/Sys_resamp.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_jags_stage_one.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_smc_stage_one.R")


debug(dc_melding_smc_sub1)

out_sub1 <- dc_melding_smc_sub1(data = data$y1,
                                mu_prior = 0,
                                sigma_prior = 2,
                                alpha_prior = 2,
                                beta_prior = 2,
                                N = 10000,
                                Ntotal = 5,
                                n_chains = 10,
                                alpha_incre = .5)



#---------------------------sub3-----------------------------------#
sourceCpp("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_util.cpp")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/Sys_resamp.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_jags_stage_one.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_smc_stage_one.R")


debug(dc_melding_smc_sub3)

out_sub3 <- dc_melding_smc_sub3(data = data$y3,
                                mu_mean_prior = 0,
                                mu_sd_prior = 20,
                                tau_alpha_prior = 10,
                                tau_beta_prior = 2,
                                nu_p_prior = rep(1/20, 20),
                                N = 1000,
                                Ntotal = 5,
                                n_chains = 10,
                                alpha_incre = .5)




#---------------------------sub5-----------------------------------#
sourceCpp("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_util.cpp")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/Sys_resamp.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_jags_stage_one.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_smc_stage_one.R")


debug(dc_melding_smc_sub5)

out_sub5 <- dc_melding_smc_sub5(data = data$y5,
                                alpha1_prior = 2,
                                beta1_prior = 2,
                                alpha2_prior = 2,
                                beta2_prior = 2,
                                N = 10000,
                                Ntotal = 5,
                                n_chains = 10,
                                alpha_incre = .5)





#---------------------------sub11-----------------------------------#
sourceCpp("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_util.cpp")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/Sys_resamp.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_jags_stage_one.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_smc_stage_one.R")


debug(dc_melding_smc_sub11)

out_sub11 <- dc_melding_smc_sub11(data = data$y11,
                                  mu_prior = 0,
                                  sigma_prior = 2,
                                  alpha_prior = 2,
                                  beta_prior = 2,
                                  N = 10000,
                                  Ntotal = 5,
                                  n_chains = 10,
                                  alpha_incre = .5)





#---------------------------sub9-----------------------------------#
sourceCpp("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_util.cpp")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/Sys_resamp.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_jags_stage_one.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_smc_stage_one.R")


debug(dc_melding_smc_sub9)

out_sub9 <- dc_melding_smc_sub9(data = data$y9,
                                mu_mean_prior = 0,
                                mu_sd_prior = 20,
                                tau_alpha_prior = 10,
                                tau_beta_prior = 2,
                                nu_p_prior = rep(1/20, 20),
                                N = 1000,
                                Ntotal = 5,
                                n_chains = 10,
                                alpha_incre = .5)





#---------------------------sub7-----------------------------------#
sourceCpp("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_util.cpp")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/Sys_resamp.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_jags_stage_one.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_smc_stage_one.R")


debug(dc_melding_smc_sub7)

out_sub7 <- dc_melding_smc_sub7(data = y7,
                                mu1_prior = 0,
                                sigma1_prior = 2,
                                mu2_prior = 0,
                                sigma2_prior = 2,
                                N = 500,
                                Ntotal = 5,
                                n_chains = 10,
                                alpha_incre = .5)







#---------------------------sub2-----------------------------------#
sourceCpp("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_util.cpp")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/Sys_resamp.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_util.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_jags_stage_two.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_smc_stage_two.R")

out1 <- out$mu[,6]
out3 <- out$tau[,6]

debug(dc_melding_smc_sub2)

out_sub2 <- dc_melding_smc_sub2(data = data$y2,
                                out_sub1 = out1,
                                out_sub3 = out3,
                                N = 10000,
                                n_sub = 2,
                                alpha_incre = .2,
                                Ntotal = 5,
                                n_chains = 10,
                                mu_prior = 0,
                                sigma_prior = 2,
                                alpha_prior = 10,
                                beta_prior = 2,
                                nu_p_prior = rep(1/20, 20),
                                lambda = c(1/2, 1/2, 1/2))






#---------------------------sub10-----------------------------------#
sourceCpp("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_util.cpp")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/Sys_resamp.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_util.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_jags_stage_two.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_smc_stage_two.R")

out9 <- out$mu[,6]
out11 <- out$sigma[,6]

debug(dc_melding_smc_sub10)

out_sub10 <- dc_melding_smc_sub10(data = data$y10,
                                  out_sub9 = out9,
                                  out_sub11 = out11,
                                  N = 10000,
                                  n_sub = 2,
                                  alpha_incre = .2,
                                  Ntotal = 5,
                                  n_chains = 10,
                                  mu_prior = 0,
                                  sigma_prior = 2,
                                  lower_prior = 1e-10,
                                  upper_prior = 50,
                                  nu_p_prior = rep(1/30, 30),
                                  lambda = c(1/2, 1/2, 1/2))






#---------------------------sub4-----------------------------------#
sourceCpp("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_util.cpp")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/Sys_resamp.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_util.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_jags_stage_three.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_smc_stage_three.R")


A2 <- out$phi_index[,2,6]
out3 <- out$mu[A2,6]
out5 <- out$sigma1[,6]

debug(dc_melding_smc_sub4)

out_sub4 <- dc_melding_smc_sub4(data = data$y4,
                                out_sub3 = out3,
                                out_sub5 = out5,
                                N = 10000,
                                n_sub = 2,
                                alpha_incre = .2,
                                Ntotal = 5,
                                n_chains = 10,
                                mu_prior = 0,
                                sigma_prior = 20,
                                alpha_prior = 2,
                                beta_prior = 2,
                                balpha_prior = 9,
                                bbeta_prior = 1,
                                lambda = c(1/2, 1/2, 1/2))






#---------------------------sub8-----------------------------------#
sourceCpp("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_util.cpp")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/Sys_resamp.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_util.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_jags_stage_three.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_smc_stage_three.R")


A10 <- out$phi_index[,1,6]
out9 <- out$tau[A10,6]
out7 <- out$mu2[,6]

debug(dc_melding_smc_sub8)

out_sub8 <- dc_melding_smc_sub8(data = data$y8,
                                out_sub7 = out7,
                                out_sub9 = out9,
                                N = 10000,
                                n_sub = 2,
                                alpha_incre = .2,
                                Ntotal = 5,
                                n_chains = 10,
                                mu_prior = 0,
                                sigma_prior = 2,
                                alpha_prior = 10,
                                beta_prior = 2,
                                balpha_prior = 9,
                                bbeta_prior = 1,
                                lambda = c(1/2, 1/2, 1/2))






#---------------------------sub6----------------------------------#
sourceCpp("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_util.cpp")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/Sys_resamp.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_util.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_jags_stage_four.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/dc_melding_11sub_smc_stage_four.R")


A4 <- out$phi_index[,2,6]
out5 <- out$sigma2[A4,6]
A8 <- out$phi_index[,1,6]
out7 <- out$mu1[A8,6]

debug(dc_melding_smc_sub6)

out_sub6 <- dc_melding_smc_sub6(data = data$y6,
                                out_sub5 = out5,
                                out_sub7 = out7,
                                N = 10000,
                                n_sub = 2,
                                alpha_incre = .2,
                                Ntotal = 5,
                                n_chains = 10,
                                mu_prior = 0,
                                sigma_prior = 2,
                                alpha_prior = 0,
                                beta_prior = 2,
                                balpha_prior = 9,
                                bbeta_prior = 1,
                                lambda = c(1/2, 1/2, 1/2))












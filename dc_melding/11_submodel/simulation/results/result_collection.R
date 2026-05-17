#-----------------------------------phi56---------------------------------#
mse_phi56_stage_one <- rep(NA, 500)
q05_phi56_stage_one <- rep(NA, 500)
q95_phi56_stage_one <- rep(NA, 500)
cover_phi56_stage_one <- 0
for (i in 1:500) {
  filename_sub5 <- paste0("RES_n50_submodel5_repN", i, ".RData")
  fileaddress_sub5 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub5/",
                             filename_sub5)
  load(fileaddress_sub5)
  phi56 <- data$phi56
  mse_phi56_stage_one[i] <- sum((out$sigma2[,6] - phi56)^2)/10000
  q05_phi56_stage_one[i] <- quantile(out$sigma2[,6], prob = .05)
  q95_phi56_stage_one[i] <- quantile(out$sigma2[,6], prob = .95)
  cover_phi56_stage_one <- cover_phi56_stage_one + 
    (q05_phi56_stage_one[i] <= phi56 && phi56 <= q95_phi56_stage_one[i])
}
mse_phi56_stage_one_avg <- sum(mse_phi56_stage_one)/500
q90_phi56_stage_one_width <- mean(q95_phi56_stage_one - q05_phi56_stage_one)
cover_phi56_stage_one <- cover_phi56_stage_one/500


mse_phi56_dcm <- rep(NA, 500)
q05_phi56_dcm <- rep(NA, 500)
q95_phi56_dcm <- rep(NA, 500)
cover_phi56_dcm <- 0
for (i in 1:500) {
  filename_sub6 <- paste0("RES_n50_submodel6_repN", i, ".RData")
  fileaddress_sub6 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6/",
                             filename_sub6)
  load(fileaddress_sub6)
  phi56 <- data$phi56
  mse_phi56_dcm[i] <- sum((out$phi[,1,21] - phi56)^2)/10000
  q05_phi56_dcm[i] <- quantile(out$phi[,1,21], prob = .05)
  q95_phi56_dcm[i] <- quantile(out$phi[,1,21], prob = .95)
  cover_phi56_dcm <- cover_phi56_dcm + 
    (q05_phi56_dcm[i] <= phi56 && phi56 <= q95_phi56_dcm[i])
}
mse_phi56_dcm_avg <- sum(mse_phi56_dcm)/500
q90_phi56_dcm_width <- mean(q95_phi56_dcm - q05_phi56_dcm)
cover_phi56_dcm <- cover_phi56_dcm/500


mse_phi56_mcmc <- rep(NA, 500)
q05_phi56_mcmc <- rep(NA, 500)
q95_phi56_mcmc <- rep(NA, 500)
cover_phi56_mcmc <- 0
for (i in 1:500) {
  filename_sub6 <- paste0("RES_n50_full_mcmc_repN", i, ".RData")
  fileaddress_sub6 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/full_mcmc/",
                             filename_sub6)
  load(fileaddress_sub6)
  phi56 <- data$phi56
  mse_phi56_mcmc[i] <- sum((out_jags[1:5000,"phi56"] - phi56)^2)/5000
  q05_phi56_mcmc[i] <- quantile(out_jags[1:5000,"phi56"], prob = .05)
  q95_phi56_mcmc[i] <- quantile(out_jags[1:5000,"phi56"], prob = .95)
  cover_phi56_mcmc <- cover_phi56_mcmc + 
    (q05_phi56_mcmc[i] <= phi56 && phi56 <= q95_phi56_mcmc[i])
}
mse_phi56_mcmc_avg <- sum(mse_phi56_mcmc)/500
q90_phi56_mcmc_width <- mean(q95_phi56_mcmc - q05_phi56_mcmc)
cover_phi56_mcmc <- cover_phi56_mcmc/500


mse_phi56_dcm_smc2 <- rep(NA, 500)
q05_phi56_dcm_smc2 <- rep(NA, 500)
q95_phi56_dcm_smc2 <- rep(NA, 500)
cover_phi56_dcm_smc2 <- 0
for (i in 1:500) {
  filename_sub6_smc2 <- paste0("RES_n50_submodel6_smc2_repN", i, ".RData")
  fileaddress_sub6_smc2 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6_smc2/",
                                  filename_sub6_smc2)
  load(fileaddress_sub6_smc2)
  phi56 <- data$phi56
  mse_phi56_dcm_smc2[i] <- sum((out$sigma[,11] - phi56)^2)/10000
  q05_phi56_dcm_smc2[i] <- quantile(out$sigma[,11], prob = .05)
  q95_phi56_dcm_smc2[i] <- quantile(out$sigma[,11], prob = .95)
  cover_phi56_dcm_smc2 <- cover_phi56_dcm_smc2 + 
    (q05_phi56_dcm_smc2[i] <= phi56 && phi56 <= q95_phi56_dcm_smc2[i])
}
mse_phi56_dcm_smc2_avg <- sum(mse_phi56_dcm_smc2)/500
q90_phi56_dcm_smc2_width <- mean(q95_phi56_dcm_smc2 - q05_phi56_dcm_smc2)
cover_phi56_dcm_smc2 <- cover_phi56_dcm_smc2/500





#-----------------------------------phi67---------------------------------#
mse_phi67_stage_one <- rep(NA, 500)
q05_phi67_stage_one <- rep(NA, 500)
q95_phi67_stage_one <- rep(NA, 500)
cover_phi67_stage_one <- 0
for (i in 1:500) {
  filename_sub7 <- paste0("RES_n50_submodel7_repN", i, ".RData")
  fileaddress_sub7 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub7/",
                             filename_sub7)
  load(fileaddress_sub7)
  phi67 <- data$phi67
  mse_phi67_stage_one[i] <- sum((out$mu1[,6] - phi67)^2)/10000
  q05_phi67_stage_one[i] <- quantile(out$mu1[,6], prob = .05)
  q95_phi67_stage_one[i] <- quantile(out$mu1[,6], prob = .95)
  cover_phi67_stage_one <- cover_phi67_stage_one + 
    (q05_phi67_stage_one[i] <= phi67 && phi67 <= q95_phi67_stage_one[i])
}
mse_phi67_stage_one_avg <- sum(mse_phi67_stage_one)/500
q90_phi67_stage_one_width <- mean(q95_phi67_stage_one - q05_phi67_stage_one)
cover_phi67_stage_one <- cover_phi67_stage_one/500


mse_phi67_dcm <- rep(NA, 500)
q05_phi67_dcm <- rep(NA, 500)
q95_phi67_dcm <- rep(NA, 500)
cover_phi67_dcm <- 0
for (i in 1:500) {
  filename_sub6 <- paste0("RES_n50_submodel6_repN", i, ".RData")
  fileaddress_sub6 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6/",
                             filename_sub6)
  load(fileaddress_sub6)
  phi67 <- data$phi67
  mse_phi67_dcm[i] <- sum((out$phi[,2,21] - phi67)^2)/10000
  q05_phi67_dcm[i] <- quantile(out$phi[,2,21], prob = .05)
  q95_phi67_dcm[i] <- quantile(out$phi[,2,21], prob = .95)
  cover_phi67_dcm <- cover_phi67_dcm + 
    (q05_phi67_dcm[i] <= phi67 && phi67 <= q95_phi67_dcm[i])
}
mse_phi67_dcm_avg <- sum(mse_phi67_dcm)/500
q90_phi67_dcm_width <- mean(q95_phi67_dcm - q05_phi67_dcm)
cover_phi67_dcm <- cover_phi67_dcm/500


mse_phi67_mcmc <- rep(NA, 500)
q05_phi67_mcmc <- rep(NA, 500)
q95_phi67_mcmc <- rep(NA, 500)
cover_phi67_mcmc <- 0
for (i in 1:500) {
  filename_sub6 <- paste0("RES_n50_full_mcmc_repN", i, ".RData")
  fileaddress_sub6 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/full_mcmc/",
                             filename_sub6)
  load(fileaddress_sub6)
  phi67 <- data$phi67
  mse_phi67_mcmc[i] <- sum((out_jags[1:5000,"phi67"] - phi67)^2)/5000
  q05_phi67_mcmc[i] <- quantile(out_jags[1:5000,"phi67"], prob = .05)
  q95_phi67_mcmc[i] <- quantile(out_jags[1:5000,"phi67"], prob = .95)
  cover_phi67_mcmc <- cover_phi67_mcmc + 
    (q05_phi67_mcmc[i] <= phi67 && phi67 <= q95_phi67_mcmc[i])
}
mse_phi67_mcmc_avg <- sum(mse_phi67_mcmc)/500
q90_phi67_mcmc_width <- mean(q95_phi67_mcmc - q05_phi67_mcmc)
cover_phi67_mcmc <- cover_phi67_mcmc/500


mse_phi67_dcm_smc2 <- rep(NA, 500)
q05_phi67_dcm_smc2 <- rep(NA, 500)
q95_phi67_dcm_smc2 <- rep(NA, 500)
cover_phi67_dcm_smc2 <- 0
for (i in 1:500) {
  filename_sub6_smc2 <- paste0("RES_n50_submodel6_smc2_repN", i, ".RData")
  fileaddress_sub6_smc2 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6_smc2/",
                                  filename_sub6_smc2)
  load(fileaddress_sub6_smc2)
  phi67 <- data$phi67
  mse_phi67_dcm_smc2[i] <- sum((out$mu[,11] - phi67)^2)/10000
  q05_phi67_dcm_smc2[i] <- quantile(out$mu[,11], prob = .05)
  q95_phi67_dcm_smc2[i] <- quantile(out$mu[,11], prob = .95)
  cover_phi67_dcm_smc2 <- cover_phi67_dcm_smc2 + 
    (q05_phi67_dcm_smc2[i] <= phi67 && phi67 <= q95_phi67_dcm_smc2[i])
}
mse_phi67_dcm_smc2_avg <- sum(mse_phi67_dcm_smc2)/500
q90_phi67_dcm_smc2_width <- mean(q95_phi67_dcm_smc2 - q05_phi67_dcm_smc2)
cover_phi67_dcm_smc2 <- cover_phi67_dcm_smc2/500






#-----------------------------------phi45---------------------------------#
mse_phi45_stage_one <- rep(NA, 500)
q05_phi45_stage_one <- rep(NA, 500)
q95_phi45_stage_one <- rep(NA, 500)
cover_phi45_stage_one <- 0
for (i in 1:500) {
  filename_sub5 <- paste0("RES_n50_submodel5_repN", i, ".RData")
  fileaddress_sub5 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub5/",
                             filename_sub5)
  load(fileaddress_sub5)
  phi45 <- data$phi45
  mse_phi45_stage_one[i] <- sum((out$sigma1[,6] - phi45)^2)/10000
  q05_phi45_stage_one[i] <- quantile(out$sigma1[,6], prob = .05)
  q95_phi45_stage_one[i] <- quantile(out$sigma1[,6], prob = .95)
  cover_phi45_stage_one <- cover_phi45_stage_one + 
    (q05_phi45_stage_one[i] <= phi45 && phi45 <= q95_phi45_stage_one[i])
}
mse_phi45_stage_one_avg <- sum(mse_phi45_stage_one)/500
q90_phi45_stage_one_width <- mean(q95_phi45_stage_one - q05_phi45_stage_one)
cover_phi45_stage_one <- cover_phi45_stage_one/500


mse_phi45_dcm <- rep(NA, 500)
q05_phi45_dcm <- rep(NA, 500)
q95_phi45_dcm <- rep(NA, 500)
cover_phi45_dcm <- 0
for (i in 1:500) {
  filename_sub6 <- paste0("RES_n50_submodel6_repN", i, ".RData")
  fileaddress_sub6 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6/",
                             filename_sub6)
  filename_sub4 <- paste0("RES_n50_submodel4_repN", i, ".RData")
  fileaddress_sub4 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub4/",
                             filename_sub4)
  load(fileaddress_sub6)
  A6 <- out$phi_index[,1,21]
  load(fileaddress_sub4)
  phi45 <- data$phi45
  mse_phi45_dcm[i] <- sum((out$phi[A6,2,6] - phi45)^2)/10000
  q05_phi45_dcm[i] <- quantile(out$phi[A6,2,6], prob = .05)
  q95_phi45_dcm[i] <- quantile(out$phi[A6,2,6], prob = .95)
  cover_phi45_dcm <- cover_phi45_dcm + 
    (q05_phi45_dcm[i] <= phi45 && phi45 <= q95_phi45_dcm[i])
}
mse_phi45_dcm_avg <- sum(mse_phi45_dcm)/500
q90_phi45_dcm_width <- mean(q95_phi45_dcm - q05_phi45_dcm)
cover_phi45_dcm <- cover_phi45_dcm/500


mse_phi45_mcmc <- rep(NA, 500)
q05_phi45_mcmc <- rep(NA, 500)
q95_phi45_mcmc <- rep(NA, 500)
cover_phi45_mcmc <- 0
for (i in 1:500) {
  filename_sub6 <- paste0("RES_n50_full_mcmc_repN", i, ".RData")
  fileaddress_sub6 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/full_mcmc/",
                             filename_sub6)
  load(fileaddress_sub6)
  phi45 <- data$phi45
  mse_phi45_mcmc[i] <- sum((out_jags[1:5000,"phi45"] - phi45)^2)/5000
  q05_phi45_mcmc[i] <- quantile(out_jags[1:5000,"phi45"], prob = .05)
  q95_phi45_mcmc[i] <- quantile(out_jags[1:5000,"phi45"], prob = .95)
  cover_phi45_mcmc <- cover_phi45_mcmc + 
    (q05_phi45_mcmc[i] <= phi45 && phi45 <= q95_phi45_mcmc[i])
}
mse_phi45_mcmc_avg <- sum(mse_phi45_mcmc)/500
q90_phi45_mcmc_width <- mean(q95_phi45_mcmc - q05_phi45_mcmc)
cover_phi45_mcmc <- cover_phi45_mcmc/500


mse_phi45_dcm_smc2 <- rep(NA, 500)
q05_phi45_dcm_smc2 <- rep(NA, 500)
q95_phi45_dcm_smc2 <- rep(NA, 500)
cover_phi45_dcm_smc2 <- 0
for (i in 1:500) {
  filename_sub6_smc2 <- paste0("RES_n50_submodel6_smc2_repN", i, ".RData")
  fileaddress_sub6_smc2 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6_smc2/",
                             filename_sub6_smc2)
  filename_sub4 <- paste0("RES_n50_submodel4_repN", i, ".RData")
  fileaddress_sub4 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub4/",
                             filename_sub4)
  load(fileaddress_sub6_smc2)
  A6 <- out$A5[,11]
  load(fileaddress_sub4)
  phi45 <- data$phi45
  mse_phi45_dcm_smc2[i] <- sum((out$phi[A6,2,6] - phi45)^2)/10000
  q05_phi45_dcm_smc2[i] <- quantile(out$phi[A6,2,6], prob = .05)
  q95_phi45_dcm_smc2[i] <- quantile(out$phi[A6,2,6], prob = .95)
  cover_phi45_dcm_smc2 <- cover_phi45_dcm_smc2 + 
    (q05_phi45_dcm_smc2[i] <= phi45 && phi45 <= q95_phi45_dcm_smc2[i])
}
mse_phi45_dcm_smc2_avg <- sum(mse_phi45_dcm_smc2)/500
q90_phi45_dcm_smc2_width <- mean(q95_phi45_dcm_smc2 - q05_phi45_dcm_smc2)
cover_phi45_dcm_smc2 <- cover_phi45_dcm_smc2/500






#-----------------------------------phi78---------------------------------#

mse_phi78_stage_one <- rep(NA, 500)
q05_phi78_stage_one <- rep(NA, 500)
q95_phi78_stage_one <- rep(NA, 500)
cover_phi78_stage_one <- 0
for (i in 1:500) {
  filename_sub7 <- paste0("RES_n50_submodel7_repN", i, ".RData")
  fileaddress_sub7 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub7/",
                             filename_sub7)
  load(fileaddress_sub7)
  phi78 <- data$phi78
  mse_phi78_stage_one[i] <- sum((out$mu2[,6] - phi78)^2)/10000
  q05_phi78_stage_one[i] <- quantile(out$mu2[,6], prob = .05)
  q95_phi78_stage_one[i] <- quantile(out$mu2[,6], prob = .95)
  cover_phi78_stage_one <- cover_phi78_stage_one + 
    (q05_phi78_stage_one[i] <= phi78 && phi78 <= q95_phi78_stage_one[i])
}
mse_phi78_stage_one_avg <- sum(mse_phi78_stage_one)/500
q90_phi78_stage_one_width <- mean(q95_phi78_stage_one - q05_phi78_stage_one)
cover_phi78_stage_one <- cover_phi78_stage_one/500


mse_phi78_dcm <- rep(NA, 500)
q05_phi78_dcm <- rep(NA, 500)
q95_phi78_dcm <- rep(NA, 500)
cover_phi78_dcm <- 0
for (i in 1:500) {
  filename_sub6 <- paste0("RES_n50_submodel6_repN", i, ".RData")
  fileaddress_sub6 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6/",
                             filename_sub6)
  filename_sub8 <- paste0("RES_n50_submodel8_repN", i, ".RData")
  fileaddress_sub8 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub8/",
                             filename_sub8)
  load(fileaddress_sub6)
  A6 <- out$phi_index[,2,21]
  load(fileaddress_sub8)
  phi78 <- data$phi78
  mse_phi78_dcm[i] <- sum((out$phi[A6,1,6] - phi78)^2)/10000
  q05_phi78_dcm[i] <- quantile(out$phi[A6,1,6], prob = .05)
  q95_phi78_dcm[i] <- quantile(out$phi[A6,1,6], prob = .95)
  cover_phi78_dcm <- cover_phi78_dcm + 
    (q05_phi78_dcm[i] <= phi78 && phi78 <= q95_phi78_dcm[i])
}
mse_phi78_dcm_avg <- sum(mse_phi78_dcm)/500
q90_phi78_dcm_width <- mean(q95_phi78_dcm - q05_phi78_dcm)
cover_phi78_dcm <- cover_phi78_dcm/500


mse_phi78_mcmc <- rep(NA, 500)
q05_phi78_mcmc <- rep(NA, 500)
q95_phi78_mcmc <- rep(NA, 500)
cover_phi78_mcmc <- 0
for (i in 1:500) {
  filename_sub6 <- paste0("RES_n50_full_mcmc_repN", i, ".RData")
  fileaddress_sub6 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/full_mcmc/",
                             filename_sub6)
  load(fileaddress_sub6)
  phi78 <- data$phi78
  mse_phi78_mcmc[i] <- sum((out_jags[1:5000,"phi78"] - phi78)^2)/5000
  q05_phi78_mcmc[i] <- quantile(out_jags[1:5000,"phi78"], prob = .05)
  q95_phi78_mcmc[i] <- quantile(out_jags[1:5000,"phi78"], prob = .95)
  cover_phi78_mcmc <- cover_phi78_mcmc + 
    (q05_phi78_mcmc[i] <= phi78 && phi78 <= q95_phi78_mcmc[i])
}
mse_phi78_mcmc_avg <- sum(mse_phi78_mcmc)/500
q90_phi78_mcmc_width <- mean(q95_phi78_mcmc - q05_phi78_mcmc)
cover_phi78_mcmc <- cover_phi78_mcmc/500


mse_phi78_dcm_smc2 <- rep(NA, 500)
q05_phi78_dcm_smc2 <- rep(NA, 500)
q95_phi78_dcm_smc2 <- rep(NA, 500)
cover_phi78_dcm_smc2 <- 0
for (i in 1:500) {
  filename_sub6_smc2 <- paste0("RES_n50_submodel6_smc2_repN", i, ".RData")
  fileaddress_sub6_smc2 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6_smc2/",
                                  filename_sub6_smc2)
  filename_sub8 <- paste0("RES_n50_submodel8_repN", i, ".RData")
  fileaddress_sub8 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub8/",
                             filename_sub8)
  load(fileaddress_sub6_smc2)
  A6 <- out$A7[,11]
  load(fileaddress_sub8)
  phi78 <- data$phi78
  mse_phi78_dcm_smc2[i] <- sum((out$phi[A6,1,6] - phi78)^2)/10000
  q05_phi78_dcm_smc2[i] <- quantile(out$phi[A6,1,6], prob = .05)
  q95_phi78_dcm_smc2[i] <- quantile(out$phi[A6,1,6], prob = .95)
  cover_phi78_dcm_smc2 <- cover_phi78_dcm_smc2 + 
    (q05_phi78_dcm_smc2[i] <= phi78 && phi78 <= q95_phi78_dcm_smc2[i])
}
mse_phi78_dcm_smc2_avg <- sum(mse_phi78_dcm_smc2)/500
q90_phi78_dcm_smc2_width <- mean(q95_phi78_dcm_smc2 - q05_phi78_dcm_smc2)
cover_phi78_dcm_smc2 <- cover_phi78_dcm_smc2/500






#-----------------------------------phi34---------------------------------#
mse_phi34_stage_one <- rep(NA, 500)
q05_phi34_stage_one <- rep(NA, 500)
q95_phi34_stage_one <- rep(NA, 500)
cover_phi34_stage_one <- 0
for (i in 1:500) {
  filename_sub3 <- paste0("RES_n50_submodel3_repN", i, ".RData")
  fileaddress_sub3 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub3/",
                             filename_sub3)
  load(fileaddress_sub3)
  phi34 <- data$phi34
  mse_phi34_stage_one[i] <- sum((out$mu[,6] - phi34)^2)/10000
  q05_phi34_stage_one[i] <- quantile(out$mu[,6], prob = .05)
  q95_phi34_stage_one[i] <- quantile(out$mu[,6], prob = .95)
  cover_phi34_stage_one <- cover_phi34_stage_one + 
    (q05_phi34_stage_one[i] <= phi34 && phi34 <= q95_phi34_stage_one[i])
}
mse_phi34_stage_one_avg <- sum(mse_phi34_stage_one)/500
q90_phi34_stage_one_width <- mean(q95_phi34_stage_one - q05_phi34_stage_one)
cover_phi34_stage_one <- cover_phi34_stage_one/500


mse_phi34_dcm <- rep(NA, 500)
q05_phi34_dcm <- rep(NA, 500)
q95_phi34_dcm <- rep(NA, 500)
cover_phi34_dcm <- 0
for (i in 1:500) {
  filename_sub6 <- paste0("RES_n50_submodel6_repN", i, ".RData")
  fileaddress_sub6 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6/",
                             filename_sub6)
  filename_sub4 <- paste0("RES_n50_submodel4_repN", i, ".RData")
  fileaddress_sub4 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub4/",
                             filename_sub4)
  load(fileaddress_sub6)
  A6 <- out$phi_index[,1,21]
  load(fileaddress_sub4)
  phi34 <- data$phi34
  mse_phi34_dcm[i] <- sum((out$phi[A6,1,6] - phi34)^2)/10000
  q05_phi34_dcm[i] <- quantile(out$phi[A6,1,6], prob = .05)
  q95_phi34_dcm[i] <- quantile(out$phi[A6,1,6], prob = .95)
  cover_phi34_dcm <- cover_phi34_dcm + 
    (q05_phi34_dcm[i] <= phi34 && phi34 <= q95_phi34_dcm[i])
}
mse_phi34_dcm_avg <- sum(mse_phi34_dcm)/500
q90_phi34_dcm_width <- mean(q95_phi34_dcm - q05_phi34_dcm)
cover_phi34_dcm <- cover_phi34_dcm/500


mse_phi34_mcmc <- rep(NA, 500)
q05_phi34_mcmc <- rep(NA, 500)
q95_phi34_mcmc <- rep(NA, 500)
cover_phi34_mcmc <- 0
for (i in 1:500) {
  filename_sub6 <- paste0("RES_n50_full_mcmc_repN", i, ".RData")
  fileaddress_sub6 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/full_mcmc/",
                             filename_sub6)
  load(fileaddress_sub6)
  phi34 <- data$phi34
  mse_phi34_mcmc[i] <- sum((out_jags[1:5000,"phi34"] - phi34)^2)/5000
  q05_phi34_mcmc[i] <- quantile(out_jags[1:5000,"phi34"], prob = .05)
  q95_phi34_mcmc[i] <- quantile(out_jags[1:5000,"phi34"], prob = .95)
  cover_phi34_mcmc <- cover_phi34_mcmc + 
    (q05_phi34_mcmc[i] <= phi34 && phi34 <= q95_phi34_mcmc[i])
}
mse_phi34_mcmc_avg <- sum(mse_phi34_mcmc)/500
q90_phi34_mcmc_width <- mean(q95_phi34_mcmc - q05_phi34_mcmc)
cover_phi34_mcmc <- cover_phi34_mcmc/500


mse_phi34_dcm_smc2 <- rep(NA, 500)
q05_phi34_dcm_smc2 <- rep(NA, 500)
q95_phi34_dcm_smc2 <- rep(NA, 500)
cover_phi34_dcm_smc2 <- 0
for (i in 1:500) {
  filename_sub6_smc2 <- paste0("RES_n50_submodel6_smc2_repN", i, ".RData")
  fileaddress_sub6_smc2 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6_smc2/",
                                  filename_sub6_smc2)
  filename_sub4 <- paste0("RES_n50_submodel4_repN", i, ".RData")
  fileaddress_sub4 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub4/",
                             filename_sub4)
  load(fileaddress_sub6_smc2)
  A6 <- out$A5[,11]
  load(fileaddress_sub4)
  phi34 <- data$phi34
  mse_phi34_dcm_smc2[i] <- sum((out$phi[A6,1,6] - phi34)^2)/10000
  q05_phi34_dcm_smc2[i] <- quantile(out$phi[A6,1,6], prob = .05)
  q95_phi34_dcm_smc2[i] <- quantile(out$phi[A6,1,6], prob = .95)
  cover_phi34_dcm_smc2 <- cover_phi34_dcm_smc2 + 
    (q05_phi34_dcm_smc2[i] <= phi34 && phi34 <= q95_phi34_dcm_smc2[i])
}
mse_phi34_dcm_smc2_avg <- sum(mse_phi34_dcm_smc2)/500
q90_phi34_dcm_smc2_width <- mean(q95_phi34_dcm_smc2 - q05_phi34_dcm_smc2)
cover_phi34_dcm_smc2 <- cover_phi34_dcm_smc2/500






#-----------------------------------phi89---------------------------------#
mse_phi89_stage_one <- rep(NA, 500)
q05_phi89_stage_one <- rep(NA, 500)
q95_phi89_stage_one <- rep(NA, 500)
cover_phi89_stage_one <- 0
for (i in 1:500) {
  filename_sub9 <- paste0("RES_n50_submodel9_repN", i, ".RData")
  fileaddress_sub9 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub9/",
                             filename_sub9)
  load(fileaddress_sub9)
  phi89 <- data$phi89
  mse_phi89_stage_one[i] <- sum((out$tau[,6] - phi89)^2)/10000
  q05_phi89_stage_one[i] <- quantile(out$tau[,6], prob = .05)
  q95_phi89_stage_one[i] <- quantile(out$tau[,6], prob = .95)
  cover_phi89_stage_one <- cover_phi89_stage_one + 
    (q05_phi89_stage_one[i] <= phi89 && phi89 <= q95_phi89_stage_one[i])
}
mse_phi89_stage_one_avg <- sum(mse_phi89_stage_one)/500
q90_phi89_stage_one_width <- mean(q95_phi89_stage_one - q05_phi89_stage_one)
cover_phi89_stage_one <- cover_phi89_stage_one/500


mse_phi89_dcm <- rep(NA, 500)
q05_phi89_dcm <- rep(NA, 500)
q95_phi89_dcm <- rep(NA, 500)
cover_phi89_dcm <- 0
for (i in 1:500) {
  filename_sub6 <- paste0("RES_n50_submodel6_repN", i, ".RData")
  fileaddress_sub6 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6/",
                             filename_sub6)
  filename_sub8 <- paste0("RES_n50_submodel8_repN", i, ".RData")
  fileaddress_sub8 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub8/",
                             filename_sub8)
  load(fileaddress_sub6)
  A6 <- out$phi_index[,2,21]
  load(fileaddress_sub8)
  phi89 <- data$phi89
  mse_phi89_dcm[i] <- sum((out$phi[A6,2,6] - phi89)^2)/10000
  q05_phi89_dcm[i] <- quantile(out$phi[A6,2,6], prob = .05)
  q95_phi89_dcm[i] <- quantile(out$phi[A6,2,6], prob = .95)
  cover_phi89_dcm <- cover_phi89_dcm + 
    (q05_phi89_dcm[i] <= phi89 && phi89 <= q95_phi89_dcm[i])
}
mse_phi89_dcm_avg <- sum(mse_phi89_dcm)/500
q90_phi89_dcm_width <- mean(q95_phi89_dcm - q05_phi89_dcm)
cover_phi89_dcm <- cover_phi89_dcm/500


mse_phi89_mcmc <- rep(NA, 500)
q05_phi89_mcmc <- rep(NA, 500)
q95_phi89_mcmc <- rep(NA, 500)
cover_phi89_mcmc <- 0
for (i in 1:500) {
  filename_mcmc <- paste0("RES_n50_full_mcmc_repN", i, ".RData")
  fileaddress_mcmc <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/full_mcmc/",
                             filename_mcmc)
  load(fileaddress_mcmc)
  phi89 <- data$phi89
  mse_phi89_mcmc[i] <- sum((out_jags[1:5000,"phi89"] - phi89)^2)/5000
  q05_phi89_mcmc[i] <- quantile(out_jags[1:5000,"phi89"], prob = .05)
  q95_phi89_mcmc[i] <- quantile(out_jags[1:5000,"phi89"], prob = .95)
  cover_phi89_mcmc <- cover_phi89_mcmc + 
    (q05_phi89_mcmc[i] <= phi89 && phi89 <= q95_phi89_mcmc[i])
}
mse_phi89_mcmc_avg <- sum(mse_phi89_mcmc)/500
q90_phi89_mcmc_width <- mean(q95_phi89_mcmc - q05_phi89_mcmc)
cover_phi89_mcmc <- cover_phi89_mcmc/500


mse_phi89_dcm_smc2 <- rep(NA, 500)
q05_phi89_dcm_smc2 <- rep(NA, 500)
q95_phi89_dcm_smc2 <- rep(NA, 500)
cover_phi89_dcm_smc2 <- 0
for (i in 1:500) {
  filename_sub6_smc2 <- paste0("RES_n50_submodel6_smc2_repN", i, ".RData")
  fileaddress_sub6_smc2 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6_smc2/",
                                  filename_sub6_smc2)
  filename_sub8 <- paste0("RES_n50_submodel8_repN", i, ".RData")
  fileaddress_sub8 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub8/",
                             filename_sub8)
  load(fileaddress_sub6_smc2)
  A6 <- out$A7[,11]
  load(fileaddress_sub8)
  phi89 <- data$phi89
  mse_phi89_dcm_smc2[i] <- sum((out$phi[A6,2,6] - phi89)^2)/10000
  q05_phi89_dcm_smc2[i] <- quantile(out$phi[A6,2,6], prob = .05)
  q95_phi89_dcm_smc2[i] <- quantile(out$phi[A6,2,6], prob = .95)
  cover_phi89_dcm_smc2 <- cover_phi89_dcm_smc2 + 
    (q05_phi89_dcm_smc2[i] <= phi89 && phi89 <= q95_phi89_dcm_smc2[i])
}
mse_phi89_dcm_smc2_avg <- sum(mse_phi89_dcm_smc2)/500
q90_phi89_dcm_smc2_width <- mean(q95_phi89_dcm_smc2 - q05_phi89_dcm_smc2)
cover_phi89_dcm_smc2 <- cover_phi89_dcm_smc2/500






#-----------------------------------phi23---------------------------------#
mse_phi23_stage_one <- rep(NA, 500)
q05_phi23_stage_one <- rep(NA, 500)
q95_phi23_stage_one <- rep(NA, 500)
cover_phi23_stage_one <- 0
for (i in 1:500) {
  filename_sub3 <- paste0("RES_n50_submodel3_repN", i, ".RData")
  fileaddress_sub3 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub3/",
                             filename_sub3)
  load(fileaddress_sub3)
  phi23 <- data$phi23
  mse_phi23_stage_one[i] <- sum((out$tau[,6] - phi23)^2)/10000
  q05_phi23_stage_one[i] <- quantile(out$tau[,6], prob = .05)
  q95_phi23_stage_one[i] <- quantile(out$tau[,6], prob = .95)
  cover_phi23_stage_one <- cover_phi23_stage_one + 
    (q05_phi23_stage_one[i] <= phi23 && phi23 <= q95_phi23_stage_one[i])
}
mse_phi23_stage_one_avg <- sum(mse_phi23_stage_one)/500
q90_phi23_stage_one_width <- mean(q95_phi23_stage_one - q05_phi23_stage_one)
cover_phi23_stage_one <- cover_phi23_stage_one/500


mse_phi23_dcm <- rep(NA, 500)
q05_phi23_dcm <- rep(NA, 500)
q95_phi23_dcm <- rep(NA, 500)
cover_phi23_dcm <- 0
for (i in 1:500) {
  filename_sub6 <- paste0("RES_n50_submodel6_repN", i, ".RData")
  fileaddress_sub6 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6/",
                             filename_sub6)
  filename_sub4 <- paste0("RES_n50_submodel4_repN", i, ".RData")
  fileaddress_sub4 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub4/",
                             filename_sub4)
  filename_sub2 <- paste0("RES_n50_submodel2_repN", i, ".RData")
  fileaddress_sub2 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub2/",
                             filename_sub2)
  load(fileaddress_sub6)
  A6 <- out$phi_index[,1,21]
  load(fileaddress_sub4)
  A4 <- out$phi_index[A6,1,6]
  load(fileaddress_sub2)
  phi23 <- data$phi23
  mse_phi23_dcm[i] <- sum((out$phi[A4,2,6] - phi23)^2)/10000
  q05_phi23_dcm[i] <- quantile(out$phi[A4,2,6], prob = .05)
  q95_phi23_dcm[i] <- quantile(out$phi[A4,2,6], prob = .95)
  cover_phi23_dcm <- cover_phi23_dcm + 
    (q05_phi23_dcm[i] <= phi23 && phi23 <= q95_phi23_dcm[i])
}
mse_phi23_dcm_avg <- sum(mse_phi23_dcm)/500
q90_phi23_dcm_width <- mean(q95_phi23_dcm - q05_phi23_dcm)
cover_phi23_dcm <- cover_phi23_dcm/500


mse_phi23_mcmc <- rep(NA, 500)
q05_phi23_mcmc <- rep(NA, 500)
q95_phi23_mcmc <- rep(NA, 500)
cover_phi23_mcmc <- 0
for (i in 1:500) {
  filename_mcmc <- paste0("RES_n50_full_mcmc_repN", i, ".RData")
  fileaddress_mcmc <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/full_mcmc/",
                             filename_mcmc)
  load(fileaddress_mcmc)
  phi23 <- data$phi23
  mse_phi23_mcmc[i] <- sum((out_jags[1:5000,"phi23"] - phi23)^2)/5000
  q05_phi23_mcmc[i] <- quantile(out_jags[1:5000,"phi23"], prob = .05)
  q95_phi23_mcmc[i] <- quantile(out_jags[1:5000,"phi23"], prob = .95)
  cover_phi23_mcmc <- cover_phi23_mcmc + 
    (q05_phi23_mcmc[i] <= phi23 && phi23 <= q95_phi23_mcmc[i])
}
mse_phi23_mcmc_avg <- sum(mse_phi23_mcmc)/500
q90_phi23_mcmc_width <- mean(q95_phi23_mcmc - q05_phi23_mcmc)
cover_phi23_mcmc <- cover_phi23_mcmc/500


mse_phi23_dcm_smc2 <- rep(NA, 500)
q05_phi23_dcm_smc2 <- rep(NA, 500)
q95_phi23_dcm_smc2 <- rep(NA, 500)
cover_phi23_dcm_smc2 <- 0
for (i in 1:500) {
  filename_sub6_smc2 <- paste0("RES_n50_submodel6_smc2_repN", i, ".RData")
  fileaddress_sub6_smc2 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6_smc2/",
                                  filename_sub6_smc2)
  filename_sub4 <- paste0("RES_n50_submodel4_repN", i, ".RData")
  fileaddress_sub4 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub4/",
                             filename_sub4)
  filename_sub2 <- paste0("RES_n50_submodel2_repN", i, ".RData")
  fileaddress_sub2 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub2/",
                             filename_sub2)
  load(fileaddress_sub6_smc2)
  A6 <- out$A5[,11]
  load(fileaddress_sub4)
  A4 <- out$phi_index[A6,1,6]
  load(fileaddress_sub2)
  phi23 <- data$phi23
  mse_phi23_dcm_smc2[i] <- sum((out$phi[A4,2,6] - phi23)^2)/10000
  q05_phi23_dcm_smc2[i] <- quantile(out$phi[A4,2,6], prob = .05)
  q95_phi23_dcm_smc2[i] <- quantile(out$phi[A4,2,6], prob = .95)
  cover_phi23_dcm_smc2 <- cover_phi23_dcm_smc2 + 
    (q05_phi23_dcm_smc2[i] <= phi23 && phi23 <= q95_phi23_dcm_smc2[i])
}
mse_phi23_dcm_smc2_avg <- sum(mse_phi23_dcm_smc2)/500
q90_phi23_dcm_smc2_width <- mean(q95_phi23_dcm_smc2 - q05_phi23_dcm_smc2)
cover_phi23_dcm_smc2 <- cover_phi23_dcm_smc2/500






#-----------------------------------phi12---------------------------------#
mse_phi12_stage_one <- rep(NA, 500)
q05_phi12_stage_one <- rep(NA, 500)
q95_phi12_stage_one <- rep(NA, 500)
cover_phi12_stage_one <- 0
for (i in 1:500) {
  filename_sub1 <- paste0("RES_n50_submodel1_repN", i, ".RData")
  fileaddress_sub1 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub1/",
                             filename_sub1)
  load(fileaddress_sub1)
  phi12 <- data$phi12
  mse_phi12_stage_one[i] <- sum((out$mu[,6] - phi12)^2)/10000
  q05_phi12_stage_one[i] <- quantile(out$mu[,6], prob = .05)
  q95_phi12_stage_one[i] <- quantile(out$mu[,6], prob = .95)
  cover_phi12_stage_one <- cover_phi12_stage_one + 
    (q05_phi12_stage_one[i] <= phi12 && phi12 <= q95_phi12_stage_one[i])
}
mse_phi12_stage_one_avg <- sum(mse_phi12_stage_one)/500
q90_phi12_stage_one_width <- mean(q95_phi12_stage_one - q05_phi12_stage_one)
cover_phi12_stage_one <- cover_phi12_stage_one/500


mse_phi12_dcm <- rep(NA, 500)
q05_phi12_dcm <- rep(NA, 500)
q95_phi12_dcm <- rep(NA, 500)
cover_phi12_dcm <- 0
for (i in 1:500) {
  filename_sub6 <- paste0("RES_n50_submodel6_repN", i, ".RData")
  fileaddress_sub6 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6/",
                             filename_sub6)
  filename_sub4 <- paste0("RES_n50_submodel4_repN", i, ".RData")
  fileaddress_sub4 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub4/",
                             filename_sub4)
  filename_sub2 <- paste0("RES_n50_submodel2_repN", i, ".RData")
  fileaddress_sub2 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub2/",
                             filename_sub2)
  load(fileaddress_sub6)
  A6 <- out$phi_index[,1,21]
  load(fileaddress_sub4)
  A4 <- out$phi_index[A6,1,6]
  load(fileaddress_sub2)
  phi12 <- data$phi12
  mse_phi12_dcm[i] <- sum((out$phi[A4,1,6] - phi12)^2)/10000
  q05_phi12_dcm[i] <- quantile(out$phi[A4,1,6], prob = .05)
  q95_phi12_dcm[i] <- quantile(out$phi[A4,1,6], prob = .95)
  cover_phi12_dcm <- cover_phi12_dcm + 
    (q05_phi12_dcm[i] <= phi12 && phi12 <= q95_phi12_dcm[i])
}
mse_phi12_dcm_avg <- sum(mse_phi12_dcm)/500
q90_phi12_dcm_width <- mean(q95_phi12_dcm - q05_phi12_dcm)
cover_phi12_dcm <- cover_phi12_dcm/500


mse_phi12_mcmc <- rep(NA, 500)
q05_phi12_mcmc <- rep(NA, 500)
q95_phi12_mcmc <- rep(NA, 500)
cover_phi12_mcmc <- 0
for (i in 1:500) {
  filename_mcmc <- paste0("RES_n50_full_mcmc_repN", i, ".RData")
  fileaddress_mcmc <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/full_mcmc/",
                             filename_mcmc)
  load(fileaddress_mcmc)
  phi12 <- data$phi12
  mse_phi12_mcmc[i] <- sum((out_jags[1:5000,"phi12"] - phi12)^2)/5000
  q05_phi12_mcmc[i] <- quantile(out_jags[1:5000,"phi12"], prob = .05)
  q95_phi12_mcmc[i] <- quantile(out_jags[1:5000,"phi12"], prob = .95)
  cover_phi12_mcmc <- cover_phi12_mcmc + 
    (q05_phi12_mcmc[i] <= phi12 && phi12 <= q95_phi12_mcmc[i])
}
mse_phi12_mcmc_avg <- sum(mse_phi12_mcmc)/500
q90_phi12_mcmc_width <- mean(q95_phi12_mcmc - q05_phi12_mcmc)
cover_phi12_mcmc <- cover_phi12_mcmc/500


mse_phi12_dcm_smc2 <- rep(NA, 500)
q05_phi12_dcm_smc2 <- rep(NA, 500)
q95_phi12_dcm_smc2 <- rep(NA, 500)
cover_phi12_dcm_smc2 <- 0
for (i in 1:500) {
  filename_sub6_smc2 <- paste0("RES_n50_submodel6_smc2_repN", i, ".RData")
  fileaddress_sub6_smc2 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6_smc2/",
                                  filename_sub6_smc2)
  filename_sub4 <- paste0("RES_n50_submodel4_repN", i, ".RData")
  fileaddress_sub4 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub4/",
                             filename_sub4)
  filename_sub2 <- paste0("RES_n50_submodel2_repN", i, ".RData")
  fileaddress_sub2 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub2/",
                             filename_sub2)
  load(fileaddress_sub6_smc2)
  A6 <- out$A5[,11]
  load(fileaddress_sub4)
  A4 <- out$phi_index[A6,1,6]
  load(fileaddress_sub2)
  phi12 <- data$phi12
  mse_phi12_dcm_smc2[i] <- sum((out$phi[A4,1,6] - phi12)^2)/10000
  q05_phi12_dcm_smc2[i] <- quantile(out$phi[A4,1,6], prob = .05)
  q95_phi12_dcm_smc2[i] <- quantile(out$phi[A4,1,6], prob = .95)
  cover_phi12_dcm_smc2 <- cover_phi12_dcm_smc2 + 
    (q05_phi12_dcm_smc2[i] <= phi12 && phi12 <= q95_phi12_dcm_smc2[i])
}
mse_phi12_dcm_smc2_avg <- sum(mse_phi12_dcm_smc2)/500
q90_phi12_dcm_smc2_width <- mean(q95_phi12_dcm_smc2 - q05_phi12_dcm_smc2)
cover_phi12_dcm_smc2 <- cover_phi12_dcm_smc2/500






#-----------------------------------phi9_10---------------------------------#
mse_phi9_10_stage_one <- rep(NA, 500)
q05_phi9_10_stage_one <- rep(NA, 500)
q95_phi9_10_stage_one <- rep(NA, 500)
cover_phi9_10_stage_one <- 0
for (i in 1:500) {
  filename_sub9 <- paste0("RES_n50_submodel9_repN", i, ".RData")
  fileaddress_sub9 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub9/",
                             filename_sub9)
  load(fileaddress_sub9)
  phi9_10 <- data$phi9_10
  mse_phi9_10_stage_one[i] <- sum((out$mu[,6] - phi9_10)^2)/10000
  q05_phi9_10_stage_one[i] <- quantile(out$mu[,6], prob = .05)
  q95_phi9_10_stage_one[i] <- quantile(out$mu[,6], prob = .95)
  cover_phi9_10_stage_one <- cover_phi9_10_stage_one + 
    (q05_phi9_10_stage_one[i] <= phi9_10 && phi9_10 <= q95_phi9_10_stage_one[i])
}
mse_phi9_10_stage_one_avg <- sum(mse_phi9_10_stage_one)/500
q90_phi9_10_stage_one_width <- mean(q95_phi9_10_stage_one - q05_phi9_10_stage_one)
cover_phi9_10_stage_one <- cover_phi9_10_stage_one/500


mse_phi9_10_dcm <- rep(NA, 500)
q05_phi9_10_dcm <- rep(NA, 500)
q95_phi9_10_dcm <- rep(NA, 500)
cover_phi9_10_dcm <- 0
for (i in 1:500) {
  filename_sub6 <- paste0("RES_n50_submodel6_repN", i, ".RData")
  fileaddress_sub6 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6/",
                             filename_sub6)
  filename_sub8 <- paste0("RES_n50_submodel8_repN", i, ".RData")
  fileaddress_sub8 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub8/",
                             filename_sub8)
  filename_sub10 <- paste0("RES_n50_submodel10_repN", i, ".RData")
  fileaddress_sub10 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub10/",
                             filename_sub10)
  load(fileaddress_sub6)
  A6 <- out$phi_index[,2,21]
  load(fileaddress_sub8)
  A8 <- out$phi_index[A6,2,6]
  load(fileaddress_sub10)
  phi9_10 <- data$phi9_10
  mse_phi9_10_dcm[i] <- sum((out$phi[A8,1,6] - phi9_10)^2)/10000
  q05_phi9_10_dcm[i] <- quantile(out$phi[A8,1,6], prob = .05)
  q95_phi9_10_dcm[i] <- quantile(out$phi[A8,1,6], prob = .95)
  cover_phi9_10_dcm <- cover_phi9_10_dcm + 
    (q05_phi9_10_dcm[i] <= phi9_10 && phi9_10 <= q95_phi9_10_dcm[i])
}
mse_phi9_10_dcm_avg <- sum(mse_phi9_10_dcm)/500
q90_phi9_10_dcm_width <- mean(q95_phi9_10_dcm - q05_phi9_10_dcm)
cover_phi9_10_dcm <- cover_phi9_10_dcm/500


mse_phi9_10_mcmc <- rep(NA, 500)
q05_phi9_10_mcmc <- rep(NA, 500)
q95_phi9_10_mcmc <- rep(NA, 500)
cover_phi9_10_mcmc <- 0
for (i in 1:500) {
  filename_mcmc <- paste0("RES_n50_full_mcmc_repN", i, ".RData")
  fileaddress_mcmc <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/full_mcmc/",
                             filename_mcmc)
  load(fileaddress_mcmc)
  phi9_10 <- data$phi9_10
  mse_phi9_10_mcmc[i] <- sum((out_jags[1:5000,"phi9_10"] - phi9_10)^2)/5000
  q05_phi9_10_mcmc[i] <- quantile(out_jags[1:5000,"phi9_10"], prob = .05)
  q95_phi9_10_mcmc[i] <- quantile(out_jags[1:5000,"phi9_10"], prob = .95)
  cover_phi9_10_mcmc <- cover_phi9_10_mcmc + 
    (q05_phi9_10_mcmc[i] <= phi9_10 && phi9_10 <= q95_phi9_10_mcmc[i])
}
mse_phi9_10_mcmc_avg <- sum(mse_phi9_10_mcmc)/500
q90_phi9_10_mcmc_width <- mean(q95_phi9_10_mcmc - q05_phi9_10_mcmc)
cover_phi9_10_mcmc <- cover_phi9_10_mcmc/500


mse_phi9_10_dcm_smc2 <- rep(NA, 500)
q05_phi9_10_dcm_smc2 <- rep(NA, 500)
q95_phi9_10_dcm_smc2 <- rep(NA, 500)
cover_phi9_10_dcm_smc2 <- 0
for (i in 1:500) {
  filename_sub6_smc2 <- paste0("RES_n50_submodel6_smc2_repN", i, ".RData")
  fileaddress_sub6_smc2 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6_smc2/",
                                  filename_sub6_smc2)
  filename_sub8 <- paste0("RES_n50_submodel8_repN", i, ".RData")
  fileaddress_sub8 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub8/",
                             filename_sub8)
  filename_sub10 <- paste0("RES_n50_submodel10_repN", i, ".RData")
  fileaddress_sub10 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub10/",
                              filename_sub10)
  load(fileaddress_sub6_smc2)
  A6 <- out$A7[,11]
  load(fileaddress_sub8)
  A8 <- out$phi_index[A6,2,6]
  load(fileaddress_sub10)
  phi9_10 <- data$phi9_10
  mse_phi9_10_dcm_smc2[i] <- sum((out$phi[A8,1,6] - phi9_10)^2)/10000
  q05_phi9_10_dcm_smc2[i] <- quantile(out$phi[A8,1,6], prob = .05)
  q95_phi9_10_dcm_smc2[i] <- quantile(out$phi[A8,1,6], prob = .95)
  cover_phi9_10_dcm_smc2 <- cover_phi9_10_dcm_smc2 + 
    (q05_phi9_10_dcm_smc2[i] <= phi9_10 && phi9_10 <= q95_phi9_10_dcm_smc2[i])
}
mse_phi9_10_dcm_smc2_avg <- sum(mse_phi9_10_dcm_smc2)/500
q90_phi9_10_dcm_smc2_width <- mean(q95_phi9_10_dcm_smc2 - q05_phi9_10_dcm_smc2)
cover_phi9_10_dcm_smc2 <- cover_phi9_10_dcm_smc2/500






#-----------------------------------phi10_11---------------------------------#
mse_phi10_11_stage_one <- rep(NA, 500)
q05_phi10_11_stage_one <- rep(NA, 500)
q95_phi10_11_stage_one <- rep(NA, 500)
cover_phi10_11_stage_one <- 0
for (i in 1:500) {
  filename_sub11 <- paste0("RES_n50_submodel11_repN", i, ".RData")
  fileaddress_sub11 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub11/",
                             filename_sub11)
  load(fileaddress_sub11)
  phi10_11 <- data$phi10_11
  mse_phi10_11_stage_one[i] <- sum((out$sigma[,6] - phi10_11)^2)/10000
  q05_phi10_11_stage_one[i] <- quantile(out$sigma[,6], prob = .05)
  q95_phi10_11_stage_one[i] <- quantile(out$sigma[,6], prob = .95)
  cover_phi10_11_stage_one <- cover_phi10_11_stage_one + 
    (q05_phi10_11_stage_one[i] <= phi10_11 && phi10_11 <= q95_phi10_11_stage_one[i])
}
mse_phi10_11_stage_one_avg <- sum(mse_phi10_11_stage_one)/500
q90_phi10_11_stage_one_width <- mean(q95_phi10_11_stage_one - q05_phi10_11_stage_one)
cover_phi10_11_stage_one <- cover_phi10_11_stage_one/500


mse_phi10_11_dcm <- rep(NA, 500)
q05_phi10_11_dcm <- rep(NA, 500)
q95_phi10_11_dcm <- rep(NA, 500)
cover_phi10_11_dcm <- 0
for (i in 1:500) {
  filename_sub6 <- paste0("RES_n50_submodel6_repN", i, ".RData")
  fileaddress_sub6 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6/",
                             filename_sub6)
  filename_sub8 <- paste0("RES_n50_submodel8_repN", i, ".RData")
  fileaddress_sub8 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub8/",
                             filename_sub8)
  filename_sub10 <- paste0("RES_n50_submodel10_repN", i, ".RData")
  fileaddress_sub10 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub10/",
                              filename_sub10)
  load(fileaddress_sub6)
  A6 <- out$phi_index[,2,21]
  load(fileaddress_sub8)
  A8 <- out$phi_index[A6,2,6]
  load(fileaddress_sub10)
  phi10_11 <- data$phi10_11
  mse_phi10_11_dcm[i] <- sum((out$phi[A8,2,6] - phi10_11)^2)/10000
  q05_phi10_11_dcm[i] <- quantile(out$phi[A8,2,6], prob = .05)
  q95_phi10_11_dcm[i] <- quantile(out$phi[A8,2,6], prob = .95)
  cover_phi10_11_dcm <- cover_phi10_11_dcm + 
    (q05_phi10_11_dcm[i] <= phi10_11 && phi10_11 <= q95_phi10_11_dcm[i])
}
mse_phi10_11_dcm_avg <- sum(mse_phi10_11_dcm)/500
q90_phi10_11_dcm_width <- mean(q95_phi10_11_dcm - q05_phi10_11_dcm)
cover_phi10_11_dcm <- cover_phi10_11_dcm/500


mse_phi10_11_mcmc <- rep(NA, 500)
q05_phi10_11_mcmc <- rep(NA, 500)
q95_phi10_11_mcmc <- rep(NA, 500)
cover_phi10_11_mcmc <- 0
for (i in 1:500) {
  filename_mcmc <- paste0("RES_n50_full_mcmc_repN", i, ".RData")
  fileaddress_mcmc <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/full_mcmc/",
                             filename_mcmc)
  load(fileaddress_mcmc)
  phi10_11 <- data$phi10_11
  mse_phi10_11_mcmc[i] <- sum((out_jags[1:5000,"phi10_11"] - phi10_11)^2)/5000
  q05_phi10_11_mcmc[i] <- quantile(out_jags[1:5000,"phi10_11"], prob = .05)
  q95_phi10_11_mcmc[i] <- quantile(out_jags[1:5000,"phi10_11"], prob = .95)
  cover_phi10_11_mcmc <- cover_phi10_11_mcmc + 
    (q05_phi10_11_mcmc[i] <= phi10_11 && phi10_11 <= q95_phi10_11_mcmc[i])
}
mse_phi10_11_mcmc_avg <- sum(mse_phi10_11_mcmc)/500
q90_phi10_11_mcmc_width <- mean(q95_phi10_11_mcmc - q05_phi10_11_mcmc)
cover_phi10_11_mcmc <- cover_phi10_11_mcmc/500


mse_phi10_11_dcm_smc2 <- rep(NA, 500)
q05_phi10_11_dcm_smc2 <- rep(NA, 500)
q95_phi10_11_dcm_smc2 <- rep(NA, 500)
cover_phi10_11_dcm_smc2 <- 0
for (i in 1:500) {
  filename_sub6_smc2 <- paste0("RES_n50_submodel6_smc2_repN", i, ".RData")
  fileaddress_sub6_smc2 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6_smc2/",
                                  filename_sub6_smc2)
  filename_sub8 <- paste0("RES_n50_submodel8_repN", i, ".RData")
  fileaddress_sub8 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub8/",
                             filename_sub8)
  filename_sub10 <- paste0("RES_n50_submodel10_repN", i, ".RData")
  fileaddress_sub10 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub10/",
                              filename_sub10)
  load(fileaddress_sub6_smc2)
  A6 <- out$A7[,11]
  load(fileaddress_sub8)
  A8 <- out$phi_index[A6,2,6]
  load(fileaddress_sub10)
  phi10_11 <- data$phi10_11
  mse_phi10_11_dcm_smc2[i] <- sum((out$phi[A8,2,6] - phi10_11)^2)/10000
  q05_phi10_11_dcm_smc2[i] <- quantile(out$phi[A8,2,6], prob = .05)
  q95_phi10_11_dcm_smc2[i] <- quantile(out$phi[A8,2,6], prob = .95)
  cover_phi10_11_dcm_smc2 <- cover_phi10_11_dcm_smc2 + 
    (q05_phi10_11_dcm_smc2[i] <= phi10_11 && phi10_11 <= q95_phi10_11_dcm_smc2[i])
}
mse_phi10_11_dcm_smc2_avg <- sum(mse_phi10_11_dcm_smc2)/500
q90_phi10_11_dcm_smc2_width <- mean(q95_phi10_11_dcm_smc2 - q05_phi10_11_dcm_smc2)
cover_phi10_11_dcm_smc2 <- cover_phi10_11_dcm_smc2/500


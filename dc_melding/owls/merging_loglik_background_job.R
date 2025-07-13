library(nimble)

count <- read.table("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/data/count.dat")

load("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/recap_result.RData")
load("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/fecundity_result.RData")

out_recap <- out
out_fec <- out_fec 


merging_model_code <- nimbleCode({
  ##------------------------------------------------------------##
  ## alpha = (alpha0, alpha2)
  ##------------------------------------------------------------##
  logit(deltaJF) <- alpha0
  logit(deltaAF) <- alpha0 + alpha2
  log(eta_t) <- 0 # mean of the prior of alpha6
  
  xJ[1] <- 25
  sur[1] <- 12
  imm[1] <- 12
  
  for (tt in 2:ti) {
    rateJ[tt-1] <- .5 * rho * deltaJF * x[tt-1]
    xJ[tt] <- rateJ[tt-1]
    sur[tt] <- deltaAF * x[tt-1]
    rate_imm[tt-1] <- x[tt-1] * eta_t
    imm[tt] <- rate_imm[tt-1]
  }
  
  for (tt in 1:ti) {
    x[tt] <- xJ[tt] + sur[tt] + imm[tt]
    y[tt] ~ dpois(x[tt])
  }
  
})


log_lik <- rep(NA, 8000)

alpha0 <- out_recap$alpha[50,1,]
alpha2 <- out_recap$alpha[50,3,]

for (i in 1:8000) {
  invisible(capture.output(merging_model <- nimbleModel(code = merging_model_code,
                                                        constants = list(ti=26),
                                                        data = list(y=as.numeric(unlist(count)),
                                                                    alpha=out_recap$alpha[50,c(1,3),i],
                                                                    rho=out_fec$rho[i, 50]))))
  log_lik[i] <- merging_model$calculate()
}




save(log_lik, "C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/log_lik.RData")







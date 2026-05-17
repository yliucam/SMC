library(rjags)

model <- "model {
    for (i in 1:50) {
      p[i] <- 1/50
    }
    
    logit(deltaJF) <- alpha0
    logit(deltaAF) <- alpha0 + alpha2
    
    xJ[1] ~ dcat(p[1:50])
    sur[1] ~ dcat(p[1:50])
    imm[1] ~ dcat(p[1:50])
    
    alpha0 ~ dnorm(0, 1/4) T(-10, 10)
    alpha2 ~ dnorm(0, 1/4) T(-10, 10)
    alpha6 ~ dnorm(0, 1/4) T(-10, 10)
    eta_t <- exp(alpha6)
    
    rho ~ dunif(0, 10)
    
    for (tt in 2:t) {
      rateJ[tt-1] <- .5 * rho * deltaJF * x[tt-1]
      xJ[tt] ~ dpois(rateJ[tt-1])
      sur[tt] ~ dbin(deltaAF, x[tt-1])
      rate_imm[tt-1] <- x[tt-1] * eta_t
      imm[tt] ~ dpois(rate_imm[tt-1])
    }
    
    for (tt in 1:t) {
      x[tt] <- xJ[tt] + sur[tt] + imm[tt]
      y[tt] ~ dpois(x[tt])
    }
  }"

foo <- jags.model(textConnection(model),
                  data = list(y = data,
                              t = 26),
                  inits = list(xJ = rep(10, 26),
                               sur = rep(10, 26),
                               imm = rep(10, 26)),
                  n.chains = 6)

update(foo, 6000)

out_jags <- coda.samples(foo,
                         variable.names = c("alpha0", "alpha2", "rho", "alpha6"),
                         n.iter = 16000)

out_jags <- as.matrix(out_jags)

save(out_jags, file = "C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/results/submodel2_mcmc.RData")


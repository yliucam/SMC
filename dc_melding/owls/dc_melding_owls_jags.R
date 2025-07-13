library(rjags)


model <- "model {

  logit(deltaJF) <- alpha0
  logit(deltaAF) <- alpha0 + alpha2
  
  for (i in 1:50) {
    p[i] <- 1/50
  }
  
  xJ[1] ~ dcat(p[1:50])
  sur[1] ~ dcat(p[1:50])
  imm[1] ~ dcat(p[1:50])
  
  alpha6 ~ dnorm(0, 1/4) T(-10, 10)
  eta_t <- exp(alpha6)
  
  for (tt in 2:t) {
    rateJ[tt-1] <- .5 * rho * deltaJF * x[tt-1]
    xJ[tt] ~ dpois(rateJ[tt-1])
    sur[tt] ~ dbin(deltaJF, x[tt-1])
    rate_imm[tt-1] <- x[tt-1] * eta_t
    imm[tt] ~ dpois(rate_imm[tt-1])
  }
  
  for (tt in 1:t) {
    x[tt] <- xJ[tt] + sur[tt] + imm[tt]
    y[tt] ~ dpois(x[tt])
  }

}"


foo <- jags.model(textConnection(model),
                  data = list(y = count,
                              t = length(count),
                              alpha0 = -2.5,
                              alpha2 = 2.5,
                              rho = 2.5))

update(foo, 3000)
jags <- coda.samples(model = foo,
                     variable.names = c("x", "alpha6"),
                     n.iter = 10000)
jags_out <- as.matrix(jags)




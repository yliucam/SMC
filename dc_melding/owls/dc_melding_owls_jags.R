library(rjags)
library(doParallel)
library(foreach)
library(peakRAM)


N <- 8000
n_chains <- 10

count <- count$V1

alpha0 <- out_recap$alpha[50,1,1:N]
alpha2 <- out_recap$alpha[50,3,1:N]

rho <- out_fec$rho[1:N,50]

DATA <- matrix(rep(count, N/n_chains), ncol = N/n_chains)
alpha0 <- matrix(alpha0, ncol = n_chains)
alpha2 <- matrix(alpha2, ncol = n_chains)
rho <- matrix(rho, ncol = n_chains)

model <- "model {
  
  for (i in 1:50) {
    p[i] <- 1/50
  }
  
  for (i in 1:N) {
   logit(deltaJF[i]) <- alpha0[i]
   logit(deltaAF[i]) <- alpha0[i] + alpha2[i]  
  
   xJ[1,i] ~ dcat(p[1:50])
   sur[1,i] ~ dcat(p[1:50])
   imm[1,i] ~ dcat(p[1:50])
   
   alpha6[i] ~ dnorm(0, 1/4) T(-10, 10)
   eta_t[i] <- exp(alpha6[i])
   
   for (tt in 2:t) {
    rateJ[tt-1,i] <- .5 * rho[i] * deltaJF[i] * x[tt-1,i]
    xJ[tt,i] ~ dpois(rateJ[tt-1,i])
    sur[tt,i] ~ dbin(deltaJF[i], x[tt-1,i])
    rate_imm[tt-1,i] <- x[tt-1,i] * eta_t[i]
    imm[tt,i] ~ dpois(rate_imm[tt-1,i])
   }
   
   for (tt in 1:t) {
    x[tt,i] <- xJ[tt,i] + sur[tt,i] + imm[tt,i]
    logL[tt,i] <- logdensity.pois(y[tt,i], x[tt,i])
    phi[tt,i] <- -alpha * logL[tt,i] + 100
    zeros[tt,i] ~ dpois(phi[tt,i])
   }
  }

}"


tim_init <- Sys.time()

cl <- makeCluster(n_chains)
registerDoParallel(cl)

mcmc_list <- foreach(i=1:n_chains, .packages = c("peakRAM","rjags")) %dopar% {
  peakRAM({
    foo <- jags.model(textConnection(model),
                      data = list(y = DATA,
                                  t = dim(DATA)[1],
                                  N = N/n_chains,
                                  alpha0 = alpha0[,i],
                                  alpha2 = alpha2[,i],
                                  rho = rho[,i],
                                  alpha = .1,
                                  zeros = matrix(rep(0, 26*N/n_chains), ncol=N/n_chains)))
    
    update(foo, 100)
    jags <- coda.samples(model = foo,
                         variable.names = c("eta_t"),
                         n.iter = 100)
  })
  
}

stopCluster(cl)

tim_end <- Sys.time()
tim_run <- tim_end - tim_init


jags_out <- as.matrix(mcmc_list)
jags_out <- unlist(jags_out)

imm_rate <- matrix(jags_out, nrow = 100)



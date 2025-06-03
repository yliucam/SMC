library(rjags)

melding_jags <- function(data,
                         Ntotal,
                         burnin) {
  n <- dim(data)[1]
  d <- dim(data)[2]
  
  model <- "model {
    phi12 ~ dunif(-1000, 1000)
    phi23 ~ dunif(-1000, 1000)
    
    psi1 ~ dgamma(1, 1)
    psi3 ~ dgamma(1, 1)
    tau1 <- pow(psi1, -2)
    tau3 <- pow(psi3, -2)
    
    for (i in 1:n) {
    y1[i] ~ dnorm(phi12, tau1)
    y2[i] ~ dnorm(phi23, tau3)
    }
  }"
  
  parameters <- c("phi12", "phi23",
                  "psi1", "psi3")
  
  foo <- jags.model(textConnection(model),
                    data = list(y1=data[,1], y2=data[,2], n=n),
                    quiet = T)
  update(foo, n.iter = burnin)
  jags <- coda.samples(model = foo, variable.names = parameters,
                       n.iter = Ntotal)
  out_bugs <- as.matrix(jags)
  
  return(out_bugs)
}


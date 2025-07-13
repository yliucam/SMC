fecundity_model <- function(data,
                            N,
                            Ntotal) {
  time_begin <- Sys.time()
  
  model <- "model{
    for (i in 1:ti) {
      rate[i] <- Nt[i] * rho
      nt[i] ~ dpois(rate[i])
    }
    
    rho ~ dunif(0, 10)
  }"
  
  parameters <- c("rho")
  
  foo <- jags.model(textConnection(model),
                    data = list(Nt=data[1:25, 1],
                                nt=data[1:25, 2],
                                ti=25),
                    quiet = T)
  
  update(foo)
  
  rho <- array(rep(NA, N*Ntotal), dim = c(N, Ntotal))
  
  for (i in 1:N) {
    jags_res <- coda.samples(model = foo, variable.names = parameters,
                             n.iter = Ntotal, thin = 1)
    out_jags <- as.matrix(jags_res)
    rho[i,] <- out_jags
  }
  
  time_end <- Sys.time()
  time_run <- time_end - time_begin
  
  return(list(rho=rho, tim=time_run))
  
}

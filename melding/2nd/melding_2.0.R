melding_2.0 <- function(data,
                        lambda,
                        Ntotal,
                        burnin) {
  data1 <- data$y1
  data3 <- data$y3
  data2 <- data$y2
    
  out_sub <- melding_jags(data = cbind(data1, data3),
                          Ntotal = Ntotal,
                          burnin = burnin)
  phi12 <- out_sub[,1]
  phi23 <- out_sub[,2]
  psi1 <- out_sub[,3]
  psi2 <- out_sub[,4]
  
  N <- length(phi12)
  
  out <- melding_2.0_mcmc(data = data2,
                          Ntotal = Ntotal,
                          burnin = burnin,
                          sub = cbind(phi12, phi23),
                          lambda = lambda,
                          N = N)
  phi12 <- out$phi12
  phi23 <- out$phi23
  psi2 <- out$psi2
  
  return(list(phi12=phi12, phi23=phi23, psi2=psi2))
}
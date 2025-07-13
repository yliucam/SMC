merging_loglik <- function(data,
                           alpha,
                           rho) {
  N <- dim(alpha)[1]
  
  alpha0 <- alpha[,1]
  alpha2 <- alpha[,2]
  ti <- length(data)
  
  deltaJF <- exp(alpha0) / (1 + exp(alpha0))
  deltaAF <- exp(alpha0+alpha2) / (1 + exp(alpha0+alpha2))
  
  eta_t <- 1
  
  y <- data
  
  rateJ <- array(rep(NA, N*(ti-1)), dim = c(N, ti-1))
  rate_imm <- array(rep(NA, N*(ti-1)), dim = c(N, ti-1))
  
  xJ <- array(rep(NA, N*ti), dim = c(N, ti))
  sur <- array(rep(NA, N*ti), dim = c(N, ti))
  imm <- array(rep(NA, N*ti), dim = c(N, ti))
  x <- array(rep(NA, N*ti), dim = c(N, ti))
  
  #log_lik <- rep(NA, N)
  
  xJ[,1] <- sample(1:50, N, replace = T)
  sur[,1] <- sample(1:25, N, replace = T)
  imm[,1] <- sample(1:25, N, replace = T)
  x[,1] <- sur[,1] + imm[,1] + xJ[,1]
  
  log_lik <- pois_log(data = y[1], rate = x[,1])
  
  for (tt in 2:ti) {
    rateJ[,tt-1] <- .5 * rho * deltaJF * x[,tt-1]
    xJ[,tt] <- sapply(rateJ[,tt-1], function(x) rpois(1, x))
    #xJ[,tt] <- rateJ[,tt-1]
    sur[,tt] <- sapply(x[,tt-1], function(x) rbinom(1, x, deltaAF))
    #sur[,tt] <- deltaAF * x[,tt-1]
    rate_imm[,tt-1] <- x[,tt-1] * eta_t
    imm[,tt] <- sapply(rate_imm[,tt-1], function(x) rpois(1, x))
    #imm[,tt] <- rate_imm[,tt-1]
    x[,tt] <- sur[,tt] + imm[,tt] + xJ[,tt]
    log_lik_inc <- pois_log(data = y[tt], rate = x[,tt])
    log_lik <- log_lik + log_lik_inc
  }
  
  return(log_lik)
  
}





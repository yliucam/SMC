#---------simple smc-sampler of some toy examples------------#



y <- rnorm(100, 1500, 10)




# Systematic resampling - Algorithm 9.6 in Chopin & Papaspiliopoulous (2020)

Sys_resamp <- function(W, P, U) {
  A <- rep(NA, P)
  v <- P * cumsum(W)
  s <- U
  m <- 1
  for (i in 1:P) {
    while (v[m] < s) {
      m <- m + 1
      #s <- s + 1
    }
    A[i] <- m
    s <- s + 1
  }
  
  return(A)
}


smc_sampler <- function(data, N, p=1) {
  if (length(data)%%p != 0) stop("p should be a divisor of the length of data")
  nt <- length(data)/p
  
  # Storage
  x <- array(rep(NA, N*nt), dim = c(N, nt))
  x_re <- array(rep(NA, N*nt), dim = c(N, nt))
  w_log <- array(rep(NA, N*nt), dim = c(N, nt))
  W <- array(rep(NA, N*nt), dim = c(N, nt))
  
  # Initialization
  x_0 <- rep(0, N)
  w_log_0 <- rep(0, N)
  L_log <- 0
  
  for (i in 1:(nt-1)) {
    if (i == 1) {
      x[,i] <- x_0 + rnorm(N, 0, 10000)
      for (j in 1:N) {
        w_log[j,i] <- w_log_0[j] + sum(sapply(data[1:p], function(y) dnorm(y, x[j,i], 10, log = T)))
        if (w_log[j,i] < log(.Machine$double.xmin)) w_log[j,i] <- log(.Machine$double.xmin)
      }
      W[,i] <- exp(w_log[,i]) / sum(exp(w_log[,i]))
      
      U <- runif(1, 0, 1)
      A <- Sys_resamp(W=W[,i], P=N, U=U)
      x_re[,i] <- x[A,i]
      
      L_log <- L_log + log(sum(exp(w_log[,i]))) - log(N)
    }
    
    x[,i+1] <- x_re[,i] + rnorm(N, 0, 10000)
    for (j in 1:N) {
      w_log[j,i+1] <- sum(sapply(data[(i*p+1):((i+1)*p)], function(y) dnorm(y, x[j,i+1], 10, log = T)))
      if (w_log[j,i+1] < log(.Machine$double.xmin)) w_log[j,i+1] <- log(.Machine$double.xmin)
    }
    W[,i+1] <- exp(w_log[,i+1]) / sum(exp(w_log[,i+1]))
    
    U <- runif(1, 0, 1)
    A <- Sys_resamp(W=W[,i+1], P=N, U=U)
    x_re[,i+1] <- x[A,i+1]
    
    L_log <- L_log + log(sum(exp(w_log[,i+1]))) - log(N)
    
    print(paste('iteration', i))
  }
  
  return(list(x=x_re[,nt], L_log=L_log))
}











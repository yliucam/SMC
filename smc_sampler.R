#---------simple smc-sampler of some toy examples------------#


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



#---------------Toy examples-------------------------------#


set.seed(240603)

##-------------Example 1 normal distribution-----------------------##

mu <- 150
sigma <- 40

data <- rnorm(100, mu, sigma)



smc_sampler <- function(data, N, p=1, Ntotal) {
  if (length(data)%%p != 0) stop("p should be a divisor of the length of data")
  nt <- length(data)/p
  data_bar <- mean(data)
  
  # Storage
  mu <- array(rep(NA, N*Ntotal*nt), dim = c(N, Ntotal, nt))
  mu_re <- array(rep(NA, N*nt), dim = c(N, nt))
  sigma <- array(rep(NA, N*Ntotal*nt), dim = c(N, Ntotal, nt))
  sigma_re <- array(rep(NA, N*nt), dim = c(N, nt))
  
  w_log <- array(rep(NA, N*nt), dim = c(N, nt))
  W <- array(rep(NA, N*nt), dim = c(N, nt))
  
  # Initialization
  mu_0 <- rep(0, N)
  sigma_0 <- rep(1/1000, N)
  w_log_0 <- rep(0, N)
  W_0 <- rep(1/N, N)
  #L_log <- 0
  
  for (i in 1:(nt-1)) {
    if (i == 1) {
      
      tim0 <- proc.time()
      
      U <- runif(1, 0, 1)
      A <- Sys_resamp(W=W_0, P=N, U=U)
      mu_0_re <- mu_0[A]
      sigma_0_re <- sigma_0[A]
      
      ## MCMC kernel
      for (j in 1:(Ntotal-1)) {
        if (j == 1) {
          mu_star <- rmvnorm(1, mu_0_re, 100*diag(N))
          sigma_star <- exp(rmvnorm(1, log(sigma_0_re), 100*diag(N)))
          
          for (jj in 1:N) {
            post_log <- dnorm(mu_0_re[jj], data_bar, 1, log = T) + dgamma(sigma_0_re[jj], 20, 1, log = T)
            + sum(sapply(data[1:p], function(y) dnorm(y, mu_0_re[jj], sigma_0_re[jj], log = T)))
            post_star_log <- dnorm(mu_star[jj], data_bar, 1, log = T) + dgamma(sigma_star[jj], 20, 1, log = T)
            + sum(sapply(data[1:p], function(y) dnorm(y, mu_star[jj], sigma_star[jj], log = T)))
            
            ### Accept/reject
            alpha <- min(0, post_star_log - post_log + log(sigma_star[jj]) - log(sigma_0_re[jj]))
            if(log(runif(1, 0, 1)) < alpha) {
              mu[jj,j,i] <- mu_star[jj]
              sigma[jj,j,i] <- sigma_star[jj]
            } else {
              mu[jj,j,i] <- mu_0_re[jj]
              sigma[jj,j,i] <- sigma_0_re[jj]
            }
          }
        }
        
        mu_star <- rmvnorm(1, mu[,j,i], 100*diag(N))
        sigma_star <- exp(rmvnorm(1, log(sigma[,j,i]), 100*diag(N)))
        
        for (jj in 1:N) {
          post_log <- dnorm(mu[jj,j,i], data_bar, 1, log = T) + dgamma(sigma[jj,j,i], 20, 1, log = T)
          + sum(sapply(data[1:p], function(y) dnorm(y, mu[jj,j,i], sigma[jj,j,i], log = T)))
          post_star_log <- dnorm(mu_star[jj], data_bar, 1, log = T) + dgamma(sigma_star[jj], 20, 1, log = T)
          + sum(sapply(data[1:p], function(y) dnorm(y, mu_star[jj], sigma_star[jj], log = T)))
          
          ### Accept/reject
          alpha <- min(0, post_star_log - post_log + log(sigma_star[jj]) - log(sigma[jj,j,i]))
          if(log(runif(1, 0, 1)) < alpha) {
            mu[jj,j+1,i] <- mu_star[jj]
            sigma[jj,j+1,i] <- sigma_star[jj]
          } else {
            mu[jj,j+1,i] <- mu[jj,j,i]
            sigma[jj,j+1,i] <- sigma[jj,j,i]
          }
        }
      }
      
      ## Update weights
      for (jj in 1:N) {
        w_log[jj,i] <- sum(sapply(data[1:p], function(y) dnorm(y, mu[jj,Ntotal,i], sigma[jj,Ntotal,i], log = T)))
        if (w_log[jj,i] < log(.Machine$double.xmin)) w_log[jj,i] <- log(.Machine$double.xmin)
      }
      if (max(w_log[,i]) == log(.Machine$double.xmin)) {
        W[,i] <- rep(1/N, N)
      } else {
        W[,i] <- exp(w_log[,i]) / sum(exp(w_log[,i]))  
      }
      
      
      #L_log <- L_log + log(sum(exp(w_log[,i]))) - log(N)
      
      tim1 <- proc.time()
      timT <- tim1-tim0
      print(paste0('iteration ', i, ', CPU time is ', round((tim1-tim0)[[3]], 2), 's, ', 'total CPU time is ', round(timT[[3]], digits=2), 's'))
    }
    
    tim0 <- proc.time()
    
    U <- runif(1, 0, 1)
    A <- Sys_resamp(W=W[,i], P=N, U=U)
    mu_re[,i] <- mu[A,Ntotal,i]
    sigma_re[,i] <- sigma[A,Ntotal,i]
    
    
    ## MCMC kernel
    for (j in 1:(Ntotal-1)) {
      if (j == 1) {
        mu_star <- rmvnorm(1, mu_re[,i], 1000*diag(N))
        sigma_star <- exp(rmvnorm(1, log(sigma_re[,i]), 1000*diag(N)))
        
        for (jj in 1:N) {
          post_log <- dnorm(mu_re[jj,i], data_bar, 1, log = T) + dgamma(sigma_re[jj,i], 20, 1, log = T)
          + sum(sapply(data[(i*p+1):((i+1)*p)], function(y) dnorm(y, mu_re[jj,i], sigma_re[jj,i], log = T)))
          post_star_log <- dnorm(mu_star[jj], data_bar, 1, log = T) + dgamma(sigma_star[jj], 20, 1, log = T)
          + sum(sapply(data[(i*p+1):((i+1)*p)], function(y) dnorm(y, mu_star[jj], sigma_star[jj], log = T)))
          
          ### Accept/reject
          alpha <- min(0, post_star_log - post_log + log(sigma_star[jj]) - log(sigma_re[jj,i]))
          if(log(runif(1, 0, 1)) < alpha) {
            mu[jj,j,i+1] <- mu_star[jj]
            sigma[jj,j,i+1] <- sigma_star[jj]
          } else {
            mu[jj,j,i+1] <- mu_re[jj,i]
            sigma[jj,j,i+1] <- sigma_re[jj,i]
          }
        }
      }
      
      mu_star <- rmvnorm(1, mu[,j,i+1], 1000*diag(N))
      sigma_star <- exp(rmvnorm(1, log(sigma[,j,i+1]), 1000*diag(N)))
      
      for (jj in 1:N) {
        post_log <- dnorm(mu[jj,j,i+1], data_bar, 1, log = T) + dgamma(sigma[jj,j,i+1], 20, 1, log = T)
        + sum(sapply(data[(i*p+1):((i+1)*p)], function(y) dnorm(y, mu[jj,j,i+1], sigma[jj,j,i+1], log = T)))
        post_star_log <- dnorm(mu_star[jj], data_bar, 1, log = T) + dgamma(sigma_star[jj], 20, 1, log = T)
        + sum(sapply(data[(i*p+1):((i+1)*p)], function(y) dnorm(y, mu_star[jj], sigma_star[jj], log = T)))
        
        ### Accept/reject
        alpha <- min(0, post_star_log - post_log + log(sigma_star[jj]) - log(sigma[jj,j,i+1]))
        if(log(runif(1, 0, 1)) < alpha) {
          mu[jj,j+1,i+1] <- mu_star[jj]
          sigma[jj,j+1,i+1] <- sigma_star[jj]
        } else {
          mu[jj,j+1,i+1] <- mu[jj,j,i+1]
          sigma[jj,j+1,i+1] <- sigma[jj,j,i+1]
        }
      }
    }
    
    
    ## Update weigths
    for (jj in 1:N) {
      w_log[jj,i+1] <- sum(sapply(data[(i*p+1):((i+1)*p)], function(y) dnorm(y, mu[jj,Ntotal,i+1], sigma[jj,Ntotal,i+1], log = T)))
      if (w_log[jj,i+1] < log(.Machine$double.xmin)) w_log[jj,i+1] <- log(.Machine$double.xmin)
    }
    W[,i+1] <- exp(w_log[,i+1]) / sum(exp(w_log[,i+1]))
    
    #L_log <- L_log + log(sum(exp(w_log[,i+1]))) - log(N)
    
    tim1 <- proc.time()
    timT <- timT + tim1 - tim0
    print(paste0('iteration ', i+1, ', CPU time is ', round((tim1-tim0)[[3]], 2), 's, ', 'total CPU time is ', round(timT[[3]], digits=2), 's'))
  }
  
  U <- runif(1, 0, 1)
  A <- Sys_resamp(W=W[,nt], P=N, U=U)
  
  return(list(mu_post=mu, sigma_post=sigma, A=A))
}



##-------------Example 2 mixture normal distribution-----------------##

omega1 <- .1
omega2 <- .6
omega3 <- .3
mu1 <- 100
mu2 <- 1
mu3 <- .1
sigma1 <- 1
sigma2 <- .1
sigma3 <- 100


data <- omega1*rnorm(100, mu1, sigma1) + omega2*rnorm(100, mu2, sigma2) + omega3*rnorm(100, mu3, sigma3)







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












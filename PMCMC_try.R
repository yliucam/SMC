u <- rnorm(100, 0, 1)
v <- rnorm(100, 0, .2)

x <- rep(NA, 100)
x[1] <- rnorm(1, 0, 1)
for (i in 2:100) {
  x[i] <- x[i-1] * .9 + u[i]
}

y <- rep(NA, 100)
for (i in 1:100) {
  y[i] <- x[i] + v[i]
}



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




gibbs_pmmh <- function(data, P, Ntotal, burnin, thin) {
  TT <- length(data)
  ESS_min <- P / 2
  
  # Storage
  x <- array(rep(NA, Ntotal*TT*P), dim = c(TT, P, Ntotal))
  rho <- rep(NA, Ntotal)
  sigma_x <- rep(NA, Ntotal)
  sigma_y <- rep(NA, Ntotal)
  w <- matrix(rep(NA, Ntotal*P), ncol = Ntotal)
  v <- rep(NA, P)
  
  # Initialisation
  x[,,1] <- 0
  rho[1] <- 0
  sigma_x[1] <- 1
  sigma_x[1] <- 1
  w[,1] <- rep(1, P)
  W <- rep(1/P, P)
  
  for (i in 1:Ntotal) {
    ## Normalisation constant initialisation
    if (i == 1) {
      # Bootstrap PF - Algorithm 10.4 in Chopin & Papaspiliopoulous (2020)
      L_current <- 1
      for (j in 1:TT) {
        ESS <- (sum(w[,i]))^2 / sum((w[,i])^2)
        if (ESS < ESS_min) {
          ## Resampling
          U <- runif(1, 0, 1)
          A <- Sys_resamp(W, P, U)
          w_star <- rep(1, P)  
        } else {
          A <- 1:P
          w_star <- w[,j]
        }
        w_sum <- sum(w[,j+1])
        W[,j+1] <- W[,j+1] / w_sum
        L_current <- L_current * mean(w[,j+1])
      } 
    }
    
    ## Proposals
    rho_star <- rho[i] + rnorm(0, 1)
    sigma_x_star <- sigma_x[i] + rnorm(0, 1)
    sigma_y_star <- sigma_y[i] + rnorm(0, 1)
    
    x_star <- matrix(rep(NA, TT*P), ncol = P)
    
    L_prop <- 1
    
    # Bootstrap PF
    for (j in 1:TT) {
      ESS <- (sum(w[,i]))^2 / sum((w[,i])^2)
      if (ESS < ESS_min) {
        ## Resampling
        U <- runif(1, 0, 1)
        A <- Sys_resamp(W, P, U)
        w_star <- rep(1, P)
      } else {
        A <- 1:P
        w_star <- w[,j]
      }
      x_star[j,] <- apply(x[j,,i][A], 1, function(x) x*rho_star + rnorm(1, 0, sigma_x_star))
      for (jj in 1:P) {
        w[jj,j+1] <- w_star[jj] * dnorm(x_star[j,jj], sigma_y_star)
      }
      w_sum <- sum(w[,j+1])
      W[,j+1] <- W[,j+1] / w_sum
      L_prop <- L_prop * mean(w[,j+1])
    }
    
    ## Posteriors
    p_num_log <- dnorm(rho_star,0,1,log=T) + dnorm(sigma_x_star,0,1,log=T) + dnorm(sigma_y_star,0,1,log=T)
              + log(L_prop) + dnorm(rho_[i],rho_star,1,log=T) + dnorm(sigma_x[i], sigma_x_star,1,log=T)
              + dnorm(sigma_y[i],sigma_y_star,1,log=T)
    p_den_log <- dnorm(rho[i],0,1,log=T) + dnorm(sigma_x[i],0,1,log=T) + dnorm(sigma_y[i],0,1,log=T)
              + log(L_current) + dnorm(rho_star,rho_[i],1,log=T) + dnorm(sigma_x_star,sigma_x[i],1,log=T)
              + dnorm(sigma_y_star,sigma_y[i],1,log=T)
    ## Accept/reject
    alpha <- min(0, p_num_log - p_den_log)
    if (log(runif(1, 0, 1)) < alpha) {
      x
    }
  }
  
  
}








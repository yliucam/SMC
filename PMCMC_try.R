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




gibbs_pmmc <- function(data, P, Ntotal, burnin, thin) {
  TT <- length(data)
  ESS_min <- P / 2
  
  # Storage
  x <- array(rep(NA, Ntotal*TT*P), dim = c(P, TT, Ntotal))
  rho <- rep(NA, Ntotal)
  sigma_x <- rep(NA, Ntotal)
  sigma_y <- rep(NA, Ntotal)
  w <- matrix(rep(NA, TT*P), ncol = TT)
  v <- rep(NA, P)
  
  # Initialisation
  x_0 <- rep(0, P)
  rho_0 <- 0
  sigma_x_0 <- 1
  sigma_y_0 <- 1
  w_0 <- rep(1, P)
  W_0 <- rep(1/P, P)
  
  W_bug <- matrix(rep(0, P*TT), ncol=TT)
  A_bug <- matrix(rep(NA, P*TT), ncol=TT)
  
  for (i in 1:Ntotal) {
    ## Normalisation constant initialisation
    if (i == 1) {
      # Bootstrap PF - Algorithm 10.4 in Chopin & Papaspiliopoulous (2020)
      L_current <- 1
        for (j in 1:TT) {
          if (j == 1) {
            x[,j,i] <- x_0 * rho_0 + rnorm(P, 0, sigma_x_0)
            for (jj in 1:P) {
              w[jj,j] <- w_0[jj] * dnorm(x[jj,j,i], sigma_y_0)
            }
            W <- w[,j] / sum(w[,j])
            ## The likelihood is updated following (10.3) in Chopin & Papaspiliopoulous (2020)
            ## No resampling occurs here
            L_current <- L_current * sum(w[,j]) / sum(w_0)
            W_bug[,j] <- W
          } else {
            ESS <- (sum(w[,j-1]))^2 / sum((w[,j-1])^2)  
            if (ESS < ESS_min) {
              ## Resampling
              U <- runif(1, 0, 1)
              A <- Sys_resamp(W=W, P=P, U=U)
              w_star <- rep(1, P)  
            } else {
              A <- 1:P
              w_star <- w[,j-1]
            }
            x[,j,i] <- x[,j-1,i][A] * rho_0 + rnorm(P, 0, sigma_x_0)
            for (jj in 1:P) {
              w[jj,j] <- w_star[jj] * dnorm(x[jj,j,i], sigma_y_0) ## (10.4b) in Chopin & Papaspiliopoulous (2020)
            }
            W <- w[,j] / sum(w[,j])
            W_bug[,j] <- W
            ## The likelihood is updated following (10.3) in Chopin & Papaspiliopoulous (2020)
            ## The update depends on whether resampling occurs at time j
            if (ESS < ESS_min) {
              L_current <- L_current * mean(w[,j]) 
            } else {
              L_current <- L_current * sum(w[,j]) / sum(w_star)
            }
          }
          
        }
      
        rho[i] <- rho_0
        sigma_x[i] <- sigma_x_0
        sigma_y[i] <- sigma_y_0
      
    } else {
      ## Proposals 
      rho_star <- rho[i-1] + rnorm(1, 0, 1)
      sigma_x_star <- exp(rnorm(1, log(sigma_x[i-1]), 1)) ## Log-normal proposal for scales
      sigma_y_star <- exp(rnorm(1, log(sigma_y[i-1]), 1)) ##
        
      x_star <- matrix(rep(NA, TT*P), ncol = TT)
        
      L_prop <- 1
      
      ESS_bug <- rep(NA, TT)
      
      # Bootstrap PF
      for (j in 1:TT) {
        if (j == 1) {
          ESS <- (sum(w[,j]))^2 / sum((w[,j])^2)
          ESS_bug[j] <- ESS
          if (ESS < ESS_min) {
            ## Resampling
            U <- runif(1, 0, 1)
            A <- Sys_resamp(W, P, U)
            w_star <- rep(1, P)
          } else {
            A <- 1:P
            w_star <- w[,j]
          }
          x_star[,j] <- x[A,j,i-1] * rho_star + rnorm(P, 0, sigma_x_star)
          for (jj in 1:P) {
            w[jj,j] <- w_star[jj] * dnorm(x_star[jj,j], sigma_y_star)
          }
          W <- w[,j] /sum(w[,j])
          if (ESS < ESS_min) {
            L_prop <- L_prop * mean(w[,j])
          } else {
            L_prop <- L_prop * sum(w[,j]) / sum(w_star)  
          }
        } else {
          ESS <- (sum(w[,j-1]))^2 / sum((w[,j-1])^2)
          ESS_bug[j] <- ESS
          if (ESS < ESS_min) {
            ## Resampling
            U <- runif(1, 0, 1)
            A <- Sys_resamp(W, P, U)
            w_star <- rep(1, P)
          } else {
            A <- 1:P
            w_star <- w[,j-1]
          }
          x_star[,j] <- x[,j-1,i][A] * rho_star + rnorm(P, 0, sigma_x_star)
          for (jj in 1:P) {
            w[jj,j] <- w_star[jj] * dnorm(x_star[jj,j], sigma_y_star)
          }
          W <- w[,j] / sum(w[,j])
          if (ESS < ESS_min) {
            L_prop <- L_prop * mean(w[,j])
          } else {
            L_prop <- L_prop * sum(w[,j]) / sum(w_star)
          }
        }
      }
      
      # MH sampling
      ## Posteriors
      p_num_log <- log(L_prop) + dnorm(rho_[i],rho_star,1,log=T) + dnorm(sigma_x[i], sigma_x_star,1,log=T)
      + dnorm(sigma_y[i],sigma_y_star,1,log=T)
      p_den_log <- dlog(L_current) + dnorm(rho_star,rho_[i],1,log=T) + dnorm(sigma_x_star,sigma_x[i],1,log=T)
      + dnorm(sigma_y_star,sigma_y[i],1,log=T)
      
      ## Jacobians for sigma_x and sigma_y
      Jacob_x_log <- -log(sigma_x[i]) - (-log(sigma_x_star))
      Jacob_y_log <- -log(sigma_y[i]) - (-log(sigma_y_star))
      
      ## Accept/reject
      alpha <- min(0, p_num_log - p_den_log + Jacob_x_log + Jacob_y_log)
      if (log(runif(1, 0, 1)) < alpha) {
        x[,,i+1] <- x_star
        rho[i+1] <- rho_star
        sigma_x[i+1] <- sigma_x_star
        sigma_y[i+1] <- sigma_y_star
        L_current <- L_prop
      } else {
        x[,,i+1] <- x[,,i]
        rho[i+1] <- rho[i]
        sigma_x[i+1] <- sigma_x[i]
        sigma_y[i+1] <- sigma_y[i]
      }
    }
  }
    
    
  keep <- seq((burnin+1), Ntotal, by = thin)
  
  return(list(x=x[,,keep],
              rho=rho[keep],
              sigma_x=sigma_x[keep],
              sigma_y=sigma_y[keep]))
}

















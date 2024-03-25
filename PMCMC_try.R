# To do -- Try a (much) smaller T


x <- rep(NA, 100)
x[1] <- rnorm(1, 0, 1)
for (i in 2:100) {
  x[i] <- rnorm(1, x[i-1]*.9, 1)
}

y <- rep(NA, 100)
for (i in 1:100) {
  y[i] <- rnorm(1, x[i], .2)
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




gibbs_pmmc <- function(data, P, Ntotal, burnin, thin, print_interval=100) {
  TT <- length(data)
  ESS_min <- P / 2
  
  # Storage
  x <- array(rep(NA, Ntotal*TT*P), dim = c(P, TT, Ntotal))
  rho <- rep(NA, Ntotal)
  sigma_x <- rep(NA, Ntotal)
  sigma_y <- rep(NA, Ntotal)
  w <- matrix(rep(NA, TT*P), ncol = TT)
  
  # Initialisation
  x_0 <- rep(0, P)
  rho_0 <- 0
  sigma_x_0 <- 1
  sigma_y_0 <- 1
  w_0 <- rep(1, P)
  
  
  for (i in 1:(Ntotal-1)) {
    ## Normalisation constant initialisation
    if (i == 1) {
      # Bootstrap PF - Algorithm 10.4 in Chopin & Papaspiliopoulous (2020)
      L_current <- 1
      for (j in 1:TT) {
        if (j == 1) {
          x[,j,i] <- x_0 * rho_0 + rnorm(P, 0, sigma_x_0)
          for (jj in 1:P) {
            w[jj,j] <- w_0[jj] * dnorm(data[j], x[jj,j,i], sigma_y_0) * 1e2
          }
          ## Ensure W is not NaN -- dangerous, need to change
          if (sum(w[,j]) == 0) {
            W <- rep(1/P, P)
          } else {
            W <- w[,j] / sum(w[,j])  
          }
          ## The likelihood is updated following (10.3) in Chopin & Papaspiliopoulous (2020)
          ## No resampling occurs here
          L_current <- L_current * sum(w[,j]) / sum(w_0)
        } else {
          ## Ensure ESS is not NaN -- dangerous, need to change
          if (((sum(w[,j-1]))^2==0)&&(sum((w[,j-1])^2)==0)) {
            ESS <- 0
          } else {
            ESS <- (sum(w[,j-1]))^2 / sum((w[,j-1])^2)  
          }
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
            w[jj,j] <- w_star[jj] * dnorm(data[j], x[jj,j,i], sigma_y_0) * 1e2 ## (10.4b) in Chopin & Papaspiliopoulous (2020)
          }
          ## Ensure W is not NaN -- dangerous, need to change
          if (sum(w[,j]) == 0) {
            W <- rep(1/P, P)
          } else {
            W <- w[,j] / sum(w[,j])  
          }
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
      
    }
      
    ## Proposals 
    rho_star <- rho[i] + rnorm(1, 0, .05)
    sigma_x_star <- rlnorm(1, log(sigma_x[i]), .05) ## Log-normal proposal for scales
    sigma_y_star <- rlnorm(1, log(sigma_y[i]), .05) ##
      
    x_star <- matrix(rep(NA, TT*P), ncol = TT)
      
    L_prop <- 1
      
      
    # Bootstrap PF
    for (j in 1:(TT-1)) {
      if (j == 1) {
        ## Ensure ESS is not NaN -- dangerous, need to change
        if (((sum(w[,j]))^2==0)&&(sum((w[,j])^2)==0)) {
          ESS <- 0
        } else {
          ESS <- (sum(w[,j]))^2 / sum((w[,j])^2)  
        }
        if (ESS < ESS_min) {
          ## Resampling
          U <- runif(1, 0, 1)
          A <- Sys_resamp(W, P, U)
          w_star <- rep(1, P)
        } else {
          A <- 1:P
          w_star <- w[,j]
        }
        x_star[,j] <- x[A,j,i] * rho_star + rnorm(P, 0, sigma_x_star)
        for (jj in 1:P) {
          w[jj,j] <- w_star[jj] * dnorm(data[j], x_star[jj,j], sigma_y_star) * 1e2
        }
        ## Ensure W is not NaN -- dangerous, need to change
        if (sum(w[,j+1]) == 0) {
          W <- rep(1/P, P)
        } else {
          W <- w[,j+1] / sum(w[,j+1])  
        }
        if (ESS < ESS_min) {
          L_prop <- L_prop * mean(w[,j])
        } else {
          L_prop <- L_prop * sum(w[,j]) / sum(w_star)  
        }
      }
      ## Ensure ESS is not NaN -- dangerous, need to change
      if (((sum(w[,j]))^2==0)&&(sum((w[,j])^2)==0)) {
        ESS <- 0
      } else {
        ESS <- (sum(w[,j]))^2 / sum((w[,j])^2)  
      }
      if (ESS < ESS_min) {
        ## Resampling
        U <- runif(1, 0, 1)
        A <- Sys_resamp(W=W, P=P, U=U)
        w_star <- rep(1, P)
      } else {
        A <- 1:P
        w_star <- w[,j]
      }
      x_star[,j+1] <- x_star[A,j] * rho_star + rnorm(P, 0, sigma_x_star)
      for (jj in 1:P) {
        w[jj,j+1] <- w_star[jj] * dnorm(data[j], x_star[jj,j+1], sigma_y_star) * 1e2
      }
      ## Ensure W is not NaN -- dangerous, need to change
      if (sum(w[,j+1]) == 0) {
        W <- rep(1/P, P)
      } else {
        W <- w[,j+1] / sum(w[,j+1])  
      }
      
      if (ESS < ESS_min) {
        L_prop <- L_prop * mean(w[,j+1])
      } else {
        L_prop <- L_prop * sum(w[,j+1]) / sum(w_star)
      }
  
    }
      
    # MH sampling
    ## Posteriors
    p_num_log <- log(L_prop) + dunif(rho_star, -1, 1, log = T) 
                    + dgamma(sigma_x_star^2, .5, .5, log = T)
                      + dgamma(sigma_y_star^2, .5, .5, log = T)
    p_den_log <- log(L_current) + dunif(rho[i], -1, 1, log = T) 
                    + dgamma(sigma_x[i]^2, .5, .5, log = T)
                      + dgamma(sigma_y[i]^2, .5, .5, log = T)
      
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
    
    if (i%%print_interval == 0) print(paste0('MCMC iterations ', i, '/', Ntotal))
  }
  
  
  keep <- seq((burnin+1), Ntotal, by = thin)
  
  return(list(x=x[,,keep],
              rho=rho[keep],
              sigma_x=sigma_x[keep],
              sigma_y=sigma_y[keep]))
}



out_pmcmc <- gibbs_pmmc(data = y, P = 100, Ntotal = 10000, burnin = 5000, thin = 2, print_interval = 1000)


hist(out_pmcmc$rho, breaks = 30)
hist(out_pmcmc$sigma_x, breaks = 30)
hist(out_pmcmc$sigma_y, breaks = 30)

library(coda)

traceplot(as.mcmc(out_pmcmc$rho))
traceplot(as.mcmc(out_pmcmc$sigma_x))
traceplot(as.mcmc(out_pmcmc$sigma_y))


mean(out_pmcmc$rho)
median(out_pmcmc$rho)

mean(out_pmcmc$sigma_x)
median(out_pmcmc$sigma_x)

mean(out_pmcmc$sigma_y)
median(out_pmcmc$sigma_y)









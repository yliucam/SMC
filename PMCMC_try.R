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




gibbs_pmmh <- function(data, Ntotal, burnin, thin) {
  TT <- length(data)
  
  # Storage
  x <- matrix(rep(NA, Ntotal*TT), ncol = Ntotal)
  rho <- rep(NA, Ntotal)
  sigma_x <- rep(NA, Ntotal)
  sigma_y <- rep(NA, Ntotal)
  w <- matrix(rep(NA, Ntotal*TT), ncol = Ntotal)
  
  # Initialisation
  x[,1] <- matrix(rep(0, Ntotal*TT), ncol = Ntotal)
  rho[1] <- 0
  sigma_x[1] <- 1
  sigma_x[1] <- 1
  w_0 <- matrix(rep(1, Ntotal*TT), ncol = Ntotal)
  W_0 <- matrix(rep(1/Ntotal, Ntotal*TT), ncol = Ntotal)
  L_0 <- rowMeans(w_0)
  
  for (i in 1:Ntotal) {
    rho_star <- rnorm(1, rho[i])
    sigma_x_star <- rnorm(1, sigma_x[i])
    sigma_y_star <- rnorm(1, sigma_y[i])
    
    for (j in 1:TT) {
      if (j == 1) {
        
      }
    }
  }
  
  
}

theta <- rep(NA, 11)

theta[1] <- rexp(1)
for (i in 2:11) {
  theta[i] <- rexp(1, 1/theta[i-1])
}

x <- rep(NA, 10)
y <- rep(NA, 10)
for (i in 1:10) {
  x[i] <- rnorm(1, 0, theta[i+1])
  y[i] <- rnorm(1, x[i], 1)
}



dc_smc_A1 <- function(data,
                      N,
                      alpha,
                      Ntotal) {
  n <- length(data)
  
  # Storage
  theta <- array(rep(NA, N*n), dim = c(N, n))
  w_log <- array(rep(NA, N*n), dim = c(N, n))
  W <- array(rep(NA, N*n), dim = c(N, n)) 
  
  ## Leaf storage
  x_leaf <- array(rep(NA, N*n), dim = c(N, n))
  w_leaf_log <- array(rep(NA, N*n), dim = c(N, n))
  W_leaf <- array(rep(NA, N*n), dim = c(N, n))
  
  
  #Initialization
  theta_0 <- rexp(1, 1)
  L_prod_log <- 0
  
  for (t in 0:(n-1)) {
    
    # Leaf node
    if (t == 0) {
      for (t_leaf in 1:n) {
        
      }
    }
  }
  
}































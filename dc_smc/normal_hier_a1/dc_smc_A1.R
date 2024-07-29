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



dc_smc_A1 <- function(data,
                      N,
                      alpha,
                      m,
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
  
  
  ## Meriging storage
  x_merge <- rep(NA, m*N)
  W_x_merge <- rep(NA, m*N)
  theta_merge <- rep(NA, m*N)
  W_theta_merge <- rep(NA, m*N)
  v_t_log <- rep(NA, m*N)
  p_check_log <- rep(NA, m*N)
  x_merge_resamp <- rep(NA, N)
  p_check_log_resamp <- rep(NA, N)
  
  
  # Initialization
  theta_0 <- rexp(1, 1)
  L_prod_log <- 0
  
  for (t in 1:n) {
    
    # Leaf node
    if (t == 1) {
      out_leaf <- dc_leaf(y = data, n = n, prior_mu = prior_mu, prior_sigma = prior_sigma, N = N, min_lim = min_lim)
      x_leaf <- out_leaf$x
      w_leaf_log <- out_leaf$w_log
      W_leaf <- out_leaf$W
      theta[,t] <- out_leaf$theta_0
      w_log[,t] <- out_leaf$w_theta0_log
      W[,t] <- out_leaf$W_theta0
    }
    
    
    
    # No merging steps as p_check() is choosen to be prod(p_c())
    
    # SMC sampler
    ## Initialization
    theta_tilde <- sapply(theta[,t], function(x) rexp(1, x))
    w_log_0 <- 0
    
    for(i in 1:(nt-1)) {
      if (i == 1) {
        ## MCMC kernel
      }
    }
    
    
    
  }
  
}

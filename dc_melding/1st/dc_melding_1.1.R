# Used for psi3 being Gamma distributed




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




dc_melding_1.1 <- function(data,
                           N,
                           n_leaf,
                           m,
                           alpha,
                           Ntotal) {
  ESS_min <- N / 2
  min_lim <- log(.Machine$double.xmin)
  
  n <- dim(data)[1]
  
  nt <- 1 / alpha
  
  # Storage
  ## Leaf storage
  x_leaf <- array(rep(NA, N*n_leaf), dim = c(N, n_leaf))
  w_leaf_log <- array(rep(NA, N*n_leaf), dim = c(N, n_leaf))
  W_leaf <- array(rep(NA, N*n_leaf), dim = c(N, n_leaf))
  L_leaf_log <- rep(NA, n_leaf)
  ### New leaf storage after merging
  x_leaf_merge <- array(rep(NA, N*n_leaf), dim = c(N, n_leaf))
  
  
  ## Root storage
  x_root <- rep(NA, N)
  x_root_array <- array(rep(NA, N*nt), dim = c(N, nt))
  w_root_log <- array(rep(NA, N*nt), dim = c(N, nt))
  W_root <- rep(NA, N)
  W_root_array <- array(rep(NA, N*nt), dim = c(N, nt))
  
  for (leaf_i in 1:n_leaf) {
    out_leaf <- dc_melding_leaf_mcmc(data = data[,leaf_i],
                                     N = N,
                                     Ntotal = Ntotal)
    x_leaf[,leaf_i] <- out_leaf$mu
    w_leaf_log[,leaf_i] <- out_leaf$w_log
    W_leaf[,leaf_i] <- out_leaf$W
  }
  
  
  # Merging
  ## Storage for mN matchings
  x_leaf_mN <- array(rep(NA, m*N*n_leaf), dim = c(m*N, n_leaf))
  w_leaf_log_mN <- array(rep(NA, m*N*n_leaf), dim = c(m*N, n_leaf))
  W_leaf_mN <- array(rep(NA, m*N*n_leaf), dim = c(m*N, n_leaf))
  
  ## Resampling mN particles -- using systematic resampling may be incorrect, try using generating indices from Unif()!!!
  for (leaf_i in 1:n_leaf) {
    U <- runif(1, 0, 1)
    A <- Sys_resamp(W=W_leaf[,leaf_i], P=m*N, U=U)
    x_leaf_mN[,leaf_i] <- x_leaf[A, leaf_i]
    w_leaf_log_mN[,leaf_i] <- w_leaf_log[A, leaf_i]
    W_leaf_mN[,leaf_i] <- W_leaf[A, leaf_i]
  }
  
  
  ## Computing weights for all mN matchings
  x_leaf_mN_means <- colMeans(x_leaf_mN)
  u_c_log_mN <- apply(x_leaf_mN, 1, function(x) dnorm(x[1], 0, 1, log = T) + dnorm(x[2], 0, 1, log = T))
  p_c_marginal_log_mN <- apply(x_leaf_mN, 1, function(x) dnorm(x[1], x_leaf_mN_means[1], 1, log = T) + dnorm(x[2], x_leaf_mN_means[2], 2, log = T))
  V_log <- alpha * (p_c_marginal_log_mN - u_c_log_mN)
  V_log[which(V_log < min_lim)] <- min_lim
  V <- exp(V_log)
  V <- V / sum(V)
  
  ## Resampling N matchings from mN ones using V
  U <- runif(1, 0, 1)
  A <- Sys_resamp(W=V, P=N, U=U)
  x_leaf_merge <- x_leaf_mN[A,]
  p_c_marginal_log <- p_c_marginal_log_mN[A]
  
  # Root
  #x_root_0 <- apply(x_leaf_merge, 1, function(x) rnorm(1, (4*x[1]+x[2])/9, 2/3))
  x_root_0 <- rgamma(N, 1, 1)
  w_root_log_0 <- rep(0, N)
  
  A_root <- array(rep(NA, N*nt), dim = c(N, nt))
  
  ## SMC sampler
  for (i in 1:(nt-1)) {
    if (i == 1) {
      ### Update weights
      p_root_log <- apply(cbind(x_leaf_merge, x_root_0), 1, function(x) dnorm(x[1], x[3], 1, log = T) + dnorm(x[2], x[3], 2, log = T))
      u_root_log <- sapply(x_root_0, function(x) dgamma(x, 1, 1, log = T))
      u_c_log <- apply(x_leaf_merge, 1, function(x) dnorm(x[1], 0, 1, log = T) + dnorm(x[2], 0, 1, log = T))
      q_root_log <- u_root_log ## Set q_t to be the same as u_t
      w_root_log[,i] <- w_root_log_0 + alpha * (p_root_log + u_root_log - (1-alpha) * u_c_log - alpha * p_c_marginal_log - q_root_log)
      w_root_log[which(w_root_log[,i] < min_lim), i] <- min_lim
      W_root <- exp(w_root_log[,i] - matrixStats::logSumExp(w_root_log[,i]))
      W_root_array[,i] <- W_root
      
      ### Resampling - optionally
      ESS <- 1 / sum(W_root^2)
      if (ESS < ESS_min) {
        U <- runif(1, 0, 1)
        A <- Sys_resamp(W=W_root, P=N, U=U)
        x_root_resamp <- x_root_0[A]
        x_leaf_merge <- x_leaf_merge[A,]
        u_c_log <- u_c_log[A]
        p_c_marginal_log <- p_c_marginal_log[A]
        w_root_log[,i] <- rep(0, N)
        A_root[,i] <- A
      } else {
        x_root_resamp <- x_root_0
      }
      
      ### Update alpha
      alpha_update <- alpha
      
      ### MCMC kernel
      out_root <- dc_melding_root_mcmc(data=x_leaf_merge, Ntotal=Ntotal, mu_init=x_root_resamp, alpha_j=alpha_update)
      x_root <- out_root$mu[,Ntotal]
      x_root_array[,i] <- x_root
    }
    
    ### Update weights
    p_root_log <- apply(cbind(x_leaf_merge, x_root), 1, function(x) dnorm(x[1], x[3], 1, log = T) + dnorm(x[2], x[3], 2, log = T))
    u_root_log <- sapply(x_root, function(x) dgamma(x, 1, 1, log = T))
    u_c_log <- apply(x_leaf_merge, 1, function(x) dnorm(x[1], 0, 1, log = T) + dnorm(x[2], 0, 1, log = T))
    q_root_log <- u_root_log ## Set q_t to be the same as u_t
    w_root_log[,i+1] <- w_root_log[,i] + alpha * (p_root_log + u_root_log - (1-alpha) * u_c_log - alpha * p_c_marginal_log - q_root_log)
    w_root_log[which(w_root_log[,i+1] < min_lim), i+1] <- min_lim
    W_root <- exp(w_root_log[,i+1] - matrixStats::logSumExp(w_root_log[,i+1]))
    W_root_array[,i+1] <- W_root
    
    ### Resampling - optionally
    ESS <- 1 / sum(W_root^2)
    if (ESS < ESS_min) {
      U <- runif(1, 0, 1)
      A <- Sys_resamp(W=W_root, P=N, U=U)
      x_root_resamp <- x_root[A]
      x_leaf_merge <- x_leaf_merge[A,]
      u_c_log <- u_c_log[A]
      p_c_marginal_log <- p_c_marginal_log[A]
      w_root_log[,i+1] <- rep(0, N)
      A_root[,i+1] <- A
    } else {
      x_root_resamp <- x_root
    }
    
    ### Update alpha
    alpha_update <- alpha_update + alpha
    
    ### MCMC kernel
    out_root <- dc_melding_root_mcmc(data=x_leaf_merge, Ntotal=Ntotal, mu_init=x_root_resamp, alpha_j=alpha_update)
    x_root <- out_root$mu[,Ntotal]
    x_root_array[,i+1] <- x_root
  }
  
  return(list(x=x_root_array, W=W_root_array, x_leaf=x_leaf, x_leaf_merge=x_leaf_merge))
}

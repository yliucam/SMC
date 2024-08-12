theta <- rep(NA, 11)

theta[1] <- rexp(1)
for (i in 2:11) {
  theta[i] <- rexp(1, 1/theta[i-1])
}

x <- rep(NA, 11)
y <- rep(NA, 11)
for (i in 1:11) {
  x[i] <- rnorm(1, 0, theta[i])
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



library(matrixStats)


dc_smc_A1 <- function(data,
                      N,
                      alpha,
                      Ntotal,
                      prior_mu,
                      prior_sigma) {
  n <- length(data)
  
  nt <- 1 / alpha
  
  # Storage
  theta <- array(rep(NA, N*n), dim = c(N, n))
  w_log <- array(rep(NA, N*n), dim = c(N, n))
  W <- array(rep(NA, N*n), dim = c(N, n)) 
  
  ## Leaf storage
  x_leaf <- array(rep(NA, N*n), dim = c(N, n))
  w_leaf_log <- array(rep(NA, N*n), dim = c(N, n))
  W_leaf <- array(rep(NA, N*n), dim = c(N, n))
  theta_0 <- rep(NA, N)
  w_theta0_log <- rep(NA, N)
  W_theta0 <- rep(NA, N)
  
  
  ## Meriging storage
  #x_merge <- rep(NA, m*N)
  #W_x_merge <- rep(NA, m*N)
  #theta_merge <- rep(NA, m*N)
  #W_theta_merge <- rep(NA, m*N)
  #v_t_log <- rep(NA, m*N)
  #p_check_log <- rep(NA, m*N)
  #x_merge_resamp <- rep(NA, N)
  #p_check_log_resamp <- rep(NA, N)
  
  
  # Initialization
  theta_0 <- rexp(1, 1)
  L_prod_log <- 0
  theta_density_forward_log <- 0 # f(theta_t|theta_t-1)
  theta_like_log <- 0 # f(x_t|theta_t)
  x_like_log <- 0 # f(y_t|x_t)
  
  
  # Leaf node
  ## Laten variables x's
  out_leaf <- dc_leaf(y = data, n = n, prior_mu = prior_mu, prior_sigma = prior_sigma, N = N)
  x_leaf <- out_leaf$x
  w_leaf_log <- out_leaf$w_log
  W_leaf <- out_leaf$W
  
  ## theta_0
  theta_0 <- rexp(N, 1)
  w_theta0_log <- apply(cbind(x_leaf[,1], theta_0), 1, function(x) dnorm(x[1], 0, x[2], log = T))
  #if (sum(exp(w_theta0_log)) == Inf) w_theta0_log <- w_theta0_log - max(w_theta0_log)
  #W_theta0 <- exp(w_theta0_log) / sum(exp(w_theta0_log))
  W_theta0 <- exp(w_theta0_log - matrixStats::logSumExp(w_theta0_log))
  
  # Some densities
  theta_0_prior_log <- sapply(theta_0, function(x) dexp(x, 1, log = T))
  theta_density_forward_log <- 0
  theta_like_log <- apply(cbind(x_leaf[,1], theta_0), 1, function(x) dnorm(x[1], 0, x[2], log = T))
  x_like_log <- sapply(x_leaf[,1], function(x) dnorm(data[1], x, 1, log = T))
  
  
  
  for (t in 1:(n-1)) {
    
    # No merging steps as p_check() is choosen to be prod(p_c())
    
    # SMC sampler
    if (t == 1) {
      alpha_update <- alpha
      
      ## Proposal
      theta_tilde <- sapply(theta_0, function(x) rexp(1, 1/x))
      
      w_log_0 <- 0
      for (i in 1:nt) {
        if (i == 1) {
          ## Update weights
          q_t_log <- apply(cbind(theta_tilde, theta_0), 1, function(x) dexp(x[1], 1/x[2], log = T))
          p_c_prod_log <- w_leaf_log[,t] + w_theta0_log
          gamma_t_0_log <- p_c_prod_log + q_t_log
          p_prev_log <- gamma_t_0_log
          
          gamma_t_nt_log <- theta_0_prior_log + theta_density_forward_log + theta_like_log + x_like_log 
          
          p_current_log <- (1-alpha_update)*gamma_t_0_log + alpha_update*gamma_t_nt_log
          
          w_log[,t] <- w_log_0 + p_current_log - p_prev_log
          #w_log[,t] <- weight_check(w_log = w_log[,t], min_lim = min_lim)
          #if (sum(exp(w_log[,t])) == Inf) w_log[,t] <- w_log[,t] - max(w_log[,t])
          #W[,t] <- exp(w_log[,t]) / sum(exp(w_log[,t]))
          
          W[,t] <- exp(w_log[,t] - matrixStats::logSumExp(w_log[,t]))
          
          
          ## Resampling
          U <- runif(1, 0, 1)
          A <- Sys_resamp(W = W[,t], P = N, U = U)
          theta_resamp <- theta_tilde[A]
          
          
          ## MCMC kernel
          theta[,t] <- mcmc_kernel(data = x_leaf[,t+1], node = t, theta_init = theta_tilde, prior_lambda = theta_0, N = N, Ntotal = Ntotal)
        }
        
        alpha_update <- alpha_update + alpha
        p_prev_log <- p_current_log
        
        ## Update weights
        q_t_log <- apply(cbind(theta[,t], theta_0), 1, function(x) dexp(x[1], 1/x[2], log = T))
        gamma_t_0_log <- p_c_prod_log + q_t_log
        
        #theta_density_forward_temp_log <- theta_density_forward_log
        #theta_density_forward_temp_log <- theta_density_forward_temp_log + sapply(cbind(theta[,t], theta_0), function(x) dexp(x[1], x[2], log = T))
        #gamma_t_nt_log <- theta_0_prior_log + theta_density_forward_temp_log + theta_like_log + x_like_log
        p_current_log <- (1-alpha_update)*gamma_t_0_log + alpha_update*gamma_t_nt_log
        
        w_log[,t] <- w_log_0 + p_current_log - p_prev_log ## Resampling for all iterations
        #w_log[,t] <- weight_check(w_log = w_log[,t], min_lim = min_lim)
        #if (sum(exp(w_log[,t])) == Inf) w_log[,t] <- w_log[,t] - max(w_log[,t])
        #W[,t] <- exp(w_log[,t]) / sum(exp(w_log[,t]))
        
        W[,t] <- exp(w_log[,t] - matrixStats::logSumExp(w_log[,t]))
        
        
        ## Resampling
        U <- runif(1, 0, 1)
        A <- Sys_resamp(W = W[,t], P = N, U = U)
        theta_resamp <- theta[A,t]
        
        
        ## MCMC kernel
        theta[,t] <- mcmc_kernel(data = x_leaf[,t+1], node = t, theta_init = theta_resamp, prior_lambda = theta_0, N = N, Ntotal = Ntotal)
      }
      
      # Update some densities
      theta_density_forward_log <- theta_density_forward_log + apply(cbind(theta[,t], theta_0), 1, function(x) dexp(x[1], x[2], log = T))
      theta_like_log <- theta_like_log + apply(cbind(x_leaf[,t+1], theta[,t]), 1, function(x) dnorm(x[1], 0, x[2], log = T))
      x_like_log <- x_like_log + sapply(x_leaf[,t+1], function(x) dnorm(data[t+1], x, 1, log = T))
    }
    
    # No merging steps as p_check() is choosen to be prod(p_c()
    
    # Reset alpha_update
    alpha_update <- alpha
    
    # Proposal
    theta_tilde <- sapply(theta[,t], function(x) rexp(1, 1/x))
    
    # Initialization
    w_log_0 <- 0
    
    for(i in 1:nt) {
      if (i == 1) {
        ## Update weights
        q_t_log <- apply(cbind(theta_tilde, theta[,t]), 1, function(x) dexp(x[1], 1/x[2], log = T))
        p_c_prod_log <- w_leaf_log[,t+1] + w_log[,t]
        gamma_t_0_log <- p_c_prod_log + q_t_log
        p_prev_log <- gamma_t_0_log
        
        gamma_t_nt_log <- theta_0_prior_log + theta_density_forward_log + theta_like_log + x_like_log 
        
        p_current_log <- (1-alpha_update)*gamma_t_0_log + alpha_update*gamma_t_nt_log
        
        w_log[,t+1] <- w_log_0 + p_current_log - p_prev_log
        #w_log[,t] <- weight_check(w_log = w_log[,t], min_lim = min_lim)
        #if (sum(exp(w_log[,t])) == Inf) w_log[,t] <- w_log[,t] - max(w_log[,t])
        #W[,t] <- exp(w_log[,t]) / sum(exp(w_log[,t]))
        W[,t+1] <- exp(w_log[,t+1] - matrixStats::logSumExp(w_log[,t+1]))
        
        ## Resampling
        U <- runif(1, 0, 1)
        A <- Sys_resamp(W = W[,t+1], P = N, U = U)
        theta_resamp <- theta_tilde[A]
        
        
        ## MCMC kernel
        theta[,t+1] <- mcmc_kernel(data = x_leaf[,t+1], node = t+1, theta_init = theta_resamp, prior_lambda = theta[,t], N = N, Ntotal = Ntotal)
      }
      
      alpha_update <- alpha_update + alpha
      p_prev_log <- p_current_log
      
      ## Update weights
      q_t_log <- apply(cbind(theta[,t+1], theta[,t]), 1, function(x) dexp(x[1], 1/x[2], log = T))
      gamma_t_0_log <- p_c_prod_log + q_t_log
      
      p_current_log <- (1-alpha_update)*gamma_t_0_log + alpha_update*gamma_t_nt_log
      
      w_log[,t+1] <- w_log_0 + p_current_log - p_prev_log ## Resampling for all iterations
      #w_log[,t] <- weight_check(w_log = w_log[,t], min_lim = min_lim)
      #if (sum(exp(w_log[,t])) == Inf) w_log[,t] <- w_log[,t] - max(w_log[,t])
      #W[,t] <- exp(w_log[,t]) / sum(exp(w_log[,t]))
      W[,t+1] <- exp(w_log[,t+1] - matrixStats::logSumExp(w_log[,t+1]))
      
      ## Resampling
      U <- runif(1, 0, 1)
      A <- Sys_resamp(W = W[,t+1], P = N, U = U)
      theta_resamp <- theta[A, t+1]
      
      
      ## MCMC kernel
      theta[,t+1] <- mcmc_kernel(data = x_leaf[,t+1], node = t+1, theta_init = theta_resamp, prior_lambda = theta[,t], N = N, Ntotal = Ntotal)
    }
    
    if (t < (n - 1)) {
      # Update some densities
      theta_density_forward_log <- theta_density_forward_log + apply(cbind(theta[,t+1], theta[,t]), 1, function(x) dexp(x[1], 1/x[2], log = T))
      theta_like_log <- theta_like_log + apply(cbind(x_leaf[t+2], theta[,t+1]), 1, function(x) dnorm(x[1], 0, x[2], log = T))
      x_like_log <- x_like_log + sapply(x_leaf[t+2], function(x) dnorm(data[t+2], x, 1, log = T))  
    }
    
  }
  
  
  return(list(theta=theta[,n], theta_resamp=theta_resamp))
  
}




debug(dc_smc_A1)

out <- dc_smc_A1(data=y, N=2000, alpha=.2, Ntotal=50, prior_mu=0, prior_sigma=1)


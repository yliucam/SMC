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






# Algorithm B2 for the log-normal hierarchical model toy example

dc_smc_algB2_new <- function(data,
                             N,
                             nodes_n,
                             #nt=nrow(data),
                             #n_batch,
                             n_trial,
                             beta_prior,
                             gamma_prior,
                             m,
                             alpha,
                             Ntotal_sub,
                             Ntotal_root) {
  
  ESS_min <- N / 2
  min_lim <- log(.Machine$double.xmin)
  
  n_leaf <- nodes_n[3] # Total number of leaf nodes
  n <- dim(data)[1] # Number of observations
  n_sub <- nodes_n[2] # Total number of subroot nodes
  
  nt <- 1 / alpha
  
  # Storage
  ## Leaf storage
  x_leaf <- array(rep(NA, N*n_leaf), dim = c(N, n_leaf))
  W_leaf <- array(rep(NA, N*n_leaf), dim = c(N, n_leaf))
  L_leaf_log <- rep(NA, n_leaf) # Store the logarithmic marginal likelihood separately
  
  ## Subroot storage
  x_sub <- array(rep(NA, N*n_sub), dim = c(N, n_sub))
  w_sub_log <- array(rep(NA, N*n_sub*nt), dim = c(N, n_sub, nt))
  W_sub <- array(rep(NA, N*n_sub), dim = c(N, n_sub))
  
  ## Root storage
  x_root <- rep(NA, N)
  x_root_array <- array(rep(NA, N*nt), dim = c(N, nt))
  w_root_log <- array(rep(NA, N*nt), dim = c(N, nt))
  W_root <- rep(NA, N)
  W_root_array <- array(rep(NA, N*nt), dim = c(N, nt))
  
  
  
  
  # Initialization
  L_leaf_prod_log <- 0 # The logarithmic marginal likelihood on the leaf node
  
  
  # Leaf nodes
  for (leaf_i in 1:n_leaf) {
    out_leaf <- dc_smc_leaf_binom(data=data[,leaf_i],
                                  alpha=beta_prior$alpha[leaf_i],
                                  beta=beta_prior$beta[leaf_i],
                                  n_trial=n_trial[leaf_i],
                                  N=N)
    x_leaf[,leaf_i] <- out_leaf$p_bin_post
    W_leaf[,leaf_i] <- out_leaf$W
    L_leaf_log[leaf_i] <- out_leaf$L_log
    
    L_leaf_prod_log <- L_leaf_prod_log + L_leaf_log[leaf_i]
  }
  
  
  
  # The merging step - leaf to subroot
  
  ## Store mN matchings for the leaf nodes
  x_leaf_mN <- array(rep(NA, m*N*n_leaf), dim = c(m*N, n_leaf))
  W_leaf_mN <- array(rep(NA, m*N*n_leaf), dim = c(m*N, n_leaf))
  
  # Initilization
  L_prod_log <- L_leaf_prod_log
  
  
  ## Reampling mN particles from each leaf node
  for (leaf_i in 1:n_leaf) {
    U <- runif(1, 0, 1)
    A <- Sys_resamp(W=W_leaf[,leaf_i], P=m*N, U=U)
    x_leaf_mN[,leaf_i] <- x_leaf[A, leaf_i]
    W_leaf_mN[,leaf_i] <- W_leaf[A, leaf_i]
  }
  
  ## No merging step -- use the results from the child nodes directly
  
  p_check_log_resamp <- array(rep(NA, N*n_sub), dim = c(N, n_sub))
  x_leaf_realloc <- array(rep(NA, N*(n_leaf/n_sub)*n_sub), dim = c(N, n_leaf/n_sub, n_sub))
  
  for (sub_i in 1:n_sub) {
    x_leaf_realloc[,,sub_i] <- x_leaf[,((sub_i-1)*n_leaf/n_sub+1):(sub_i*n_leaf/n_sub)]
    p_check_log_resamp[,sub_i] <- rowMeans(W_leaf[,((sub_i-1)*n_leaf/n_sub+1):(sub_i*n_leaf/n_sub)])
  }
  
  
  
  # Subroot
  
  
  for (sub_i in 1:n_sub) {
    ## Initialization
    x_sub_0 <- runif(N, 0.001, 10000)
    w_log_0 <- rep(0, N)
    
    
    ## SMC sampler
    alpha_update <- alpha
    for (i in 1:(nt-1)) {
      if (i == 1) {
        
        ### Update weights
        #q_sub_log <- sapply(x_sub_0, function(x) dunif(x, .001, 10000, log = T))
        #p_sub_0_log <- p_check_log_resamp[,sub_i] + q_sub_log
        p_sub_0_log <- p_check_log_resamp[,sub_i]
        p_sub_like_log <- likelihood_particle_beta(obs = x_leaf_realloc[,,sub_i],
                                                   alpha = matrix(rep(x_sub_0, n_leaf/n_sub), ncol = n_leaf/n_sub),
                                                   beta = matrix(rep(beta_prior$beta[((sub_i-1)*n_leaf/n_sub+1):(sub_i*n_leaf/n_sub)], N), ncol=n_leaf/n_sub, byrow=T))
        #p_sub_prior_log <- q_sub_log # Note: at the first iteration, the prior is the same as q_sub in this example!
        #p_sub_t_log <- p_sub_prior_log + p_sub_like_log 
        p_sub_t_log <- p_sub_like_log
        
        p_sub_current <-  (1-alpha_update)*p_sub_0_log + alpha_update*p_sub_t_log
        w_sub_log[,sub_i,i] <- w_log_0 + p_sub_current - p_sub_0_log
        w_sub_log[which(w_sub_log[,sub_i,i] < min_lim), sub_i, i] <- min_lim
        
        W_sub[,sub_i] <- exp(w_sub_log[,sub_i,i] - matrixStats::logSumExp(w_sub_log[,sub_i,i]))
        
        
        ### Resampling
        ESS <- 1 / sum(W_sub[,sub_i]^2)
        if (ESS < ESS_min) {
          U <- runif(1, 0, 1)
          A <- Sys_resamp(W=W_sub[,sub_i], P=N, U=U)
          x_sub_resamp <- x_sub_0[A]  
        } else {
          x_sub_resamp <- x_sub_0
        }
        
        ### MCMC kernel
        #var_sub <- var(x_sub_resamp)
        x_sub[,sub_i] <- sub_mcmc(data = x_leaf_realloc[,,sub_i],
                                  x_0 = x_sub_resamp,
                                  like_beta = beta_prior$beta[((sub_i-1)*n_leaf/n_sub+1):(sub_i*n_leaf/n_sub)],
                                  var_sub = 1,
                                  Ntotal = Ntotal_sub)
        
      }
      
      ### Update p_sub_prev
      #q_sub_log <- sapply(x_sub[,sub_i], function(x) dunif(x, .001, 10000, log = T))
      #p_sub_0_log <- p_check_log_resamp[,sub_i] + q_sub_log
      p_sub_0_log <- p_check_log_resamp[,sub_i]
      p_sub_like_log <- likelihood_particle_beta(obs = x_leaf_realloc[,,sub_i],
                                                 alpha = matrix(rep(x_sub[,sub_i], n_leaf/n_sub), ncol = n_leaf/n_sub),
                                                 beta = matrix(rep(beta_prior$beta[((sub_i-1)*n_leaf/n_sub+1):(sub_i*n_leaf/n_sub)], N), ncol=n_leaf/n_sub, byrow=T))
      #p_sub_prior_log <- sapply(x_sub[,sub_i], function(x) dgamma(x, gamma_prior$alpha[sub_i], gamma_prior$beta[sub_i], log = T))
      #p_sub_prior_log <- q_sub_log
      #p_sub_t_log <- p_sub_prior_log + p_sub_like_log
      p_sub_t_log <- p_sub_like_log
      p_sub_prev <-  (1-alpha_update)*p_sub_0_log + alpha_update*p_sub_t_log
      
      
      ### Update alpha_i
      alpha_update <- alpha_update + alpha
      
    
      ### Update weights
      p_sub_current <- (1-alpha_update)*p_sub_0_log + alpha_update*p_sub_t_log
      if (ESS < ESS_min) {
        w_sub_log[,sub_i,i+1] <- 0 + p_sub_current - p_sub_prev # After the resampling, the previous w_log is 0
      } else {
        w_sub_log[,sub_i,i+1] <- w_sub_log[,sub_i,i] + p_sub_current - p_sub_prev
      }
      w_sub_log[which(w_sub_log[,sub_i,i+1] < min_lim), sub_i, i+1] <- min_lim
      
      W_sub[,sub_i] <- exp(w_sub_log[,sub_i,i+1] - matrixStats::logSumExp(w_sub_log[,sub_i,i+1]))
      
      ### Resampling -- optionally
      ESS <- 1 / sum(W_sub[,sub_i]^2)
      if (ESS < ESS_min) {
        U <- runif(1, 0, 1)
        A <- Sys_resamp(W=W_sub[,sub_i], P=N, U=U)
        x_sub_resamp <- x_sub[A,sub_i]
      } else {
        x_sub_resamp <- x_sub[,sub_i]
      }
      
      
      ### Marginal likelihood update
      L_prod_log <- L_prod_log + matrixStats::logSumExp(w_sub_log[,sub_i,i+1]) - log(N)
      
      
      ### MCMC kernel
      var_sub <- var(x_sub_resamp)
      x_sub[,sub_i] <- sub_mcmc(data = x_leaf_realloc[,,sub_i],
                                x_0 = x_sub_resamp,
                                like_beta = beta_prior$beta[((sub_i-1)*n_leaf/n_sub+1):(sub_i*n_leaf/n_sub)],
                                var_sub = var_sub,
                                Ntotal = Ntotal_sub)
    }
  }
  
  
  # No merging step - use the results from the subroot nodes directly
  
  p_check_log_resamp <- rowMeans(log(W_sub))
  
  #L_prod_log <- L_prod_log + matrixStats::logSumExp(v_t_log) - log(m*N)
  
  
  
  # Root
  ## Initialization
  x_root_0 <- runif(N, .001, 10000)
  w_log_0 <- rep(0, N)
  
  
  ## SMC sampler
  alpha_update <- alpha
  for (i in 1:(nt-1)) {
    if (i == 1) {
      
      ### Update weights
      p_root_0_log <- p_check_log_resamp
      p_root_like_log <- likelihood_particle_gamma(obs=x_sub,
                                                   alpha=cbind(x_root_0, x_root_0),
                                                   beta=matrix(rep(gamma_prior$beta, N), ncol=2, byrow=T))
      #p_root_prior_log <- q_root_log
      p_root_t_log <- p_root_like_log
      p_root_current <- (1-alpha_update)*p_root_0_log + alpha_update*p_root_t_log 
      w_root_log[,i] <- w_log_0 + p_root_current - p_root_0_log
      
      w_root_log[which(w_root_log[,i] < min_lim), i] <- min_lim
      
      W_root <- exp(w_root_log[, i] - matrixStats::logSumExp(w_root_log[, i]))
      
      W_root_array[, i] <- W_root
      
      
      ### Resampling
      ESS <- 1 / sum(W_root^2)
      if (ESS < ESS_min) {
        U <- runif(1, 0, 1)
        A <- Sys_resamp(W=W_root, P=N, U=U)
        x_root_resamp <- x_root_0[A]
      } else {
        x_root_resamp <- x_root_0
      }
      
      ### MCMC kernel
      #var_root <- var(x_root_resamp)
      x_root <- root_mcmc(data = x_sub,
                          x_0 = x_root_resamp,
                          like_beta = gamma_prior$beta,
                          var_root = 1,
                          Ntotal = Ntotal_root)
      x_root_array[,i] <- x_root
    }
    
    ### Update p_root_prev
    p_root_0_log <- p_check_log_resamp
    p_root_like_log <- likelihood_particle_gamma(obs = x_sub,
                                                 alpha = cbind(x_root, x_root),
                                                 beta = matrix(rep(gamma_prior$beta, N), ncol=2, byrow=T))
    #p_root_prior_log <- q_root_log
    p_root_t_log <- p_root_like_log
    p_root_prev <- (1-alpha_update)*p_root_0_log + alpha_update*p_root_t_log
    
    ### Update alpha_i
    alpha_update <- alpha_update + alpha
    
    ### Update weights
    p_root_current <- (1-alpha_update)*p_root_0_log + alpha_update*p_root_t_log
    if (ESS < ESS_min) {
      w_root_log[,i+1] <- 0 + p_root_current - p_root_prev # The previous w_log is 0 after the resampling
    } else {
      w_root_log[,i+1] <- w_root_log[,i] + p_root_current - p_root_prev
    }
    w_root_log[which(w_root_log[,i+1] < min_lim), i+1] <- min_lim
   
    W_root <- exp(w_root_log[,i+1] - matrixStats::logSumExp(w_root_log[,i+1]))

    W_root_array[, i+1] <- W_root
    
    
    ### Resampling -- optionally
    ESS <- 1 / sum(W_root^2)
    if (ESS < ESS_min) {
      U <- runif(1, 0, 1)
      A <- Sys_resamp(W=W_root, P=N, U=U)
      x_root_resamp <- x_root[A]
    } else {
      x_root_resamp <- x_root
    }
    
    
    
    ### Marginal likelihood update
    L_prod_log <- L_prod_log + matrixStats::logSumExp(w_root_log[,i+1]) - log(N)
    
    ### MCMC kernel
    var_root <- max(50, var(x_root_resamp))
    if (isTRUE(all.equal(0, var_root))){
      browser()
    }
    x_root <- root_mcmc(data = x_sub,
                        x_0 = x_root_resamp,
                        like_beta = gamma_prior$beta,
                        var_root = 1,
                        Ntotal = Ntotal_root)
    x_root_array[,i+1] <- x_root
    
  }
  
  
  return(list(x_leaf = x_leaf, x_sub = x_sub, x_root = x_root, x_root_array = x_root_array, W_root_array = W_root_array))
}


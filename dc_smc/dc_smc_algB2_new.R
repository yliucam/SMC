dc_smc_algB2_new <- function(data,
                         n_trial,
                         N,
                         nt=nrow(data),
                         n_batch,
                         gamma_prior,
                         m,
                         alpha,
                         Ntotal_sub) {
  
  ESS_min <- P / 2
  min_lim <- log(.Machine$double.xmin)
  
  n_leaf <- dim(data)[2] # Total number of leaf nodes
  n <- dim(data)[1] # Number of observations
  n_sub <- n_leaf / 2 # Total number of subroot nodes
  
  alpha_inc <- 1 / nt
  
  # Storage
  ## Leaf storage
  x_leaf <- array(rep(NA, N*n_leaf), dim = c(N, n_leaf))
  W_leaf <- array(rep(NA, N*n_leaf), dim = c(N, n_leaf))
  L_leaf_log <- array(rep(NA, N*n_leaf), dim = c(N, n_leaf)) # Store the logarithmic marginal likelihood separately
  
  ## Subroot storage
  x_sub <- array(rep(NA, N*n_sub), dim = c(N, n_sub))
  w_sub_log <- array(rep(NA, N*n_sub*nt), dim = c(N, n_sub, nt))
  W_sub <- array(rep(NA, N*n_sub), dim = c(N, n_sub))
  
  
  
  
  # Initialization
  L_leaf_prod_log <- 0 # The logarithmic marginal likelihood on the leaf node 
  
  
  # Leaf nodes
  for (leaf_i in 1:n_leaf) {
    out_leaf <- dc_smc_leaf_binom(data=data[,leaf_i],
                                      alpha=.5,
                                      beta=.5,
                                      n_trial=n_trial[leaf_i],
                                      N=N)
    x_leaf[,leaf_i] <- out_leaf$p_bin_post
    W_leaf[,leaf_i] <- out_leaf$W
    L_leaf_log[,leaf_i] <- out_leaf$L_log
    
    L_leaf_prod_log <- L_leaf_prod_log + L_leaf_log[,leaf_i]
  }
  
  
  
  # The merging step - leaf to subroot
  
  ## mN matchings storage for leaf
  x_leaf_mN <- array(rep(NA, m*N*n_leaf), dim = c(m*N, n_leaf))
  W_leaf_mN <- array(rep(NA, m*N*n_leaf), dim = c(m*N, n_leaf))
  
  
  ## Reampling mN particles from each leaf node
  for (leaf_i in 1:n_leaf) {
    U <- runif(1, 0, 1)
    A <- Sys_resamp(W=W_leaf[,leaf_i], P=m*N, U=U)
    x_leaf_mN[,leaf_i] <- x_leaf[A, leaf_i]
    W_leaf_mN[,leaf_i] <- W_leaf[A, leaf_i]
  }
  
  ## Resampling N matchings from the mN matchings above based on V_t
  
  beta_shape2 <- c(2, 10, 2, 10) # Known parameters for the beta prior
  
  v_t_log <- array(rep(NA, m*N*n_sub), dim = c(m*N, n_sub))
  p_check_log <- array(rep(NA, m*N*n_sub), dim = c(m*N, n_sub))
  x_leaf_resamp <- array(rep(NA, N*2*n_sub), dim = c(N, 2, n_sub)) # Store the N resampling matchings
  p_check_log_resamp <- array(rep(NA, N*n_sub), dim = c(N, n_sub)) # Store the N corresponding log pi_check
  
  for (sub_i in 1:n_sub) {
    mu <- gamma_prior$alpha[sub_i] / gamma_prior$beta[sub_i]
    for (j in 1:(m*N)) {
      pc_log <- sum(log(W_leaf_mN[j,((sub_i-1)*2+1):(i*2)]))
      pc_int_log <- dbeta(x_leaf_mN[j,((sub_i-1)*2+1)], mu, beta_shape2[((sub_i-1)*2+1)], log=T) +
        dbeta(x_leaf_mN[j,(sub_i*2)], mu, beta_shape2[sub_i*2], log=T)
      p_check_log[j,sub_i] <- (1-alpha)*pc_log + alpha*pc_int_log
      v_t_log[j,sub_i] <- p_check_log[j,sub_i] - pc_log
      if (v_t_log[j,sub_i] < log(.Machine$double.xmin)) v_t_log[j,sub_i] <- log(.Machine$double.xmin)
    }
    if (sum(exp(v_t_log[,sub_i])) == Inf) v_t_log[,sub_i] <- v_t_log[,sub_i] - max(v_t_log[,sub_i])
    V_t <- exp(v_t_log[,sub_i])/sum(exp(v_t_log[,sub_i]))
    
    U <- runif(1, 0, 1)
    A <- Sys_resamp(W=V_t, P=N, U=U)
    x_leaf_resamp[,,sub_i] <- x_leaf_mN[A,((sub_i-1)*2+1):(sub_i*2)]
    p_check_log_resamp[,sub_i] <- p_check_log[A,sub_i]
    
    L_prod_log <- L_leaf_prod_log + log(sum(exp(v_t_log[,sub_i]))) - log(m*N)
  }
  
  
  
  # Subroot
  
  
  for (sub_i in 1:n_sub) {
    ## Initialization
    x_sub_0 <- rgamma(N, gamma_prior$alpha[sub_i], gamma_prior$beta[sub_i])
    w_log_0 <- rep(0, N)
    
    
    ## SMC sampler
    alpha_update <- alpha_inc
    for (i in 1:(nt-1)) {
      if (i == 1) {
        q_sub_log <- sapply(x_sub_0, function(x) dgamma(x, gamma_prior$alpha[sub_i], gamma_prior$beta[sub_i], log = T))
        p_sub_0_log <- p_check_log_resamp[,sub_i] + q_sub
        p_sub_like_log <- likelihood_particle_beta(obs = x_leaf_resamp[,,sub_i],
                                               alpha = cbind(x_sub_0, x_sub_0),
                                               beta = matrix(rep(beta_shape2[((sub_i-1)*2+1):(i*2)], N), ncol=2, byrow=T))
        p_sub_t_log <- q_sub_log + p_sub_like_log # Note: at the first iteration, the prior is the same as q_sub in this example!
        
        ### Update weights
        p_sub_current <-  (1-alpha_update)*p_sub_0_log + alpha_update*p_sub_t_log # Better to store here for the use by the next iteration!
        w_sub_log[,sub_i,i] <- w_log_0 + p_sub_current - p_sub_0_log
        w_sub_log[,sub_i,i] <- weight_check(w_log = w_sub_log[,sub_i,i], min_lim = min_lim)
        
        if (max(w_sub_log[,sub_i,i]) == min_lim) w_sub_log[,sub_i,i] <- 0
        W_sub[,sub_i] <- exp(w_sub_log[,sub_i,i]) / sum(exp(w_sub_log[,sub_i,i]))
        
        ### Resampling
        U <- runif(1, 0, 1)
        A <- Sys_resamp(W=W_sub[,sub_i], P=N, U=U)
        x_sub_resamp <- x_sub_0[A]
        
        ### Marginal likelihood update
        L_prod_log <- L_prod_log + log(sum(exp(w_sub_log[,sub_i,i]))) - log(N)
        
        ### MCMC kernel
        x_sub[,sub_i] <- sub_mcmc(data = x_leaf_resamp[,,sub_i],
                                  x_0 = x_sub_resamp,
                                  prior_alpha = gamma_prior$alpha[sub_i],
                                  prior_beta = gamma_prior$beta[sub_i],
                                  like_beta = beta_shape2[((sub_i-1)*2+1):(i*2)],
                                  Ntotal = Ntotal_sub)
        
      }
      
      alpha_update <- alpha_update + alpha_inc
      p_sub_prev <- p_sub_current
      
      ## Note: p_sub_0 keeps the same for all iterations!
      ## The prior for new samples is needed -- note: q_sub_log is not the item used here!
      p_sub_prior_log <- sapply(x_sub[,sub_i], function(x) dgamma(x, gamma_prior$alpha[sub_i], gamma_prior$beta[sub_i], log = T))
      p_sub_like_log <- likelihood_particle_beta(obs = x_leaf_resamp[,,sub_i],
                                             alpha = cbind(x_sub[,sub_i], x_sub[,sub_i]),
                                             beta = matrix(rep(beta_shape2[((sub_i-1)*2+1):(i*2)], N), ncol=2, byrow=T))
      p_sub_t_log <- p_sub_prior_log + p_sub_like_log
      
      ### Update weights
      p_sub_current <- (1-alpha_update)*p_sub_0_log + alpha_update*p_sub_t_log
      w_sub_log[,sub_i,i+1] <- 0 + p_sub_current - p_sub_prev # After the resampling, the previous w_log is 0
      w_sub_log[,sub_i,i+1] <- weight_check(w_log = w_sub[,sub_i,i+1], min_lim = min_lim)
      
      if (max(w_sub_log[,sub_i,i+1]) == min_lim) w_sub_log[,sub+i,i+1] <- 0
      W_sub[,sub_i] <- exp(w_sub_log[,sub_i,i+1]) / sum(exp(w_sub_log[,sub_i,i+1]))
      
      ### Resampling
      U <- runif(1, 0, 1)
      A <- Sys_resamp(W=W_sub[,sub_i], P=N, U=U)
      x_sub_resamp <- x_sub[A,sub_i]
      
      ### Marginal likelihood update
      L_prod_log <- L_prod_log + log(sum(exp(w_sub_log[,sub_i,i+1]))) - log(N)
      
      ### MCMC kernel
      x_sub[,sub_i] <- sub_mcmc(data = x_leaf_resamp[,,sub_i],
                                x_0 = x_sub_resamp,
                                prior_alpha = gamma_prior$alpha[sub_i],
                                prior_beta = gamma_prior$beta[sub_i],
                                like_beta = beta_shape2[((sub_i-1)*2+1):(i*2)],
                                Ntotal = Ntotal_sub)
      
    }
  }
  
  
  
  
  
  
  
  
}









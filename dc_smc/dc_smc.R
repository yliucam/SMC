dc_smc <- function(data, P, nt) {
  # Try the original Algorithm 2 first
  
  #ESS_min <- P / 2 ## Initially, do not use this adaptive implementation
  
  ## Try a toy example with TWO child nodes
  TT_c <- dim(data)[1] ## Assume the data are stored in TWO columns in the example
  n_n <- dim(data)[2]  ## The nodes have the same number of observations
  
  # Storage
  #x_1 <- x_2 <- array(rep(NA, Ntotal*P), dim = c(P, Ntotal))
  #x <- list(x1, x2)
  x <- array(rep(NA, P*n_n), dim = c(P, n_n))
  #x_1_mP <- x_2_mP <- rep(NA, m*P)
  #x_mP <- cbind(x_1_mP, x_2_mP)
  
  w_log <- array(rep(NA, P*n_n), dim = c(P, n_n))
  w_t_log <- rep(0, P)
  
  # Initialization
  x_0 <- rep(0, n_n)
  w_log_0 <- 0
  
  for (c_i in 1:n_n) {
    x_ci <- x[,n_n]
    x_ci[1] <- x_0[c_i] + rnorm(1, 0, 1)
    w_log[1,c_i] <- w_log_0 + sum(sapply(data[,c_i], function(y) dnorm(y, x_ci[1], 1, log = T)))
    for (j in 2:P) {
      x_ci[j] <- x_ci[j-1] + rnorm(1, 0, 1)
      w_log[j,c_i] <- w_log[j-1,c_i] + sum(apply(data[,c_i], 1, function(y) dnorm(y, x_ci[j], 1, log = T)))
    }
    W <- exp(w_log[j,c_i]) / sum(exp(w_log[j,c_i]))
    ### Independently resampling N particles -- 1(b)
    U <- runif(1, 0, 1)
    A <- Sys_resamp(W, P, U)
    x[,n_n,i] <- x_ci[A]
    ### Resampling mP particles -- 1(b) in Algorithm B2
    #for (jj in 1:m) {
    #  U <- runif(1, 0, 1)
    #  A <- Sys_resamp(W, P, U)
    #  x_mP[((jj-1)*P+1):(jj*P), c_i] <- x_ci[A,i]
    #}
  }
  
  
  
  # Propose x_t from q_t(|x_c1, x_c2,...). Not consider \tilde{X}_t is empty for now
  ## How? Use a conjugate distribution with mean of (x_c1, x_c2,...)?
  x_t_tilde <- rowMeans(x) + rnorm(3, 0, 1) ## Using a Normal random walk centered on the mean of x_c1, x_c2,...
  
  # Sampler steps
  
  ## Storage
  x_sampler <- array(rep(NA, P*nt), dim = c(P, nt))
  w_sampler_log <- array(rep(NA, P*nt), dim = c(P, nt))
  
  ## Initialization
  x_sampler_0 <- x_t_tilde
  w_sampler_log_0 <- 0
  
  ## Sampler iterations
  for (s_i in 1:(nt-1)) {
    if (s_i == 1) {
      for (j in 1:P) {
        w_sampler_log[j,s_i] <- w_sampler_log_0 + sum(sapply(x[j,], function(y) dnorm(y, x_sampler_0[j], 1, log = T)))
      }
      W <- exp(w_sampler_log[j,s_i]) / sum(exp(w_sampler_log[j,s_i]))
      
      ### Resampling
      U <- runif(1, 0, 1)
      A <- Sys_resamp(W, P, U)
      x_sampler_0_re <- x_sampler_0[A]
      
      ### Proposals
      x_sampler[,s_i] <- x_sampler_0 + rnorm(P, 0, 1)
    }
    
    for (j in 1:P) {
      w_sampler_log[j,s_i+1] <- w_sampler_log[j,s_i] + sum(sapply(x[j,], function(y) dnorm(y, x_sampler[j,s_i], 1, log = T)))
    }
    W <- exp(w_sampler_log[j,s_i+1]) / sum(exp(w_sampler_log[j,s_i+1]))
    
    ### Resampling
    U <- runif(1, 0, 1)
    A <- Sys_resamp(W, P, U)
    x_sampler_re <- x_sampler[A,s_i]
    
    ### Proposals
    x_sampler[,s_i+1] <- x_sampler_re + rnorm(P, 0, 1)
  }
}

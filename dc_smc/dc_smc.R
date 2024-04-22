dc_smc <- function(data, P, nt=nrow(data)) {
  # Try the original Algorithm 2 first
  
  #ESS_min <- P / 2 ## Initially, do not use this adaptive implementation
  
  ## Try a toy example with TWO child nodes
  TT_c <- dim(data)[1] ## Assume the data are stored in TWO columns in the example
  n_n <- dim(data)[2]  ## The nodes have the same number of observations
  
  # Storage
  x_c <- array(rep(NA, P*nt), dim = c(P, nt))
  x_c_nt <- matrix(rep(NA, P*n_n), ncol = n_n)
  x <- array(rep(NA, P*n_n), dim = c(P, n_n))
  #x_1_mP <- x_2_mP <- rep(NA, m*P)
  #x_mP <- cbind(x_1_mP, x_2_mP)
  
  w_log <- array(rep(NA, P*nt), dim = c(P, nt))
  w_t_log <- array(rep(NA, P*nt), dim = c(P, nt))
  
  # Initialization
  x_c_0 <- rep(0, P)
  w_log_0 <- 0
  L_c_prod_log <- 0
  
  # Child nodes
  for (c_i in 1:n_n) {
    for (i in 1:(nt-1)) {
      if (i == 1) {
        x_c[,i] <- x_c_0 + rnorm(P, 0, 1000)
        for (j in 1:P) {
          w_log[j,i] <- w_log_0 + sum(sapply(data[,c_i], function(y) dnorm(y, x_c[j,i], 1, log = T)))
          if (w_log[j,i] < log(.Machine$double.xmin)) w_log[j,i] <- log(.Machine$double.xmin)
        }
        W <- exp(w_log[,i]) / sum(exp(w_log[,i]))
        
        U <- runif(1, 0, 1)
        A <- Sys_resamp(W=W, P=P, U=U)
        x_c[,i] <- x_c[A,i]
        
        L_c_prod_log <- L_c_prod_log + log(sum(exp(w_log[,i]))) - log(P)
      }
      
      x_c[,i+1] <- x_c[,i] + rnorm(P, 0, 1000)
      for (j in 1:P) {
        w_log[j,i+1] <- sum(sapply(data[,c_i], function(y) dnorm(y, x_c[j,i+1], 1, log = T)))
        if (w_log[j,i+1] < log(.Machine$double.xmin)) w_log[j,i+1] <- log(.Machine$double.xmin)
      }
      W <- exp(w_log[,i+1]) / sum(exp(w_log[,i+1]))
      
      U <- runif(1, 0, 1)
      A <- Sys_resamp(W=W, P=P, U=U)
      x_c[,i+1] <- x_c[A,i+1]
      
      L_c_prod_log <- L_c_prod_log + log(sum(exp(w_log[,i+1]))) - log(P)
    }
    
    x_c_nt[,c_i] <- x_c[,nt]
  }
  
  
  
  # Propose x_t from q_t(|x_c1, x_c2,...). 
  ## Choose q_t to be u_t, the prior for \tilde{x}_t
  ## Select u_t as a conjugate prior to p(x_c|\tilde{x}_t) ~ N(\tilde{x}_t, 1).
  ## Due to above, we choose u_t ~ N(0, 1). Then q_t ~ N((x_c1+x_c2+...)/(1+n_n), 1/(1+n_n)).
  x_t_tilde <- rowSums(x_c_nt)/(1+n_n) + sqrt(1/(1+n_n)) * rnorm(P, 0, 1)
  
  
  # Root
  
  ## Initialization
  x_0 <- x_t_tilde
  w_t_log_0 <- 0
  L_prod_log <- L_c_prod_log
  
  ## Sampler iterations
  for (i in 1:(nt-1)) {
    if (i == 1) {
      x[,i] <- x_0 + rnorm(P, 0, 1)
      for (j in 1:P) {
        w_t_log[j,i] <- w_t_log_0 + sum(sapply(x_c_nt[j,], function(y) dnorm(y, x[j,i], 1, log = T))) ## Should it be x_c_nt[j,] or x_c_nt?
        if (w_t_log[j,i] < log(.Machine$double.xmin)) w_t_log[j,i] <- log(.Machine$double.xmin)
      }
      W <- exp(w_t_log[,i]) / sum(exp(w_t_log[,i]))
      
      ### Resampling
      U <- runif(1, 0, 1)
      A <- Sys_resamp(W, P, U)
      x_sampler_0_re <- x_sampler_0[A]
      
      ### Proposals
      x_sampler[,s_i] <- x_sampler_0 + rnorm(P, 0, 1000)
      
      ### Update the marginal likelihood
      L_prod_log <- L_prod_log + L_c_prod_log
    }
    
    for (j in 1:P) {
      w_sampler_log[j,s_i+1] <- w_sampler_log[j,s_i] + sum(sapply(x[j,], function(y) dnorm(y, x_sampler[j,s_i], 1, log = T)))
      if (w_sampler_log[j,s_i+1] < log(.Machine$double.xmin)) w_sampler_log[j,s_i+1] <- log(.Machine$double.xmin)
    }
    W <- exp(w_sampler_log[,s_i+1]) / sum(exp(w_sampler_log[,s_i+1]))
    
    ### Resampling
    U <- runif(1, 0, 1)
    A <- Sys_resamp(W, P, U)
    x_sampler_re <- x_sampler[A,s_i]
    
    ### Update the marginal likelihood
    L_prod_log <- L_prod_log + L_prod_log
    
    ### Proposals
    x_sampler[,s_i+1] <- x_sampler_re + rnorm(P, 0, 1000)
  }
  
  return(list(x_sampler=x_sampler[,nt], x_c=x))
}

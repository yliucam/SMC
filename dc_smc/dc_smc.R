dc_smc <- function(data, P, nt=nrow(data)) {
  # Try the original Algorithm 2 first
  
  #ESS_min <- P / 2 ## Initially, do not use this adaptive implementation
  
  ## Try a toy example with TWO child nodes
  TT_c <- dim(data)[1] ## Assume the data are stored in TWO columns in the example
  n_n <- dim(data)[2]  ## The nodes have the same number of observations
  
  # Storage
  x_c <- array(rep(NA, P*nt), dim = c(P, nt))
  x_c_nt <- matrix(rep(NA, P*n_n), ncol = n_n)
  x <- array(rep(NA, P*nt), dim = c(P, nt))
  #x_1_mP <- x_2_mP <- rep(NA, m*P)
  #x_mP <- cbind(x_1_mP, x_2_mP)
  
  w_log <- array(rep(NA, P*nt), dim = c(P, nt))
  w_t_log <- array(rep(NA, P*nt), dim = c(P, nt))
  
  # Initialization
  x_c_0 <- rep(0, P)
  w_log_0 <- 0
  L_c_prod_log <- 0
  
  
  # NOTE: the updates follow the order in the regular partical filter (see e.g. Algorithm 10.4 in Chopin (2020)).
  # This is in a different updating order from the dc_smc sampler proposed by Lindsten (2016) (see its Algorithem B1 and B2)
  # This leads to the question on whether the final outputs of x should be resampled?
  
  # Child nodes
  for (c_i in 1:n_n) {
    ## Sampler iterations
    for (i in 1:(nt-1)) {
      if (i == 1) {
        ### Resampling
        U <- runif(1, 0, 1)
        A <- Sys_resamp(W=rep(1/P,P), P=P, U=U)
        x_c_0 <- x_c_0[A]
        
        ### Proposals
        x_c[,i] <- x_c_0 + rnorm(P, 0, 10)
        
        ### Update weights
        for (j in 1:P) {
          w_log[j,i] <- w_log_0 + sum(sapply(data[,c_i], function(y) dnorm(y, x_c[j,i], 1, log = T)))
          if (w_log[j,i] < log(.Machine$double.xmin)) w_log[j,i] <- log(.Machine$double.xmin)
        }
        W <- exp(w_log[,i]) / sum(exp(w_log[,i]))
        
        ### Update the marginal likelihood
        L_c_prod_log <- L_c_prod_log + log(sum(exp(w_log[,i]))) - log(P)
      }
      
      ### Resampling
      U <- runif(1, 0, 1)
      A <- Sys_resamp(W=W, P=P, U=U)
      x_c[,i] <- x_c[A,i]
      
      ### Proposals
      x_c[,i+1] <- x_c[,i] + rnorm(P, 0, 10)
      
      ### Update weights
      for (j in 1:P) {
        w_log[j,i+1] <- sum(sapply(data[,c_i], function(y) dnorm(y, x_c[j,i+1], 1, log = T)))
        if (w_log[j,i+1] < log(.Machine$double.xmin)) w_log[j,i+1] <- log(.Machine$double.xmin)
      }
      W <- exp(w_log[,i+1]) / sum(exp(w_log[,i+1]))
      
      ### Update the marginal likelihood
      L_c_prod_log <- L_c_prod_log + log(sum(exp(w_log[,i+1]))) - log(P)
    }
    
    ### Is resampling needed for the final output of x?
    #U <- runif(1, 0, 1)
    #A <- Sys_resamp(W=W, P=P, U=U)
    #x_c[,nt] <- x_c[A,nt]
    
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
      ### Resampling
      U <- runif(1, 0, 1)
      A <- Sys_resamp(W=rep(1/P,P), P=P, U=U)
      x_0 <- x_0[A]
      
      ### Proposals
      x[,i] <- x_0 + rnorm(P, 0, 10)
      
      ### Update weights
      for (j in 1:P) {
        w_t_log[j,i] <- w_t_log_0 + sum(sapply(x_c_nt[j,], function(y) dnorm(y, x[j,i], 1, log = T))) ## Should it be x_c_nt[j,] or x_c_nt?
        if (w_t_log[j,i] < log(.Machine$double.xmin)) w_t_log[j,i] <- log(.Machine$double.xmin)
      }
      W <- exp(w_t_log[,i]) / sum(exp(w_t_log[,i]))
      
      ### Update the marginal likelihood
      L_prod_log <- L_prod_log + log(sum(exp(w_t_log[,i]))) - log(P)
    }
    
    ### Resampling
    U <- runif(1, 0, 1)
    A <- Sys_resamp(W=W, P=P, U=U)
    x[,i] <- x[A,i]
    
    ### Proposals
    x[,i+1] <- x[,i] + rnorm(P, 0, 10)
    
    ### Update weights
    for (j in 1:P) {
      w_t_log[j,i+1] <- sum(sapply(x_c_nt[j,], function(y) dnorm(y, x[j,i+1], 1, log = T))) ## Same question as above
      if (w_t_log[j,i+1] < log(.Machine$double.xmin)) w_t_log[j,i+1] <- log(.Machine$double.xmin)
    }
    W <- exp(w_t_log[,i+1]) / sum(exp(w_t_log[,i+1]))
    
    ### Update the marginal likelihood
    L_prod_log <- L_prod_log + log(sum(exp(w_t_log[,i+1]))) - log(P)
    
  }
  
  ### Is resampling needed for the final output of x?
  #U <- runif(1, 0, 1)
  #A <- Sys_resamp(W=W, P=P, U=U)
  #x[,nt] <- x[A,nt]
  
  return(list(x=x[,nt], x_c_nt=x_c_nt, x_t_tilde=x_t_tilde, W=W))
}


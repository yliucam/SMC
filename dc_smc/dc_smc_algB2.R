#
#---- Algorithm B2------------------------#
#


# Toy example
mu <- rnorm(1, 100, 1)
mu1 <- rnorm(1, mu, 1)
mu2 <- rnorm(1, mu, 1)

y1 <- rnorm(100, mu1, 1)
y2 <- rnorm(100, mu2, 1)
data <- cbind(y1, y2)




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




dc_smc_algB2 <- function(data, P, nt=nrow(data), m, alpha) {
  ESS_min <- P / 2 ## Initially, do not use this adaptive implementation
  
  ## Try a toy example with TWO child nodes
  TT_c <- dim(data)[1] ## Assume the data are stored in TWO columns in the example
  n_n <- dim(data)[2]  ## The nodes have the same number of observations
  
  # Storage
  x_c <- array(rep(NA, P*nt), dim = c(P, nt))
  x_c_nt <- matrix(rep(NA, P*n_n), ncol = n_n)
  x <- array(rep(NA, P*nt), dim = c(P, nt))
  #x_1_mP <- x_2_mP <- rep(NA, m*P)
  #x_mP <- cbind(x_1_mP, x_2_mP)
  W_c_nt <- matrix(rep(NA, P*n_n), ncol = n_n)
  
  w_log <- array(rep(NA, P*nt), dim = c(P, nt))
  w_t_log <- array(rep(NA, P*nt), dim = c(P, nt))
  
  pc_log <- matrix(rep(NA, P*n_n), ncol = n_n) # Store the marginal likelihood for v_t
  
  # Initialization
  x_c_0 <- rep(0, P)
  w_log_0 <- 0
  L_c_prod_log <- 0
  
  
  # Child nodes
  for (c_i in 1:n_n) {
    for (i in 1:(nt-1)) {
      if (i == 1) {
        x_c[,i] <- x_c_0 + rnorm(P, 0, 1e4)
        for (j in 1:P) {
          w_log[j,i] <- w_log_0 + sum(sapply(data[,c_i], function(y) dnorm(y, x_c[j,i], 1, log = T)))
          if (w_log[j,i] < log(.Machine$double.xmin)) w_log[j,i] <- log(.Machine$double.xmin)
        }
        W <- exp(w_log[,i]) / sum(exp(w_log[,i]))
        
        ESS <- 1/sum(W^2)
        if (ESS < ESS_min) {
          U <- runif(1, 0, 1)
          A <- Sys_resamp(W=W, P=P, U=U)
          x_c[,i] <- x_c[A,i]
        }
        
        L_c_prod_log <- L_c_prod_log + log(sum(exp(w_log[,i]))) - log(P)
      }
      
      x_c[,i+1] <- x_c[,i] + rnorm(P, 0, 1e4)
      for (j in 1:P) {
        w_log[j,i+1] <- sum(sapply(data[,c_i], function(y) dnorm(y, x_c[j,i+1], 1, log = T)))
        if (w_log[j,i+1] < log(.Machine$double.xmin)) w_log[j,i+1] <- log(.Machine$double.xmin)
      }
      W <- exp(w_log[,i+1]) / sum(exp(w_log[,i+1]))
      
      ### The last iteration is not resampled --- for the following mP times resampling
      if (i != (nt-1)) {
        ESS <- 1/sum(W^2)
        if (ESS < ESS_min) {
          U <- runif(1, 0, 1)
          A <- Sys_resamp(W=W, P=P, U=U)
          x_c[,i+1] <- x_c[A,i+1]
        }
      }
      
      L_c_prod_log <- L_c_prod_log + log(sum(exp(w_log[,i+1]))) - log(P)
    }
    
    x_c_nt[,c_i] <- x_c[,nt]
    W_c_nt[,c_i] <- W
    pc_log[,c_i] <- w_log[,nt]
  }
  
  
  # mP matchings storage
  x_c_mP <- matrix(rep(NA, m*P*n_n), ncol = n_n)
  pc_log_mP <- matrix(rep(NA, m*P*n_n), ncol = n_n) 
  ## Resampling
  for (c_i in 1:n_n) {
    U <- runif(1, 0, 1)
    A <- Sys_resamp(W=W_c_nt[,c_i], P=m*P, U=U)
    x_c_mP[,c_i] <- x_c_nt[A,c_i]
    pc_log_mP[,c_i] <- pc_log[A,c_i]
  }
  
  
  
  # Resampling P tuples with weights V_t
  v_t_log <- rep(NA, m*P)
  p_check_log <- rep(NA, m*P) # Store pi_check
  for (j in 1:(m*P)) {
    pc_sum_log <- sum(pc_log_mP[j,])
    p_int_log <- sum(sapply(x_c_mP[j,], function(y) dnorm(y, sum(x_c_mP[j,])/(1+n_n), 1, log = T)))
    v_t_log[j] <- p_int_log - pc_sum_log
    v_t_log[j] <- v_t_log[j]*alpha
    if (v_t_log[j] < log(.Machine$double.xmin)) v_t_log[j] <- log(.Machine$double.xmin)
    p_check_log[j] <- (1-alpha)*pc_sum_log + alpha*p_int_log
  }
  if (sum(exp(v_t_log)) == Inf) v_t_log <- v_t_log - max(v_t_log)
  V_t <- exp(v_t_log)/sum(exp(v_t_log))
  U <- runif(1, 0, 1)
  A <- Sys_resamp(W=V_t, P=P, U=U)
  x_c_re <- x_c_mP[A,]
  p_check_log_re <- p_check_log[A]
  
  L_prod_log <- L_c_prod_log + log(sum(exp(v_t_log))) - log(m*P)
  
  
  
  # Root
  
  ## alpha increment for the annealing implementation
  alpha_inc <- 1/nt
  
  ## Initialization
  x_c_sum <- rowSums(x_c_re)
  x_0 <- x_c_sum/(1+n_n) + sqrt(1/(1+n_n)) * rnorm(P, 0, 1)
  w_t_log_0 <- 0
  
  ## Sampler iterations
  for (i in 1:(nt-1)) {
    if (i == 1) {
      
      ### Update weights
      for (j in 1:P) {
        gamma_update <- sum(sapply(x_c_re[j,], function(y) dnorm(y, x_0[j], 1, log = T)))
        - p_check_log_re[j] - dnorm(x_0[j], x_c_sum[j], sqrt(1/(1+n_n)), log = T)
        gamma_update <- gamma_update*alpha_inc
        w_t_log[j,i] <- w_t_log_0 + gamma_update
        #w_t_log[j,i] <- w_t_log_0 + sum(sapply(x_c_re[j,], function(y) dnorm(y, x[j,i], 1, log = T))) ## Should it be x_c_nt[j,] or x_c_nt?
        if (w_t_log[j,i] < log(.Machine$double.xmin)) w_t_log[j,i] <- log(.Machine$double.xmin)
      }
      W <- exp(w_t_log[,i]) / sum(exp(w_t_log[,i]))
      
      ### Update the marginal likelihood
      L_prod_log <- L_prod_log + log(sum(exp(w_t_log[,i]))) - log(P)
      
      ### Resampling
      ESS <- 1/sum(W^2)
      if (ESS < ESS_min) {
        U <- runif(1, 0, 1)
        A <- Sys_resamp(W=W, P=P, U=U)
        x_0 <- x_0[A]  
      }
      
      
      ### Proposal
      x[,i] <- x_0 + rnorm(P, 0, 1e4)
    }
    
    ### Update weights
    for (j in 1:P) {
      gamma_update <- sum(sapply(x_c_re[j,], function(y) dnorm(y, x[j,i], 1, log = T)))
      - p_check_log_re[j] - dnorm(x[j,i], x_c_sum[j], sqrt(1/(1+n_n)), log = T)
      gamma_update <- gamma_update*alpha_inc
      w_t_log[j,i+1] <- gamma_update
      #w_t_log[j,i+1] <- sum(sapply(x_c_re[j,], function(y) dnorm(y, x[j,i+1], 1, log = T))) ## Same question as above
      if (w_t_log[j,i+1] < log(.Machine$double.xmin)) w_t_log[j,i+1] <- log(.Machine$double.xmin)
    }
    W <- exp(w_t_log[,i+1]) / sum(exp(w_t_log[,i+1]))
    
    
    ### Update the marginal likelihood
    L_prod_log <- L_prod_log + log(sum(exp(w_t_log[,i+1]))) - log(P)
    
    ### Resampling
    ESS <- 1/sum(W^2)
    if (ESS < ESS_min) {
      U <- runif(1, 0, 1)
      A <- Sys_resamp(W=W, P=P, U=U)
      x[,i] <- x[A,i]  
    }
    
    
    ### Proposal
    x[,i+1] <- x[,i] + rnorm(P, 0, 1e4)
    
  }
  
  x_t <- x[,nt]
  x_t_n <- x[,nt-1]
  
  
  return(list(x_t=x_t, x_t_n=x_t_n, x_c_nt=x_c_nt, x_t_tilde=x_0, W=W))
}


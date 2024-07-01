dc_smc_leaf_binom <- function(data, alpha, beta, n_trial, N) {
  min_lim <- log(.Machine$double.xmin)
  n <- length(data)
  
  # Storage
  p_bin <- rep(NA, N)
  p_bin_re <- rep(NA, N)
  
  w_log <- rep(NA, N)
  W <- rep(NA, N)
  
  theta1 <- alpha + sum(data)
  theta2 <- beta + n_trial*n - sum(data)
  
  p_bin <- rbeta(N, theta1, theta2)
  
  ## Update weights
  w_log <- weight_update_binom(data=data,
                               p=p_bin,
                               n_trial=n_trial,
                               min_lim=min_lim)
  if (max(w_log) == min_lim) {
    W <- rep(1/N, N)
  } else {
    W <- exp(w_log) / sum(exp(w_log))
  }
  
  L_log <- sum(w_log) - log(N)
  
  return(list(p_bin_post=p_bin,
              W=W,
              L_log=L_log))
}














dc_smc <- function() {
  ESS_min <- P / 2 ## Initially, not use this adaptive implementation
  
  ## Try a toy example with TWO child nodes
  TT_c <- dim(data)[1] ## Assume the data are stored in TWO columns in the example
  n_n <- dim(data)[2]  ## The nodes have the same number of observations
  
  # Storage
  x_1 <- x_2 <- array(rep(NA, Ntotal*P), dim = c(P, Ntotal))
  x <- list(x1, x2)
  w_log <- rep(NA, P)
  
  # Initialization
  x_0 <- rep(0, n_n)
  w_log_0 <- 0
  
  for (c_i in 1:n_n) {
    L_current <- 0
    x_ci <- x[[c_i]]
    x_ci[1,i] <- x_0[c_i] + rnorm(1, 0, 1)
    w_log[1,i] <- w_log_0 + sum(apply(data[,c_i], 1, function(y) dnorm(y, x_ci[1,i], 1, log = T)))
    for (j in 2:P) {
      x_ci[j,i] <- x_ci[j-1,i] + rnorm(1, 0, 1)
      
    }
  }
}

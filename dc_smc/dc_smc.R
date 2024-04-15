dc_smc <- function() {
  ESS_min <- P / 2 ## Initially, not use this adaptive implementation

  ## Try a toy example with TWO child nodes
  TT_c <- dim(data)[1] ## Assume the data are stored in TWO columns in the example
  n_n <- dim(data)[2]  ## The nodes have the same number of observations
  
  # Storage
  x_1 <- x_2 <- array(rep(NA, Ntotal*TT_c*P), dim = c(P, TT_c, Ntotal))
  
  # Initialization
  x_0 <- matrix(rep(0, P*n_n), ncol = n_n)
  
  for (c_i in 1:2) {
    L_current <- 0
    for (t in 1:TT_c) {
      
    } 
    
  }
}
sigma <- .1
rho <- .95
mu <- -1.5
x <- rep(NA, 100)
x[1] <- rnorm(1, 0, sigma)
for (i in 2:100) {
  x[i] <- rho*(x[i-1]-mu) + mu + sigma*rnorm(1, 0, 1)
}
data <- sqrt(exp(x))*rnorm(100, 0, 1)



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


pf <- function(data, P) {
  TT <- length(data)
  
  # Storage
  x <- matrix(rep(NA, P*TT), ncol = TT)
  w_log <- rep(NA, P)
  w_t_log <- rep(NA, P)
  #W <- rep(NA, P)
  W_t <- rep(NA, P)
  mu_t <- matrix(rep(NA, P*TT), ncol = TT)
  mu_t_star <- matrix(rep(NA, P*TT), ncol = TT)
  mu_t_re <- matrix(rep(NA, P*TT), ncol = TT)
  mu_t_star_re <- matrix(rep(NA, P*TT), ncol = TT) 
  g_t_log <- rep(NA, P) ## store p(y_t|mu_t)
  
  # Fixed parameter values for the auxiliary part
  sigma_0 <- .1
  rho_0 <- .9
  mu_0 <- -1
  
  # Initialization
  x_0 <- rep(0, P)
  w_log_0 <- 0
  L_log <- 0 # mariginal likelihood (or normalization constant)
  
  for (t in 1:(TT-1)) {
    if (t == 1) {
      mu_t[,t] <- mu_0 + rho_0*(x_0-mu_0)
      mu_t_star[,t] <- mu_t[,t] + sigma_0^2/4*(data[t]^2*exp(-mu_t[,t])-2)
      for (j in 1:P) {
        g_t_log[j] <- (mu_t_star[j,t]^2 - mu_t[j,t]^2)/2/sigma_0^2 - data[t]^2/2*exp(-mu_t[j,t])*(1+mu_t[j,t])
        w_t_log[j] <- w_log_0 + g_t_log[j]
        if (w_t_log[j] < log(.Machine$double.xmin)) w_t_log[j] <- log(.Machine$double.xmin)
        x[j,t] <- mu_t_star[j,t] + rnorm(1, 0, sigma_0)
        w_log[j] <- w_log_0 + dnorm(data[t], 0, exp(x[j,t]/2), log=T) - g_t_log[j]
        if (w_log[j] < log(.Machine$double.xmin)) w_log[j] <- log(.Machine$double.xmin)
      }
      W_t <- exp(w_t_log) / sum(exp(w_t_log))
      #W <- exp(w_log) / sum(exp(w_log))
      
      ## Update the marginal likelihood
      L_log <- L_log + log(sum(exp(w_log))) - log(P)
    }
    
    mu_t[,t+1] <- mu_0 + rho_0*(x[,t]-mu_0)
    mu_t_star[,t+1] <- mu_t[,t+1] + sigma_0^2/4*(data[t+1]^2*exp(-mu_t[,t+1])-2)
    for (j in 1:P) {
      g_t_log[j] <- (mu_t_star[j,t+1]^2 - mu_t[j,t+1]^2)/2/sigma_0^2 - data[t+1]^2/2*exp(-mu_t[j,t+1])*(1+mu_t[j,t+1])
      w_t_log[j] <- w_log[j] + g_t_log[j]
      #if (w_t_log[j] < log(.Machine$double.xmin)) w_t_log[j] <- log(.Machine$double.xmin)
      if (w_t_log[j] > log(.Machine$double.xmax)) w_t_log[j] <- log(.Machine$double.xmax)
    }
    w_t_log <- w_t_log - max(w_t_log)
    W_t <- exp(w_t_log) / sum(exp(w_t_log))
    
    # Resampling with W_t
    U <- runif(1, 0, 1)
    A <- Sys_resamp(W=W_t, P=P, U=U)
    mu_t_re[,t+1] <- mu_0 + rho_0*(x[A,t]-mu_0)
    mu_t_star_re[,t+1] <- mu_t_re[,t+1] + sigma_0^2/4*(data[t+1]^2*exp(-mu_t_re[,t+1])-2)
    for (j in 1:P) {
      x[j,t+1] <- mu_t_star_re[j,t+1] + rnorm(1, 0, sigma_0)
      ## Update w
      w_log[j] <- dnorm(data[t+1], 0, exp(x[j,t+1]/2), log=T) - g_t_log[j]
      if (w_log[j] < log(.Machine$double.xmin)) w_log[j] <- log(.Machine$double.xmin)
    }
    #W <- exp(w_log) / sum(exp(w_log))
    
    ## Update the marginal likelihood
    L_log <- L_log + log(sum(exp(w_log))) - log(P)
  }
  
  return(list(x=x, MLik=exp(L_log), mu_t_star=mu_t_star, Wt=W_t))
}










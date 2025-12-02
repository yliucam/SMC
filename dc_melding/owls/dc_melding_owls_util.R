pooling_prior_log_log <- function(lambda,
                                  alpha,
                                  rho,
                                  mu,
                                  sigma,
                                  trun,
                                  a,
                                  b) {
  
  Z <- pnorm(trun[2], mu[1], sigma[1]) - pnorm(trun[1], mu[1], sigma[1])
  
  res <- (lambda[1] + lambda[2] - 1) * trunc_normal_uni_log(data = alpha[,1],
                                                            mu = mu,
                                                            sigma = sigma,
                                                            Z = Z)
  
  res <- res + (lambda[1] + lambda[2] - 1) * trunc_normal_uni_log(data = alpha[,2],
                                                                  mu = mu,
                                                                  sigma = sigma,
                                                                  Z = Z)
  
  res <- res + (lambda[2] + lambda[3] - 1) * unif_log(data = rho,
                                                      a = a,
                                                      b = b) 
  
  return(res)
}


# Generate (vectorized) random numbers from truncated Poisson distribution
# using the inverse-CDF sampling
rpois_trunc <- function(lambda, K) {
  # vector of lambda allowed
  F_K <- ppois(K, lambda)                 # P(X ≤ K)
  u <- runif(length(lambda)) * F_K        # uniform on [0, F_K]
  qpois(u, lambda)                        # inverse CDF Poisson
}



# Generate (vectorized) random numbers from truncated Binomial distribution
# using the inverse-CDF sampling
rbinom_trunc <- function(n, p, K) {
  F_K <- pbinom(K, n, p)           # P(X ≤ K)
  u <- runif(length(p)) * F_K
  qbinom(u, n, p)
}







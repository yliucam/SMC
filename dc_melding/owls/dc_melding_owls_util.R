pooling_prior_log_log <- function(lambda,
                                  alpha,
                                  rho,
                                  mu,
                                  sigma,
                                  trun,
                                  a,
                                  b) {
  
  Z <- pnorm(trun[2], mu[1], sigma[1]) - pnorm(trun[1], mu[1], sigma[1])
  
  res <- lambda[1] * lambda[2] * trunc_normal_uni_log(data = alpha[,1],
                                                      mu = mu,
                                                      sigma = sigma,
                                                      Z = Z)
  
  res <- res + lambda[1] * lambda[2] * trunc_normal_uni_log(data = alpha[,2],
                                                            mu = mu,
                                                            sigma = sigma,
                                                            Z = Z)
  
  res <- res + lambda[2] * lambda[3] * unif_log(data = rho,
                                                a = a,
                                                b = b)
  
  return(res)
}





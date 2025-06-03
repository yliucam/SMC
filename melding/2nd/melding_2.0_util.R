library(mvtnorm)

pooling_prior_log_log <- function(lambda,
                                  phi12,
                                  phi23) {
  res <- lambda[1] * dnorm(phi12, 0, 1e3, log = T) +
    lambda[3] * dnorm(phi23, 0, 1e3, log = T) +
    lambda[2] * dmvnorm(c(phi12, phi23), rep(0, 2), cbind(c(1e6, .8*1e6), c(.8*1e6, 1e6)), log = T)
  
  return(res)
}

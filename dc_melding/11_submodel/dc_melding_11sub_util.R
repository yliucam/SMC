# Pooled prior p_{pool,2} for submodel 2
pooled_2_prior_log <- function(phi12,
                               phi23,
                               mu,
                               sigma,
                               alpha,
                               beta,
                               lambda1,
                               lambda2,
                               lambda3) {
  res <- (lambda1 + lambda2 - 1) * normal_uni_log(phi12, mu, sigma)
  res <- res + (lambda2 + lambda3 - 1) * gamma_log(phi23, alpha, beta)
  
  return(res)
}



# Pooled prior p_{pool,10} for submodel 10
pooled_10_prior_log <- function(phi9_10,
                                phi10_11,
                                mu,
                                sigma,
                                alpha,
                                beta,
                                lambda1,
                                lambda2,
                                lambda3) {
  res <- (lambda1 + lambda2 - 1) * normal_uni_log(phi9_10, mu, sigma)
  res <- res + (lambda2 + lambda3 - 1) * gamma_log(phi10_11, alpha, beta)
  
  return(res)
}



# Pooled prior p_{pool,4} for submodel 4
pooled_4_prior_log <- function(phi34,
                               phi45,
                               mu,
                               sigma,
                               alpha,
                               beta,
                               lambda1,
                               lambda2,
                               lambda3) {
  res <- (lambda1 + lambda2 - 1) * normal_uni_log(phi34, mu, sigma)
  res <- res + (lambda2 + lambda3 - 1) * gamma_log(phi45, alpha, beta)
  
  return(res)
}



# Pooled prior p_{pool,8} for submodel 8
pooled_8_prior_log <- function(phi78,
                               phi89,
                               mu,
                               sigma,
                               alpha,
                               beta,
                               lambda1,
                               lambda2,
                               lambda3) {
  res <- (lambda1 + lambda2 - 1) * normal_uni_log(phi78, mu, sigma)
  res <- res + (lambda2 + lambda3 - 1) * gamma_log(phi89, alpha, beta)
  
  return(res)
}



# Pooled prior p_{pool,8} for submodel 8
pooled_6_prior_log <- function(phi56,
                               phi67,
                               mu,
                               sigma,
                               alpha,
                               beta,
                               lambda1,
                               lambda2,
                               lambda3) {
  res <- (lambda1 + lambda2 - 1) * gamma_log(phi56, alpha, beta) 
  res <- (lambda2 + lambda3 - 1) * normal_uni_log(phi67, mu, sigma)
  
  return(res)
}


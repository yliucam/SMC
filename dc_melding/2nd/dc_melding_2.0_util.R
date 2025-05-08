pooling_prior_log_log <- function(lambda,
                                  phi,
                                  mu_uni,
                                  sigma_uni,
                                  mu_mult,
                                  Sigma_mult) {
  res <- lambda[1] * normal_uni_log(phi[,1], mu_uni[,1], sigma_uni[,1]) +
    lambda[3] * normal_uni_log(phi[,2], mu_uni[,2], sigma_uni[,2]) +
    lambda[2] * normal_mult_log(phi, mu_mult, Sigma_mult)
  
  return(res)
}
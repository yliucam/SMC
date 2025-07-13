N <- 8000
m <- 100
n_sub <- 3

para_sub <- cbind(t(out_recap$alpha[50,c(1,3),]), out_fec$rho[,50])
lambda <- c(1, 1, 1)

para_sub_mN <- array(rep(NA, m*N*n_sub), dim = c(m*N, n_sub))

for (sub_i in 1:n_sub) {
  A <- sample(1:N, m*N, replace = T)
  para_sub_mN[,sub_i] <- para_sub[A, sub_i]
}

u_pooling_log_mN <- pooling_prior_log_log(lambda = lambda,
                                          alpha = para_sub_mN[,1:2],
                                          rho = para_sub_mN[,3],
                                          mu = 0,
                                          sigma = 2,
                                          trun = c(-10, 10),
                                          a = 0,
                                          b = 10)
p_approx_log <- merging_loglik(data = as.matrix(count),
                               alpha = para_sub_mN[,1:2],
                               rho = para_sub_mN[,3])
V_log <- .02 * (u_pooling_log_mN + p_approx_log)
min_lim <- log(.Machine$double.xmin)
if (quantile(V_log, .8) < min_lim) V_log <- V_log / 1e6
V_log[which(V_log < min_lim)] <- min_lim
V <- exp(V_log)
V <- V / sum(V)

A <- sample(1:(m*N), N, replace = T, prob = V)
para_sub_merged <- para_sub_mN[A,]


u_pooling_log <- u_pooling_log_mN[A]
p_approx_log <- p_approx_log[A]









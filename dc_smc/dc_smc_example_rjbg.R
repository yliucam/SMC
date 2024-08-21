set.seed(24035)

alpha <- rlnorm(1, 3, 2)

alpha1 <- rgamma(1, alpha, 1)
alpha2 <- rgamma(1, alpha, .5)

p11 <- rbeta(1, alpha1, 2)
p12 <- rbeta(1, alpha1, 2)
p13 <- rbeta(1, alpha1, 2)
p14 <- rbeta(1, alpha1, 2)

p21 <- rbeta(1, alpha2, 10)
p22 <- rbeta(1, alpha2, 10)
p23 <- rbeta(1, alpha2, 10)
p24 <- rbeta(1, alpha2, 10)

y11 <- rbinom(20, 100, p11)
y12 <- rbinom(20, 100, p12)
y13 <- rbinom(20, 100, p13)
y14 <- rbinom(20, 100, p14)

y21 <- rbinom(20, 200, p21)
y22 <- rbinom(20, 200, p22)
y23 <- rbinom(20, 200, p23)
y24 <- rbinom(20, 200, p24)

data1 <- cbind(y11, y12, y13, y14, y21, y22, y23, y24)


library(Rcpp)
library(MASS)
sourceCpp("dc_smc_util.cpp")
source("dc_smc_util.R")
source("dc_smc_leaf_binom.R")

source("dc_smc_algB2_new.R")

beta_prior <- list(alpha=rep(1, 8), beta=c(rep(2, 4), rep(10, 4)))
n_trial <- c(rep(100, 4), rep(200, 4))
gamma_prior <- list(alpha = c(1, 1), beta = c(1, .5))
LN_prior <- list(mu = 3, sigma = 2)

out_smc <- dc_smc_algB2_new(data = data1, 
                            n_trial = n_trial,
                            N = 1000,
                            nodes_n = c(1, 2, 8),
                            beta_prior = beta_prior,
                            gamma_prior = gamma_prior,
                            LN_prior = LN_prior,
                            m = 2,
                            alpha = .05,
                            Ntotal_sub = 5,
                            Ntotal_root = 5)


library(tidyverse)
library(ggplot2)

library(posterior)

source("example_BUGS.R")


# draws

library(ggdist)

# create draws object for DC-SMC results
leaf_node_names <- paste0("p[", rep(1:2, each = 4), ",", 1:4, "]")
x_leaf_draws <- out_smc$x_leaf
colnames(x_leaf_draws) <- leaf_node_names
x_leaf_draws <- as_draws_matrix(x_leaf_draws)

subroot_node_names <- paste0("alpha[", 1:2, "]")
x_subroot_draws <- out_smc$x_sub
colnames(x_subroot_draws) <- subroot_node_names
x_subroot_draws <- as_draws_matrix(x_subroot_draws)

root_node_names <- paste0("alpha0")
x_root_draws <- matrix(out_smc$x_root, ncol = 1)
colnames(x_root_draws) <- root_node_names
x_root_draws <- as_draws(x_root_draws)

dc_draws <- as_draws_matrix(bind_cols(x_leaf_draws, x_subroot_draws, x_root_draws))

# leaf nodes
summarise_leaf <- function(draws){
draws |>
  spread_draws(p[i, j]) |>
  point_interval()
}

jags_spread <- jags_draws |>
  summarise_leaf() |>
  mutate(.sampler = "jags")

dc_spread <- dc_draws |>
  summarise_leaf() |>
  mutate(.sampler = "dc")

bind_rows(jags_spread, dc_spread) |>
  ggplot(aes(y = interaction(i, j), x = p, xmin = .lower, xmax = .upper, colour = .sampler)) +
    geom_pointinterval(position = position_dodge(width = 0.2))

# subroot
summarise_subroot <- function(draws){
  draws |>
    spread_draws(alpha[i]) |>
    point_interval()
}

jags_spread <- jags_draws |>
  summarise_subroot() |>
  mutate(.sampler = "jags")

dc_spread <- dc_draws |>
  summarise_subroot() |>
  mutate(.sampler = "dc")

bind_rows(jags_spread, dc_spread) |>
  ggplot(aes(y = i, x = alpha, xmin = .lower, xmax = .upper, colour = .sampler)) +
    geom_pointinterval(position = position_dodge(width = 0.2))

# root
summarise_root <- function(draws){
  draws |>
    spread_draws(alpha0) |>
    point_interval()
}

jags_spread <- jags_draws |>
  summarise_root() |>
  mutate(.sampler = "jags")

dc_spread <- dc_draws |>
  summarise_root() |>
  mutate(.sampler = "dc")

bind_rows(jags_spread, dc_spread) |>
  ggplot(aes(y = 0, x = alpha0, xmin = .lower, xmax = .upper, colour = .sampler)) +
    geom_pointinterval(position = position_dodge(width = 0.2))

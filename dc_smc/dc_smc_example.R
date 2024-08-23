set.seed(24001)

alpha <- rlnorm(1, 3, 2)

alpha1 <- rgamma(1, alpha, 1)
alpha2 <- rgamma(1, alpha, .5)

p1 <- rbeta(30, alpha1, 2)
p2 <- rbeta(30, alpha2, 10)

y1 <- matrix(rep(NA, 20*30), ncol = 30)
y2 <- matrix(rep(NA, 20*30), ncol = 30)
for (i in 1:30) {
  y1[,i] <- rbinom(20, 100, p1[i])
  y2[,i] <- rbinom(20, 200, p2[i])
}


data1 <- cbind(y1, y2)

library(Rcpp)
library(MASS)
sourceCpp("C:/Users/Yixuan/Documents/codes/SMC/dc_smc/dc_smc_util.cpp")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_smc/dc_smc_util.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_smc/dc_smc_leaf_binom.R")

source("C:/Users/Yixuan/Documents/codes/SMC/dc_smc/dc_smc_algB2_new.R")


beta_prior <- list(alpha=rep(1, 60), beta=c(rep(2, 30), rep(10, 30)))
n_trial <- c(rep(100, 30), rep(200, 30))
gamma_prior <- list(beta = c(1, .5))

debug(dc_smc_algB2_new)

out_smc <- dc_smc_algB2_new(data = data1, 
                        n_trial = n_trial,
                        N = 5000,
                        nodes_n = c(1, 2, 60),
                        beta_prior = beta_prior,
                        gamma_prior = gamma_prior,
                        m = 2,
                        alpha = .05,
                        Ntotal_sub = 5,
                        Ntotal_root = 5)


library(tidyverse)
library(ggplot2)

# Convert SMC array to a tibble for plotting
smc_array_as_tibble <- function(x){
  colnames(x) <- paste0("iteration", 1:ncol(x))
  rownames(x) <- 1:nrow(x)
  x_df <- as_tibble(x, rownames = "particle")
  x_df |>
    pivot_longer(cols = -particle) |>
    mutate(name = as.integer(str_remove(name, "iteration")))
}

# x values
x_root_df <- smc_array_as_tibble(out_smc$x_root_array)

ggplot(x_root_df |> filter(name > 2), aes(x = name, y = value, group = particle)) +
  geom_point(alpha = 0.05) +
  geom_line(alpha = 0.05)

ggplot(x_root_df, aes(x = value)) +
  geom_histogram(binwidth = 0.5) +
  facet_wrap("name", scales = "free_x")

# weights
W_root_df <- smc_array_as_tibble(out_smc$W_root_array)

W_root_df |>
  group_by(name) |>
  summarise(min = min(value),
            max = max(value),
            mean = mean(value))

ggplot(W_root_df |> filter(name > 2), aes(x = name, y = value, group = particle)) +
  geom_point(alpha = 0.05) +
  geom_line(alpha = 0.05)

ggplot(W_root_df, aes(x = value)) +
  geom_histogram(binwidth = 0.01) +
  facet_wrap("name", scales = "free_x")

source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/1st/dc_melding_1.1_mcmc.R")
source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/1st/dc_melding_1.1.R")

set.seed(12345)
#psi3 <- rnorm(1, 1, 100)
psi3 <- rgamma(1, 8, 8)
phi12 <- rnorm(1, psi3, 1)
phi23 <- rnorm(1, psi3, 2)

psi1 <- rgamma(1, 2, 4)
psi2 <- rgamma(1, 9, 6)

y1 <- rnorm(50, phi12, psi1)
y2 <- rnorm(50, phi23, psi2)

debug(dc_melding_1.1)
out <- dc_melding_1.1(data=cbind(y1, y2), 
                      N=2000, 
                      n_leaf=2, 
                      m=10, 
                      alpha=.05, 
                      Ntotal=50)

# Convert dc_melding output to format for ggplot
library(tidyverse)
colnames(out$x) <- paste0("root_", 1:20)
root_df <- as_tibble(out$x) |>
  pivot_longer(cols = everything(),
               names_to = c("node", "step"),
               names_sep = "_") |>
  mutate(step = factor(step, levels = 1:20))

colnames(out$x_leaf) <- c("phi12", "phi23")
leaf_df1 <- as_tibble(out$x_leaf) |>
  pivot_longer(cols = everything()) |>
  mutate(type = "marginal")

colnames(out$x_leaf_merge) <- c("phi12", "phi23")
leaf_df2 <- as_tibble(out$x_leaf_merge) |>
  pivot_longer(cols = everything()) |>
  mutate(type = "postVresample")

leaf_df <- bind_rows(leaf_df1, leaf_df2)

# Pretty histogram plots
library(colorspace)
geom_histogram_outline <- function(alpha = 0.3, ...){
  list(geom_histogram(colour = NA,
                      alpha = alpha,
                      position = "identity",
                      ...),
       geom_step(aes(y = after_stat(count),
                     colour = after_scale(colorspace::darken(fill, 0.1))),
                 stat = "bin",
                 pad = TRUE,
                 direction = "mid",
                 ...))
}

# Compare original samples (x_leaf) and post merge samples (x_leaf_merge)
# for phi12 and phi23
ggplot(leaf_df, aes(x = value, fill = type)) +
  facet_grid(rows = vars(name)) +
  geom_histogram_outline(origin = 0, binwidth = 0.1)

# Assess evolution of the phi3 samples (x) for the root through the tempering
# SMC steps
ggplot(root_df, aes(x = value, fill = step)) +
  facet_grid(rows = "step") +
  geom_histogram_outline(alpha = 0.4, binwidth = 0.2)


library(rjags)
model <- " model {
  psi3 ~ dgamma(1, 1)
  
  phi12 ~ dnorm(psi3, 1)
  phi23 ~ dnorm(psi3, 1/4)
  
  psi1 ~ dgamma(1, 1)
  psi2 ~ dgamma(1, 1)
  
  for (i in 1:50) {
    y1[i] ~ dnorm(phi12, 1/psi1^2)
    y2[i] ~ dnorm(phi23, 1/psi2^2)
  }
}
"

parameters <- c("phi12", "phi23", "psi3")

burn_in <- 30000

steps <- 100000

thin <- 5

foo <- jags.model(textConnection(model), data = list(y1=y1, y2=y2))

update(foo, burn_in)

jags <- coda.samples(model=foo, variable.names = parameters, n.iter = steps, thin = thin)

out_bugs <- as.matrix(jags)


# Leaf submodel for phi12
model_leaf1 <- " model {
  phi12 ~ dnorm(0, 1)
  psi1 ~ dgamma(1, 1)
  for (i in 1:50) {
    y1[i] ~ dnorm(phi12, 1/psi1^2)
  }
}
"

parameters <- c("phi12")
burn_in <- 30000
steps <- 100000
thin <- 5
foo_leaf1 <- jags.model(textConnection(model_leaf1), data = list(y1=y1, y2=y2))
update(foo_leaf1, burn_in)
jags_leaf1 <- coda.samples(model=foo_leaf1, variable.names = parameters, n.iter = steps, thin = thin)


# Leaf submodel for phi23
model_leaf2 <- " model {
  phi23 ~ dnorm(0,1)
  psi2 ~ dgamma(1, 1)
  for (i in 1:50) {
    y2[i] ~ dnorm(phi23, 1/psi2^2)
  }
}
"
parameters <- c("phi23")
burn_in <- 30000
steps <- 100000
thin <- 5
foo_leaf2 <- jags.model(textConnection(model_leaf2), data = list(y1=y1, y2=y2))
update(foo_leaf2, burn_in)

jags_leaf2 <- coda.samples(model=foo_leaf2, variable.names = parameters, n.iter = steps, thin = thin)

# compare the various distributions for JAGS
summary(jags)
summary(jags_leaf1)
summary(jags_leaf2)

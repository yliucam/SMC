
# Examining weights etc

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

# leaf
x_leaf_df <- smc_array_as_tibble(out_smc$x_leaf)

x_leaf_df |>
  group_by(name) |>
  summarise(min = min(value),
            max = max(value),
            mean = mean(value))

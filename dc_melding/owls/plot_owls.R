library(ggplot2)
library(gridExtra)
library(grid)
library(dplyr)
library(patchwork)
library(tidyr)



original_model_samples <- readRDS("C:/Users/Yixuan/Documents/codes/SMC/melding/owls/result/original-ipm-samples.rds")
sub1_samples <- readRDS("C:/Users/Yixuan/Documents/codes/SMC/melding/owls/result/capture-recapture-subposterior-samples.rds")
sub3_samples <- readRDS("C:/Users/Yixuan/Documents/codes/SMC/melding/owls/result/fecundity-subposterior-samples.rds")
sub2_samples <- readRDS("C:/Users/Yixuan/Documents/codes/SMC/melding/owls/result/count-data-subposterior-samples.rds")
stage_two_samples <- readRDS("C:/Users/Yixuan/Documents/codes/SMC/melding/owls/result/melded-posterior-samples.rds")
vars <- c("fec", "v[1]", "v[2]", "v[6]")
orig_samples_slim <- array(
  original_model_samples[, , vars],
  dim = c(
    dim(original_model_samples)[1] * dim(original_model_samples)[2],
    4
  ), 
  dimnames = list(
    NULL,
    vars
  )
)

meld_samples_slim <- array(
  stage_two_samples[, , vars],
  dim = c(
    dim(stage_two_samples)[1] * dim(stage_two_samples)[2],
    4
  ), 
  dimnames = list(
    NULL,
    vars
  )
)

sub2_samples_slim <- array(
  sub2_samples[,,vars],
  dim = c(
    dim(sub2_samples)[1] * dim(sub2_samples)[2],
    4
  ),
  dimnames = list(
    NULL,
    vars
  )
)

sub1_samples_slim <- array(
  sub1_samples[,,vars[2:3]],
  dim = c(
    dim(sub1_samples)[1] * dim(sub1_samples)[2],
    2
  ),
  dimnames = list(
    NULL,
    vars[2:3]
  )
)

sub3_samples_slim <- array(
  sub3_samples,
  dim = c(
    dim(sub3_samples)[1] * dim(sub3_samples)[2],
    1
  ),
  dimnames = list(
    NULL,
    "fec"
  )
)


# load("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/results/submodel2_mcmc.RData")
# submodel2_mcmc <- out_jags

load("C:/Users/Yixuan/Documents/codes/SMC/melding/owls/result/pointwise_results.RData")
pointwise_sample <- out_jags2

load("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/results/out_owls_dc_melding_N_16000_alpha_0.5.RData")
dc_melding_sample <- out

load("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/results/recap_result_N_16000.RData")
load("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/results/fecundity_result_N_16000.RData")



#alpha0_dcm_stage_one <- out_recap$alpha[50,1,]
alpha0_sub1 <- sub1_samples_slim[,"v[1]"]
alpha0_sub2 <- sub2_samples_slim[,"v[1]"]
alpha0_dcm <- dc_melding_sample$phi[,1,2]
#alpha0_m_stage_one <- meld_recap_samples_slim[,"v[1]"]
alpha0_m <- meld_samples_slim[,"v[1]"]
alpha0_ipm <- orig_samples_slim[,"v[1]"]
# alpha0_sub2_mcmc <- submodel2_mcmc[,"alpha0"]
# alpha0 <- cbind(alpha0_dcm_stage_one, alpha0_dcm, alpha0_m, alpha0_ipm, alpha0_sub2_mcmc)
# alpha0 <- data.frame(alpha0)
# colnames(alpha0) <- c("dc_melding stage one", "dc_melding", "melding", "ipm", "submodel 2")

#alpha2_dcm_stage_one <- out_recap$alpha[50,3,]
alpha2_sub1 <- sub1_samples_slim[,"v[2]"]
alpha2_sub2 <- sub2_samples_slim[,"v[2]"]
alpha2_dcm <- dc_melding_sample$phi[,2,2]
#alpha2_m_stage_one <- meld_recap_samples_slim[,"v[2]"]
alpha2_m <- meld_samples_slim[,"v[2]"]
alpha2_ipm <- orig_samples_slim[,"v[2]"]
# alpha2_sub2_mcmc <- submodel2_mcmc[,"alpha2"]
# alpha2 <- cbind(alpha2_dcm_stage_one, alpha2_dcm, alpha2_m, alpha2_ipm, alpha2_sub2_mcmc)
# alpha2 <- data.frame(alpha2)
# colnames(alpha2) <- c("dc_melding stage one", "dc_melding", "melding", "ipm", "submodel 2")

#rho_dcm_stage_one <- out_fec$rho[,11]
rho_sub3 <- sub3_samples_slim[,"fec"]
rho_sub2 <- sub2_samples_slim[,"fec"]
rho_dcm <- dc_melding_sample$phi[,3,2]
#rho_m_stage_one <- meld_fec_samples_slim[,"fec"]
rho_m <- meld_samples_slim[,"fec"]
rho_ipm <- orig_samples_slim[,"fec"]
# rho_sub2_mcmc <- submodel2_mcmc[,"rho"]
# rho <- cbind(rho_dcm_stage_one, rho_dcm, rho_m, rho_ipm, rho_sub2_mcmc)
# rho <- data.frame(rho)
# colnames(rho) <- c("dc_melding stage one", "dc_melding", "melding", "ipm", "submodel 2")

alpha6_sub2 <- sub2_samples_slim[,"v[6]"]
alpha6_dcm <- dc_melding_sample$psi2$alpha6[,3]
alpha6_m <- meld_samples_slim[,"v[6]"]
alpha6_ipm <- orig_samples_slim[,"v[6]"]
# alpha6_sub2_mcmc <- submodel2_mcmc[,"alpha6"]
alpha6_pointwise <- pointwise_sample[,"alpha6"]


get_intervals <- function(samples, method_name) {
  probs <- c(0.005, 0.025, 0.10, 0.25, 0.75, 0.90, 0.975, 0.995)
  q <- quantile(samples, probs)
  
  data.frame(
    method = method_name,
    level = c("99%", "95%", "80%", "50%"),
    xmin = c(q[1], q[2], q[3], q[4]),
    xmax = c(q[8], q[7], q[6], q[5]),
    y = 1  # Dummy y position
  )
}


intervals_alpha0_dcm_stage_one <- get_intervals(out_recap$alpha[50,1,], "dc_melding stage one")
intervals_alpha0_dcm <- get_intervals(dc_melding_sample$phi[,1,51], "dc_melding")
intervals_alpha0_m <- get_intervals(meld_samples_slim[,"v[1]"], "melding")
intervals_alpha0_ipm <- get_intervals(orig_samples_slim[,"v[1]"], "ipm")
intervals_alpha0 <- rbind(intervals_alpha0_dcm_stage_one, intervals_alpha0_dcm, intervals_alpha0_m, intervals_alpha0_ipm)
intervals_alpha0$y <- ifelse(intervals_alpha0$method == "dc_melding stage one", 4, 
                             ifelse(intervals_alpha0$method == "dc_melding", 3,
                                    ifelse(intervals_alpha0$method == "melding", 2, 1)))
intervals_alpha0$alpha <- recode(intervals_alpha0$level,
                                     "99%" = 0.2,
                                     "95%" = 0.4,
                                     "80%" = 0.6,
                                     "50%" = 0.8)

alpha0_plot <- ggplot(intervals_alpha0) +
  geom_rect(aes(xmin = xmin, xmax = xmax,
                ymin = y - 0.3, ymax = y + 0.3,
                fill = method, alpha = alpha)) +
  # geom_segment(#data = data.frame(y = c(1, 2)),
  #   aes(x = ps2, xend = ps2,
  #       y = y - 0.3, yend = y + 0.3),
  #   linewidth = 1) +
  scale_y_continuous(breaks = c(1, 2), labels = NULL) +
  scale_alpha_identity() +
  labs(x = expression(alpha[0]), y = NULL) +
  #theme_minimal() +
  theme(axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        legend.position = "none", panel.grid.major.y = element_blank())



intervals_alpha2_dcm_stage_one <- get_intervals(out_recap$alpha[50,3,], "dc_melding stage one")
intervals_alpha2_dcm <- get_intervals(dc_melding_sample$phi[,2,6], "dc_melding")
intervals_alpha2_m <- get_intervals(meld_samples_slim[,"v[2]"], "melding")
intervals_alpha2_ipm <- get_intervals(orig_samples_slim[,"v[2]"], "ipm")
intervals_alpha2 <- rbind(intervals_alpha2_dcm_stage_one, intervals_alpha2_dcm, intervals_alpha2_m, intervals_alpha2_ipm)
intervals_alpha2$y <- ifelse(intervals_alpha2$method == "dc_melding stage one", 4, 
                             ifelse(intervals_alpha2$method == "dc_melding", 3,
                                    ifelse(intervals_alpha2$method == "melding", 2, 1)))
intervals_alpha2$alpha <- recode(intervals_alpha2$level,
                                 "99%" = 0.2,
                                 "95%" = 0.4,
                                 "80%" = 0.6,
                                 "50%" = 0.8)

alpha2_plot <- ggplot(intervals_alpha2) +
  geom_rect(aes(xmin = xmin, xmax = xmax,
                ymin = y - 0.3, ymax = y + 0.3,
                fill = method, alpha = alpha)) +
  # geom_segment(#data = data.frame(y = c(1, 2)),
  #   aes(x = ps2, xend = ps2,
  #       y = y - 0.3, yend = y + 0.3),
  #   linewidth = 1) +
  scale_y_continuous(breaks = c(1, 2), labels = NULL) +
  scale_alpha_identity() +
  labs(x = expression(alpha[2]), y = NULL) +
  #theme_minimal() +
  theme(axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        legend.position = "none", panel.grid.major.y = element_blank())




intervals_rho_dcm_stage_one <- get_intervals(out_fec$rho[,50], "dc_melding stage one")
intervals_rho_dcm <- get_intervals(dc_melding_sample$phi[,3,6], "dc_melding")
intervals_rho_m <- get_intervals(meld_samples_slim[,"fec"], "melding")
intervals_rho_ipm <- get_intervals(orig_samples_slim[,"fec"], "ipm")
intervals_rho <- rbind(intervals_rho_dcm_stage_one, intervals_rho_dcm, intervals_rho_m, intervals_rho_ipm)
intervals_rho$y <- ifelse(intervals_rho$method == "dc_melding stage one", 4, 
                          ifelse(intervals_rho$method == "dc_melding", 3,
                                 ifelse(intervals_rho$method == "melding", 2, 1)))
intervals_rho$alpha <- recode(intervals_rho$level,
                                 "99%" = 0.2,
                                 "95%" = 0.4,
                                 "80%" = 0.6,
                                 "50%" = 0.8)

rho_plot <- ggplot(intervals_rho) +
  geom_rect(aes(xmin = xmin, xmax = xmax,
                ymin = y - 0.3, ymax = y + 0.3,
                fill = method, alpha = alpha)) +
  # geom_segment(#data = data.frame(y = c(1, 2)),
  #   aes(x = ps2, xend = ps2,
  #       y = y - 0.3, yend = y + 0.3),
  #   linewidth = 1) +
  scale_y_continuous(breaks = c(1, 2), labels = NULL) +
  scale_alpha_identity() +
  labs(x = expression(rho), y = NULL) +
  #theme_minimal() +
  theme(axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.spacing.y = unit(.5, "cm"),
        panel.grid.major.y = element_blank())


pdf(file = "C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/results/owls_credible_intervals_alpha_0.02.pdf", width = 18, height = 12)

(alpha0_plot | alpha2_plot | rho_plot) +
  plot_annotation(title = "50%, 80%, 95% and 99% posterior credible intervals",
                  theme = theme(plot.title = element_text(size = 30, hjust = .5,
                                                          face = "bold")))

dev.off()



#---------------------------box-plot------------------------------------#

##-----------------------individual_submodels---------------------------##
alpha0_sub_long <- data.frame(
  value = c(alpha0_sub1, alpha0_sub2, alpha0_ipm),
  method = factor(rep(
    c("p1", "p2", "ipm"),
    times = c(length(alpha0_sub1),
              length(alpha0_sub2),
              length(alpha0_ipm))
  ))
)


alpha2_sub_long <- data.frame(
  value = c(alpha2_sub1, alpha0_sub2, alpha2_ipm),
  method = factor(rep(
    c("p1", "p2", "ipm"),
    times = c(length(alpha0_sub1),
              length(alpha0_sub2),
              length(alpha0_ipm))
  ))
)


rho_sub_long <- data.frame(
  value = c(rho_sub2, rho_sub3, rho_ipm),
  method = factor(rep(
    c("p2", "p3", "ipm"),
    times = c(length(rho_sub2),
              length(rho_sub3),
              length(rho_ipm))
  ))
)


# Common set of levels for ALL plots (even if some are absent in a given plot)
fill_levels1 <- c("p1","p2","ipm")
fill_levels2 <- c("p2","p3","ipm")

fill_scale_common1 <- scale_fill_manual(
  values = c("p1"="olivedrab3", "p2"="lightpink", "ipm"="grey10"),
  labels = c("p1"=expression(p[1]),
             "p2"=expression(p[2]),
             "ipm"=expression(p[ipm])),
  name = "Model"
)

fill_scale_common2 <- scale_fill_manual(
  values = c("p2"="lightpink", "p3"="cyan3", "ipm"="grey10"),
  labels = c("p2"=expression(p[2]),
             "p3"=expression(p[3]),
             "ipm"=expression(p[ipm])),
  name = NULL
)

# Ensure factor levels are consistent
alpha0_sub_long$method <- factor(alpha0_sub_long$method, levels = fill_levels1)
alpha2_sub_long$method <- factor(alpha2_sub_long$method, levels = fill_levels1)
rho_sub_long$method <- factor(rho_sub_long$method, levels = fill_levels2)

# --- Keep legend ONLY on the first plot ---
alpha0_sub_boxplot <- ggplot(alpha0_sub_long, aes(x = method, y = value, fill = method)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  coord_flip() +
  scale_x_discrete(labels = NULL) +
  fill_scale_common1 +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = NULL, y = expression(alpha[0])) +
  theme(axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 40),
        panel.grid.major.y = element_blank(),
        legend.position = "none")

# Hide legend on the others
alpha2_sub_boxplot <- ggplot(alpha2_sub_long, aes(x = method, y = value, fill = method)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  coord_flip() +
  scale_x_discrete(labels = NULL) +
  fill_scale_common1 +
  labs(x = NULL, y = expression(alpha[2])) +
  theme(axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 40),
        panel.grid.major.y = element_blank(),
        legend.position = "none")

rho_sub_boxplot <- ggplot(rho_sub_long, aes(x = method, y = value, fill = method)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  coord_flip() +
  scale_x_discrete(labels = NULL) +
  fill_scale_common2 +
  labs(x = NULL, y = expression(rho)) +
  theme(axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 40),
        panel.grid.major.y = element_blank(),
        legend.position = "none")


fill_levels <- c("p1","p2","p3","ipm")

legend_df <- data.frame(
  method = factor(fill_levels, levels = fill_levels),
  x = 1, y = 1
)

legend_plot <- ggplot(legend_df, aes(x, y, fill = method)) +
  geom_boxplot() +
  scale_fill_manual(
    limits = fill_levels,
    drop = FALSE,
    values = c("p1"="olivedrab3", "p2"="lightpink", "p3"="cyan3", "ipm"="grey50"),
    labels = c("p1"  = expression(p[1]),
               "p2"  = expression(p[2]),
               "p3"  = expression(p[3]),
               "ipm" = expression(p[ipm])),
    name = "Model "
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 40),
    legend.text  = element_text(size = 40),
    # kill any remaining grey background/panel
    panel.background = element_blank(),
    plot.background  = element_blank(),
  )


merged_plot <- legend_plot /
  (alpha0_sub_boxplot | alpha2_sub_boxplot | rho_sub_boxplot) +
  plot_layout(heights = c(.08, 1))

ggsave(filename = "C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/results/owls_boxplot_stage_one_alpha_0.5.jpg",
       plot = merged_plot,
       width = 18,
       height = 12,
       units = "in",
       dpi = 300)





##-----------------------stage_two---------------------------##
alpha0_long <- data.frame(
  value = c(alpha0_dcm, alpha0_m, alpha0_ipm),
  method = factor(rep(
    c("dc_melding", "melding", "ipm"),
    times = c(length(alpha0_dcm),
              length(alpha0_m),
              length(alpha0_ipm))
  ))
)

alpha2_long <- data.frame(
  value = c(alpha2_dcm, alpha2_m, alpha2_ipm),
  method = factor(rep(
    c("dc_melding", "melding", "ipm"),
    times = c(length(alpha2_dcm),
              length(alpha2_m),
              length(alpha2_ipm))
  ))
)

rho_long <- data.frame(
  value = c(rho_dcm, rho_m, rho_ipm),
  method = factor(rep(
    c("dc_melding", "melding", "ipm"),
    times = c(length(rho_dcm),
              length(rho_m),
              length(rho_ipm))
  ))
)


alpha0_boxplot <- ggplot(alpha0_long, aes(x = method, y = value, fill = method)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  coord_flip() +   # <-- makes boxes horizontal
  scale_fill_manual(
    values = c("dc_melding" = "olivedrab3",
               "melding"    = "cyan3",
               "ipm"        = "grey30")
  ) +
  labs(x = NULL, y = expression(alpha[0]), fill = 'Method') +
  scale_x_discrete(labels = NULL) +
  #theme_minimal(base_size = 16) +
  theme(axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 40),
        legend.position = "none", panel.grid.major.y = element_blank())

alpha2_boxplot <- ggplot(alpha2_long, aes(x = method, y = value, fill = method)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  coord_flip() +   # <-- makes boxes horizontal
  scale_fill_manual(
    values = c("dc_melding" = "olivedrab3",
               "melding"    = "cyan3",
               "ipm"        = "grey30")
  ) +
  labs(x = NULL, y = expression(alpha[2]), fill = 'Method') +
  scale_x_discrete(labels = NULL) +
  #theme_minimal(base_size = 16) +
  theme(axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 40),
        legend.position = "none", panel.grid.major.y = element_blank())

rho_boxplot <- ggplot(rho_long, aes(x = method, y = value, fill = method)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  coord_flip() +   # <-- makes boxes horizontal
  scale_fill_manual(
    values = c("dc_melding" = "olivedrab3",
               "melding"    = "cyan3",
               "ipm"        = "grey30")
  ) +
  labs(x = NULL, y = expression(rho), fill = 'Method') +
  scale_x_discrete(labels = NULL) +
  #theme_minimal(base_size = 16) +
  theme(axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 40),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 40),
        legend.spacing.y = unit(.5, "cm"),
        panel.grid.major.y = element_blank())


pdf(file = "C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/results/owls_boxplot_stage_two_alpha_0.5.pdf", width = 18, height = 12)

(alpha0_boxplot | alpha2_boxplot | rho_boxplot)
 
dev.off()





##--------------------------alpha_6------------------------------##
alpha6_sub_long <- data.frame(
  value = c(alpha6_sub2, alpha6_ipm),
  method = factor(rep(
    c("p2", "ipm"),
    times = c(length(alpha6_sub2),
              length(alpha6_ipm))
  ))
)


alpha6_sub_boxplot <- ggplot(alpha6_sub_long, aes(x = method, y = value, fill = method)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  # coord_flip() +   # <-- makes boxes horizontal
  scale_fill_manual(
    values = c("p2" = "lightpink",
               "ipm" = "grey30"),
    labels = c("p2"=expression(p[2]),
               "ipm"=expression(p[ipm]))
  ) +
  labs(x = NULL, y = expression(alpha[6]), fill = 'Model') +
  scale_y_continuous(position = "right") +
  scale_x_discrete(labels = NULL) +
  #theme_minimal(base_size = 16) +
  theme(axis.text = element_text(size = 20),
        axis.title.y.right = element_text(size = 40,
                                          angle = 0,
                                          vjust = 0.5,
                                          hjust = 0.5,
                                          margin = margin(l = 20)),
        legend.text = element_text(size = 50),
        legend.title = element_text(size = 50),
        legend.spacing.y = unit(.5, "cm"),
        legend.position = "left",
        panel.grid.major.y = element_blank())


pdf(file = "C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/results/owls_boxplot_alpha6_sub_alpha_0.5.pdf", width = 18, height = 12)

alpha6_sub_boxplot

dev.off()






alpha6_long <- data.frame(
  value = c(alpha6_dcm, alpha6_m, alpha6_ipm, alpha6_pointwise),
  method = factor(rep(
    c("dc_melding", "melding", "ipm", "pointwise"),
    times = c(length(alpha6_dcm),
              length(alpha6_m),
              length(alpha6_ipm),
              length(alpha6_pointwise))
  ))
)


alpha6_boxplot <- ggplot(alpha6_long, aes(x = method, y = value, fill = method)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  scale_fill_manual(
    values = c("dc_melding" = "olivedrab3",
               "melding"    = "cyan3",
               "ipm"        = "grey30",
               "pointwise"  = "lightpink")
  ) +
  labs(x = NULL, y = expression(alpha[6]), fill = 'Method') +
  scale_x_discrete(labels = NULL) +
  #theme_minimal(base_size = 16) +
  theme(axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 40, 
                                    angle = 0,
                                    vjust = 0.5,
                                    margin = margin(r = 15)),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 40),
        legend.spacing.y = unit(.5, "cm"),
        panel.grid.major.y = element_blank())

pdf(file = "C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/results/owls_boxplot_alpha6_alpha_0.5.pdf", width = 18, height = 12)

alpha6_boxplot

dev.off()


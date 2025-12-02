library(ggplot2)
library(gridExtra)
library(grid)
library(dplyr)
library(patchwork)
library(tidyr)



original_model_samples <- readRDS("C:/Users/Yixuan/Documents/codes/SMC/melding/owls/result/original-ipm-samples.rds")
stage_one_recap_samples <- readRDS("C:/Users/Yixuan/Documents/codes/SMC/melding/owls/result/capture-recapture-subposterior-samples.rds")
stage_one_fec_samples <- readRDS("C:/Users/Yixuan/Documents/codes/SMC/melding/owls/result/fecundity-subposterior-samples.rds")
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

meld_recap_samples_slim <- array(
  stage_one_recap_samples[,,c("v[1]", "v[2]")],
  dim = c(
    dim(stage_one_recap_samples)[1] * dim(stage_one_recap_samples)[2],
    2
  ),
  dimnames = list(
    NULL,
    c("v[1]", "v[2]")
  )
)

meld_fec_samples_slim <- array(
  stage_one_fec_samples[,,"rho"],
  dim = c(
    dim(stage_one_fec_samples)[1] * dim(stage_one_fec_samples)[2],
    1
  ),
  dimnames = list(
    NULL,
    "fec"
  )
)

load("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/results/submodel2_mcmc.RData")
submodel2_mcmc <- out_jags

load("C:/Users/Yixuan/Documents/codes/SMC/melding/owls/result/pointwise_results.RData")
pointwise_sample <- out_jags2

load("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/results/out_owls_dc_melding_N_16000_alpha_0.5.RData")
dc_melding_sample <- out

load("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/results/recap_result_N_16000.RData")
load("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/results/fecundity_result_N_16000.RData")



alpha0_dcm_stage_one <- out_recap$alpha[50,1,]
alpha0_dcm <- dc_melding_sample$phi[,1,2]
alpha0_m_stage_one <- meld_recap_samples_slim[,"v[1]"]
alpha0_m <- meld_samples_slim[,"v[1]"]
alpha0_ipm <- orig_samples_slim[,"v[1]"]
alpha0_sub2_mcmc <- submodel2_mcmc[,"alpha0"]
alpha0 <- cbind(alpha0_dcm_stage_one, alpha0_dcm, alpha0_m, alpha0_ipm, alpha0_sub2_mcmc)
alpha0 <- data.frame(alpha0)
colnames(alpha0) <- c("dc_melding stage one", "dc_melding", "melding", "ipm", "submodel 2")

alpha2_dcm_stage_one <- out_recap$alpha[50,3,]
alpha2_dcm <- dc_melding_sample$phi[,2,2]
alpha2_m_stage_one <- meld_recap_samples_slim[,"v[2]"]
alpha2_m <- meld_samples_slim[,"v[2]"]
alpha2_ipm <- orig_samples_slim[,"v[2]"]
alpha2_sub2_mcmc <- submodel2_mcmc[,"alpha2"]
alpha2 <- cbind(alpha2_dcm_stage_one, alpha2_dcm, alpha2_m, alpha2_ipm, alpha2_sub2_mcmc)
alpha2 <- data.frame(alpha2)
colnames(alpha2) <- c("dc_melding stage one", "dc_melding", "melding", "ipm", "submodel 2")

rho_dcm_stage_one <- out_fec$rho[,11]
rho_dcm <- dc_melding_sample$phi[,3,2]
rho_m_stage_one <- meld_fec_samples_slim[,"fec"]
rho_m <- meld_samples_slim[,"fec"]
rho_ipm <- orig_samples_slim[,"fec"]
rho_sub2_mcmc <- submodel2_mcmc[,"rho"]
rho <- cbind(rho_dcm_stage_one, rho_dcm, rho_m, rho_ipm, rho_sub2_mcmc)
rho <- data.frame(rho)
colnames(rho) <- c("dc_melding stage one", "dc_melding", "melding", "ipm", "submodel 2")

alpha6_dcm <- dc_melding_sample$psi2$alpha6[,3]
alpha6_m <- meld_samples_slim[,"v[6]"]
alpha6_ipm <- orig_samples_slim[,"v[6]"]
alpha6_sub2_mcmc <- submodel2_mcmc[,"alpha6"]
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

##-----------------------stage_one---------------------------##
alpha0_stage_one_long <- data.frame(
  value = c(alpha0_dcm_stage_one, alpha0_m_stage_one, alpha0_ipm, alpha0_sub2_mcmc),
  method = factor(rep(
    c("dc_melding", "melding", "ipm", "sub 2 only"),
    times = c(length(alpha0_dcm_stage_one),
              length(alpha0_m_stage_one),
              length(alpha0_ipm),
              length(alpha0_sub2_mcmc))
  ))
)


alpha2_stage_one_long <- data.frame(
  value = c(alpha2_dcm_stage_one, alpha2_m_stage_one, alpha2_ipm, alpha2_sub2_mcmc),
  method = factor(rep(
    c("dc_melding", "melding", "ipm", "sub 2 only"),
    times = c(length(alpha2_dcm_stage_one),
              length(alpha2_m_stage_one),
              length(alpha2_ipm),
              length(alpha2_sub2_mcmc))
  ))
)


rho_stage_one_long <- data.frame(
  value = c(rho_dcm_stage_one, rho_m_stage_one, rho_ipm, rho_sub2_mcmc),
  method = factor(rep(
    c("dc_melding", "melding", "ipm", "sub 2 only"),
    times = c(length(rho_dcm_stage_one),
              length(rho_m_stage_one),
              length(rho_ipm),
              length(rho_sub2_mcmc))
  ))
)


alpha0_stage_one_boxplot <- ggplot(alpha0_stage_one_long, aes(x = method, y = value, fill = method)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  coord_flip() +   # <-- makes boxes horizontal
  labs(x = NULL, y = expression(alpha[0]), fill = 'Method') +
  scale_x_discrete(labels = NULL) +
  #theme_minimal(base_size = 16) +
  theme(axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 40),
        legend.position = "none", panel.grid.major.y = element_blank())

alpha2_stage_one_boxplot <- ggplot(alpha2_stage_one_long, aes(x = method, y = value, fill = method)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  coord_flip() +   # <-- makes boxes horizontal
  labs(x = NULL, y = expression(alpha[2]), fill = 'Method') +
  scale_x_discrete(labels = NULL) +
  #theme_minimal(base_size = 16) +
  theme(axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 40),
        legend.position = "none", panel.grid.major.y = element_blank())

rho_stage_one_boxplot <- ggplot(rho_stage_one_long, aes(x = method, y = value, fill = method)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  coord_flip() +   # <-- makes boxes horizontal
  labs(x = NULL, y = expression(rho), fill = 'Method') +
  scale_x_discrete(labels = NULL) +
  #theme_minimal(base_size = 16) +
  theme(axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 40),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 40),
        legend.spacing.y = unit(.5, "cm"),
        panel.grid.major.y = element_blank())



pdf(file = "C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/results/owls_boxplot_stage_one_alpha_0.5.pdf", width = 18, height = 12)

(alpha0_stage_one_boxplot | alpha2_stage_one_boxplot | rho_stage_one_boxplot) +
  plot_annotation(theme = theme(plot.title = element_text(size = 30, hjust = .5,
                                                          face = "bold")))

dev.off()




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
  labs(x = NULL, y = expression(alpha[0]), fill = 'Method') +
  scale_x_discrete(labels = NULL) +
  #theme_minimal(base_size = 16) +
  theme(axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 40),
        legend.position = "none", panel.grid.major.y = element_blank())

alpha2_boxplot <- ggplot(alpha2_long, aes(x = method, y = value, fill = method)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  coord_flip() +   # <-- makes boxes horizontal
  labs(x = NULL, y = expression(alpha[2]), fill = 'Method') +
  scale_x_discrete(labels = NULL) +
  #theme_minimal(base_size = 16) +
  theme(axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 40),
        legend.position = "none", panel.grid.major.y = element_blank())

rho_boxplot <- ggplot(rho_long, aes(x = method, y = value, fill = method)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  coord_flip() +   # <-- makes boxes horizontal
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

(alpha0_boxplot | alpha2_boxplot | rho_boxplot) +
  plot_annotation(theme = theme(plot.title = element_text(size = 30, hjust = .5,
                                                          face = "bold")))

dev.off()





##--------------------------alpha_6------------------------------##

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
  labs(x = NULL, y = expression(alpha[6]), fill = 'Method') +
  scale_x_discrete(labels = NULL) +
  #theme_minimal(base_size = 16) +
  theme(axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 40),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 40),
        legend.spacing.y = unit(.5, "cm"),
        panel.grid.major.y = element_blank())

pdf(file = "C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/results/owls_boxplot_alpha6_alpha_0.5.pdf", width = 18, height = 12)

alpha6_boxplot +
  plot_annotation(theme = theme(plot.title = element_text(size = 30, hjust = .5,
                                                          face = "bold")))


dev.off()


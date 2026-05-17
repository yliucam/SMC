library(ggplot2)
library(gridExtra)
library(grid)
library(dplyr)
library(patchwork)
library(tidyr)

#-----------------------------------phi56---------------------------------#

i <- 10

filename_sub5 <- paste0("RES_n50_submodel5_repN", i, ".RData")
fileaddress_sub5 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub5/",
                           filename_sub5)
load(fileaddress_sub5)
out_stage_one <- out$sigma2[,6]


filename_sub6 <- paste0("RES_n50_submodel6_repN", i, ".RData")
fileaddress_sub6 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6/",
                           filename_sub6)
load(fileaddress_sub6)
out_dcm <- out$phi[,1,21]


filename_sub6 <- paste0("RES_n50_full_mcmc_repN", i, ".RData")
fileaddress_sub6 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/full_mcmc/",
                           filename_sub6)
load(fileaddress_sub6)
out_mcmc <- out_jags[1:5000,"phi56"]


filename_sub6_smc2 <- paste0("RES_n50_submodel6_smc2_repN", i, ".RData")
fileaddress_sub6_smc2 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6_smc2/",
                                filename_sub6_smc2)
load(fileaddress_sub6_smc2)
out_dcm_smc2 <- out$sigma[,11]


true_value_phi56 <- data$phi56


out_phi56_long <- data.frame(
  value = c(out_stage_one, out_dcm, out_mcmc, out_dcm_smc2),
  method = factor(rep(
    c("stage one", "dc-melding", "mcmc", "dc-melding smc2"),
    times = c(length(out_stage_one),
              length(out_dcm),
              length(out_mcmc),
              length(out_dcm_smc2))
  ),
  levels = c("stage one", "dc-melding", "mcmc", "dc-melding smc2")) # custom the order of the variables
)



phi56_boxplot <- ggplot(out_phi56_long, aes(x = method, y = value, fill = method)) +
  geom_boxplot(width = .5, alpha = .7) +
  # add the true value
  geom_hline(
    aes(yintercept = true_value_phi56, colour = "true value"),
    linetype = "dashed", linewidth = 1) +
  scale_fill_manual(
    values = c("stage one"  = "lightpink",
               "dc-melding" = "olivedrab3",
               "mcmc" = "grey30",
               "dc-melding smc2" = "cyan3"),
    labels = c("stage one"  = "stage one",
               "dc-melding" = "dc-melding",
               "mcmc" = "mcmc",
               "dc-melding smc2" = expression("dc-melding " * smc^2))
  ) +
  # Color scale JUST for the true-value line
  scale_colour_manual(
    values = c("true value" = "red"),
    name = NULL
  ) +
  labs(x = NULL, y = NULL, title = expression(phi[5 * "," * 6]), fill = 'method') +
  scale_x_discrete(labels = NULL) +
  theme(axis.text = element_text(size = 20),
        plot.title = element_text(size = 40, hjust = .5),
        legend.position = "none", panel.grid.major.y = element_blank())
  # theme(axis.text = element_text(size = 20),
  #       axis.title.y = element_text(size = 40),
  #       plot.title = element_text(size = 40, hjust = .5),
  #       legend.text = element_text(size = 40),
  #       legend.title = element_text(size = 40),
  #       legend.spacing.y = unit(.5, "cm"),
  #       panel.grid.major.y = element_blank())


pdf(file = "C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/boxplot_phi56.pdf",
    width = 18, height = 12)

phi56_boxplot

dev.off()







#-----------------------------------phi67---------------------------------#

i <- 10

filename_sub7 <- paste0("RES_n50_submodel7_repN", i, ".RData")
fileaddress_sub7 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub7/",
                           filename_sub7)
load(fileaddress_sub7)
out_stage_one <- out$mu1[,6]


filename_sub6 <- paste0("RES_n50_submodel6_repN", i, ".RData")
fileaddress_sub6 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6/",
                           filename_sub6)
load(fileaddress_sub6)
out_dcm <- out$phi[,2,21]


filename_sub6 <- paste0("RES_n50_full_mcmc_repN", i, ".RData")
fileaddress_sub6 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/full_mcmc/",
                           filename_sub6)
load(fileaddress_sub6)
out_mcmc <- out_jags[1:5000,"phi67"]


filename_sub6_smc2 <- paste0("RES_n50_submodel6_smc2_repN", i, ".RData")
fileaddress_sub6_smc2 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6_smc2/",
                                filename_sub6_smc2)
load(fileaddress_sub6_smc2)
out_dcm_smc2 <- out$mu[,11]


true_value_phi67 <- data$phi67


out_phi67_long <- data.frame(
  value = c(out_stage_one, out_dcm, out_mcmc, out_dcm_smc2),
  method = factor(rep(
    c("stage one", "dc-melding", "mcmc", "dc-melding smc2"),
    times = c(length(out_stage_one),
              length(out_dcm),
              length(out_mcmc),
              length(out_dcm_smc2))
  ),
  levels = c("stage one", "dc-melding", "mcmc", "dc-melding smc2")) # custom the order of the variables
)



phi67_boxplot <- ggplot(out_phi67_long, aes(x = method, y = value, fill = method)) +
  geom_boxplot(width = .5, alpha = .7) +
  # add the true value
  geom_hline(
    aes(yintercept = true_value_phi67, colour = "true value"),
    linetype = "dashed", linewidth = 1) +
  scale_fill_manual(
    values = c("stage one"  = "lightpink",
               "dc-melding" = "olivedrab3",
               "mcmc" = "grey30",
               "dc-melding smc2" = "cyan3"),
    labels = c("stage one"  = "stage one",
               "dc-melding" = "dc-melding",
               "mcmc" = "mcmc",
               "dc-melding smc2" = expression("dc-melding " * smc^2))
  ) +
  # Color scale JUST for the true-value line
  scale_colour_manual(
    values = c("true value" = "red"),
    name = NULL
  ) +
  labs(x = NULL, y = NULL, title = expression(phi[6 * "," * 7]), fill = 'method') +
  scale_x_discrete(labels = NULL) +
  theme(axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 40),
        plot.title = element_text(size = 40, hjust = .5),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 40),
        legend.spacing.y = unit(.5, "cm"),
        panel.grid.major.y = element_blank())


pdf(file = "C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/boxplot_phi67.pdf",
    width = 18, height = 12)

phi67_boxplot

dev.off()


## boxplot for phi56 and phi67 together
pdf(file = "C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/boxplot_phi56_67.pdf",
    width = 18, height = 12)

(phi56_boxplot | phi67_boxplot)

dev.off()







#-----------------------------------phi45---------------------------------#

i <- 10

filename_sub5 <- paste0("RES_n50_submodel5_repN", i, ".RData")
fileaddress_sub5 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub5/",
                           filename_sub5)
load(fileaddress_sub5)
out_stage_one <- out$sigma1[,6]


filename_sub6 <- paste0("RES_n50_submodel6_repN", i, ".RData")
fileaddress_sub6 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6/",
                           filename_sub6)
filename_sub4 <- paste0("RES_n50_submodel4_repN", i, ".RData")
fileaddress_sub4 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub4/",
                           filename_sub4)
load(fileaddress_sub6)
A6 <- out$phi_index[,1,21]
load(fileaddress_sub4)
out_dcm <- out$phi[A6,2,6]


filename_sub6 <- paste0("RES_n50_full_mcmc_repN", i, ".RData")
fileaddress_sub6 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/full_mcmc/",
                           filename_sub6)
load(fileaddress_sub6)
out_mcmc <- out_jags[1:5000,"phi45"]


filename_sub6_smc2 <- paste0("RES_n50_submodel6_smc2_repN", i, ".RData")
fileaddress_sub6_smc2 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6_smc2/",
                                filename_sub6_smc2)
filename_sub4 <- paste0("RES_n50_submodel4_repN", i, ".RData")
fileaddress_sub4 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub4/",
                           filename_sub4)
load(fileaddress_sub6_smc2)
A6 <- out$A5[,11]
load(fileaddress_sub4)
out_dcm_smc2 <- out$phi[A6,2,6]


true_value_phi45 <- data$phi45


out_phi45_long <- data.frame(
  value = c(out_stage_one, out_dcm, out_mcmc, out_dcm_smc2),
  method = factor(rep(
    c("stage one", "dc-melding", "mcmc", "dc-melding smc2"),
    times = c(length(out_stage_one),
              length(out_dcm),
              length(out_mcmc),
              length(out_dcm_smc2))
  ),
  levels = c("stage one", "dc-melding", "mcmc", "dc-melding smc2")) # custom the order of the variables
)



phi45_boxplot <- ggplot(out_phi45_long, aes(x = method, y = value, fill = method)) +
  geom_boxplot(width = .5, alpha = .7) +
  # add the true value
  geom_hline(
    aes(yintercept = true_value_phi45, colour = "true value"),
    linetype = "dashed", linewidth = 1) +
  scale_fill_manual(
    values = c("stage one"  = "lightpink",
               "dc-melding" = "olivedrab3",
               "mcmc" = "grey30",
               "dc-melding smc2" = "cyan3"),
    labels = c("stage one"  = "stage one",
               "dc-melding" = "dc-melding",
               "mcmc" = "mcmc",
               "dc-melding smc2" = expression("dc-melding " * smc^2))
  ) +
  # Color scale JUST for the true-value line
  scale_colour_manual(
    values = c("true value" = "red"),
    name = NULL
  ) +
  labs(x = NULL, y = NULL, title = expression(phi[4*","*5]), fill = 'method') +
  scale_x_discrete(labels = NULL) +
  theme(axis.text = element_text(size = 20),
        plot.title = element_text(size = 40, hjust = .5),
        legend.position = "none", panel.grid.major.y = element_blank())
  # theme(axis.text = element_text(size = 20),
  #       axis.title.y = element_text(size = 40),
  #       plot.title = element_text(size = 40, hjust = .5),
  #       legend.text = element_text(size = 40),
  #       legend.title = element_text(size = 40),
  #       legend.spacing.y = unit(.5, "cm"),
  #       panel.grid.major.y = element_blank())



pdf(file = "C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/boxplot_phi45.pdf",
    width = 18, height = 12)

phi45_boxplot

dev.off()







#-----------------------------------phi78---------------------------------#

i <- 10

filename_sub7 <- paste0("RES_n50_submodel7_repN", i, ".RData")
fileaddress_sub7 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub7/",
                           filename_sub7)
load(fileaddress_sub7)
out_stage_one <- out$mu2[,6]


filename_sub6 <- paste0("RES_n50_submodel6_repN", i, ".RData")
fileaddress_sub6 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6/",
                           filename_sub6)
filename_sub8 <- paste0("RES_n50_submodel8_repN", i, ".RData")
fileaddress_sub8 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub8/",
                           filename_sub8)
load(fileaddress_sub6)
A6 <- out$phi_index[,2,21]
load(fileaddress_sub8)
out_dcm <- out$phi[A6,1,6]


filename_sub6 <- paste0("RES_n50_full_mcmc_repN", i, ".RData")
fileaddress_sub6 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/full_mcmc/",
                           filename_sub6)
load(fileaddress_sub6)
out_mcmc <- out_jags[1:5000,"phi78"]


filename_sub6_smc2 <- paste0("RES_n50_submodel6_smc2_repN", i, ".RData")
fileaddress_sub6_smc2 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6_smc2/",
                                filename_sub6_smc2)
filename_sub8 <- paste0("RES_n50_submodel8_repN", i, ".RData")
fileaddress_sub8 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub8/",
                           filename_sub8)
load(fileaddress_sub6_smc2)
A6 <- out$A7[,11]
load(fileaddress_sub8)
out_dcm_smc2 <- out$phi[A6,1,6]


true_value_phi78 <- data$phi78


out_phi78_long <- data.frame(
  value = c(out_stage_one, out_dcm, out_mcmc, out_dcm_smc2),
  method = factor(rep(
    c("stage one", "dc-melding", "mcmc", "dc-melding smc2"),
    times = c(length(out_stage_one),
              length(out_dcm),
              length(out_mcmc),
              length(out_dcm_smc2))
  ),
  levels = c("stage one", "dc-melding", "mcmc", "dc-melding smc2")) # custom the order of the variables
)



phi78_boxplot <- ggplot(out_phi78_long, aes(x = method, y = value, fill = method)) +
  geom_boxplot(width = .5, alpha = .7) +
  # add the true value
  geom_hline(
    aes(yintercept = true_value_phi78, colour = "true value"),
    linetype = "dashed", linewidth = 1) +
  scale_fill_manual(
    values = c("stage one"  = "lightpink",
               "dc-melding" = "olivedrab3",
               "mcmc" = "grey30",
               "dc-melding smc2" = "cyan3"),
    labels = c("stage one"  = "stage one",
               "dc-melding" = "dc-melding",
               "mcmc" = "mcmc",
               "dc-melding smc2" = expression("dc-melding " * smc^2))
  ) +
  # Color scale JUST for the true-value line
  scale_colour_manual(
    values = c("true value" = "red"),
    name = NULL
  ) +
  labs(x = NULL, y = NULL, title = expression(phi[7*","*8]), fill = 'method') +
  scale_x_discrete(labels = NULL) +
  theme(axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 40),
        plot.title = element_text(size = 40, hjust = .5),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 40),
        legend.spacing.y = unit(.5, "cm"),
        panel.grid.major.y = element_blank())



pdf(file = "C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/boxplot_phi78.pdf",
    width = 18, height = 12)

phi78_boxplot

dev.off()



# boxplot for phi45 and phi78 together

pdf(file = "C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/boxplot_phi45_78.pdf",
    width = 18, height = 12)

(phi45_boxplot | phi78_boxplot)

dev.off()








#-----------------------------------phi34---------------------------------#

i <- 10

filename_sub3 <- paste0("RES_n50_submodel3_repN", i, ".RData")
fileaddress_sub3 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub3/",
                           filename_sub3)
load(fileaddress_sub3)
out_stage_one <- out$mu[,6]


filename_sub6 <- paste0("RES_n50_submodel6_repN", i, ".RData")
fileaddress_sub6 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6/",
                           filename_sub6)
filename_sub4 <- paste0("RES_n50_submodel4_repN", i, ".RData")
fileaddress_sub4 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub4/",
                           filename_sub4)
load(fileaddress_sub6)
A6 <- out$phi_index[,1,21]
load(fileaddress_sub4)
out_dcm <- out$phi[A6,1,6]


filename_sub6 <- paste0("RES_n50_full_mcmc_repN", i, ".RData")
fileaddress_sub6 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/full_mcmc/",
                           filename_sub6)
load(fileaddress_sub6)
out_mcmc <- out_jags[1:5000,"phi34"]


filename_sub6_smc2 <- paste0("RES_n50_submodel6_smc2_repN", i, ".RData")
fileaddress_sub6_smc2 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6_smc2/",
                                filename_sub6_smc2)
filename_sub4 <- paste0("RES_n50_submodel4_repN", i, ".RData")
fileaddress_sub4 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub4/",
                           filename_sub4)
load(fileaddress_sub6_smc2)
A6 <- out$A5[,11]
load(fileaddress_sub4)
out_dcm_smc2 <- out$phi[A6,1,6]


true_value_phi34 <- data$phi34


out_phi34_long <- data.frame(
  value = c(out_stage_one, out_dcm, out_mcmc, out_dcm_smc2),
  method = factor(rep(
    c("stage one", "dc-melding", "mcmc", "dc-melding smc2"),
    times = c(length(out_stage_one),
              length(out_dcm),
              length(out_mcmc),
              length(out_dcm_smc2))
  ),
  levels = c("stage one", "dc-melding", "mcmc", "dc-melding smc2")) # custom the order of the variables
)



phi34_boxplot <- ggplot(out_phi34_long, aes(x = method, y = value, fill = method)) +
  geom_boxplot(width = .5, alpha = .7) +
  # add the true value
  geom_hline(
    aes(yintercept = true_value_phi34, colour = "true value"),
    linetype = "dashed", linewidth = 1) +
  scale_fill_manual(
    values = c("stage one"  = "lightpink",
               "dc-melding" = "olivedrab3",
               "mcmc" = "grey30",
               "dc-melding smc2" = "cyan3"),
    labels = c("stage one"  = "stage one",
               "dc-melding" = "dc-melding",
               "mcmc" = "mcmc",
               "dc-melding smc2" = expression("dc-melding " * smc^2))
  ) +
  # Color scale JUST for the true-value line
  scale_colour_manual(
    values = c("true value" = "red"),
    name = NULL
  ) +
  labs(x = NULL, y = NULL, title = expression(phi[3*","*4]), fill = 'method') +
  scale_x_discrete(labels = NULL) +
  theme(axis.text = element_text(size = 20),
        plot.title = element_text(size = 40, hjust = .5),
        legend.position = "none", panel.grid.major.y = element_blank())
  # theme(axis.text = element_text(size = 20),
  #       axis.title.y = element_text(size = 40),
  #       plot.title = element_text(size = 40, hjust = .5),
  #       legend.text = element_text(size = 40),
  #       legend.title = element_text(size = 40),
  #       legend.spacing.y = unit(.5, "cm"),
  #       panel.grid.major.y = element_blank())



pdf(file = "C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/boxplot_phi34.pdf",
    width = 18, height = 12)

phi34_boxplot

dev.off()







#-----------------------------------phi89---------------------------------#

i <- 10

filename_sub9 <- paste0("RES_n50_submodel9_repN", i, ".RData")
fileaddress_sub9 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub9/",
                           filename_sub9)
load(fileaddress_sub9)
out_stage_one <- out$tau[,6]


filename_sub6 <- paste0("RES_n50_submodel6_repN", i, ".RData")
fileaddress_sub6 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6/",
                           filename_sub6)
filename_sub8 <- paste0("RES_n50_submodel8_repN", i, ".RData")
fileaddress_sub8 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub8/",
                           filename_sub8)
load(fileaddress_sub6)
A6 <- out$phi_index[,2,21]
load(fileaddress_sub8)
out_dcm <- out$phi[A6,2,6]


filename_mcmc <- paste0("RES_n50_full_mcmc_repN", i, ".RData")
fileaddress_mcmc <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/full_mcmc/",
                           filename_mcmc)
load(fileaddress_mcmc)
out_mcmc <- out_jags[1:5000,"phi89"]


filename_sub6_smc2 <- paste0("RES_n50_submodel6_smc2_repN", i, ".RData")
fileaddress_sub6_smc2 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6_smc2/",
                                filename_sub6_smc2)
filename_sub8 <- paste0("RES_n50_submodel8_repN", i, ".RData")
fileaddress_sub8 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub8/",
                           filename_sub8)
load(fileaddress_sub6_smc2)
A6 <- out$A7[,11]
load(fileaddress_sub8)
out_dcm_smc2 <- out$phi[A6,2,6]


true_value_phi89 <- data$phi89


out_phi89_long <- data.frame(
  value = c(out_stage_one, out_dcm, out_mcmc, out_dcm_smc2),
  method = factor(rep(
    c("stage one", "dc-melding", "mcmc", "dc-melding smc2"),
    times = c(length(out_stage_one),
              length(out_dcm),
              length(out_mcmc),
              length(out_dcm_smc2))
  ),
  levels = c("stage one", "dc-melding", "mcmc", "dc-melding smc2")) # custom the order of the variables
)



phi89_boxplot <- ggplot(out_phi89_long, aes(x = method, y = value, fill = method)) +
  geom_boxplot(width = .5, alpha = .7) +
  # add the true value
  geom_hline(
    aes(yintercept = true_value_phi89, colour = "true value"),
    linetype = "dashed", linewidth = 1) +
  scale_fill_manual(
    values = c("stage one"  = "lightpink",
               "dc-melding" = "olivedrab3",
               "mcmc" = "grey30",
               "dc-melding smc2" = "cyan3"),
    labels = c("stage one"  = "stage one",
               "dc-melding" = "dc-melding",
               "mcmc" = "mcmc",
               "dc-melding smc2" = expression("dc-melding " * smc^2))
  ) +
  # Color scale JUST for the true-value line
  scale_colour_manual(
    values = c("true value" = "red"),
    name = NULL
  ) +
  labs(x = NULL, y = NULL, title = expression(phi[8*","*9]), fill = 'method') +
  scale_x_discrete(labels = NULL) +
  theme(axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 40),
        plot.title = element_text(size = 40, hjust = .5),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 40),
        legend.spacing.y = unit(.5, "cm"),
        panel.grid.major.y = element_blank())



pdf(file = "C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/boxplot_phi89.pdf",
    width = 18, height = 12)

phi89_boxplot

dev.off()



# boxplot for phi34 and phi89 together

pdf(file = "C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/boxplot_phi34_89.pdf",
    width = 18, height = 12)

(phi34_boxplot | phi89_boxplot)

dev.off()








#-----------------------------------phi23--------------------------------#

i <- 10

filename_sub3 <- paste0("RES_n50_submodel3_repN", i, ".RData")
fileaddress_sub3 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub3/",
                           filename_sub3)
load(fileaddress_sub3)
out_stage_one <- out$tau[,6]


filename_sub6 <- paste0("RES_n50_submodel6_repN", i, ".RData")
fileaddress_sub6 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6/",
                           filename_sub6)
filename_sub4 <- paste0("RES_n50_submodel4_repN", i, ".RData")
fileaddress_sub4 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub4/",
                           filename_sub4)
filename_sub2 <- paste0("RES_n50_submodel2_repN", i, ".RData")
fileaddress_sub2 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub2/",
                           filename_sub2)
load(fileaddress_sub6)
A6 <- out$phi_index[,1,21]
load(fileaddress_sub4)
A4 <- out$phi_index[A6,1,6]
load(fileaddress_sub2)
out_dcm <- out$phi[A4,2,6]


filename_mcmc <- paste0("RES_n50_full_mcmc_repN", i, ".RData")
fileaddress_mcmc <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/full_mcmc/",
                           filename_mcmc)
load(fileaddress_mcmc)
out_mcmc <- out_jags[1:5000,"phi23"]


filename_sub6_smc2 <- paste0("RES_n50_submodel6_smc2_repN", i, ".RData")
fileaddress_sub6_smc2 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6_smc2/",
                                filename_sub6_smc2)
filename_sub4 <- paste0("RES_n50_submodel4_repN", i, ".RData")
fileaddress_sub4 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub4/",
                           filename_sub4)
filename_sub2 <- paste0("RES_n50_submodel2_repN", i, ".RData")
fileaddress_sub2 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub2/",
                           filename_sub2)
load(fileaddress_sub6_smc2)
A6 <- out$A5[,11]
load(fileaddress_sub4)
A4 <- out$phi_index[A6,1,6]
load(fileaddress_sub2)
out_dcm_smc2 <- out$phi[A4,2,6]


true_value_phi23 <- data$phi23


out_phi23_long <- data.frame(
  value = c(out_stage_one, out_dcm, out_mcmc, out_dcm_smc2),
  method = factor(rep(
    c("stage one", "dc-melding", "mcmc", "dc-melding smc2"),
    times = c(length(out_stage_one),
              length(out_dcm),
              length(out_mcmc),
              length(out_dcm_smc2))
  ),
  levels = c("stage one", "dc-melding", "mcmc", "dc-melding smc2")) # custom the order of the variables
)



phi23_boxplot <- ggplot(out_phi23_long, aes(x = method, y = value, fill = method)) +
  geom_boxplot(width = .5, alpha = .7) +
  # add the true value
  geom_hline(
    aes(yintercept = true_value_phi23, colour = "true value"),
    linetype = "dashed", linewidth = 1) +
  scale_fill_manual(
    values = c("stage one"  = "lightpink",
               "dc-melding" = "olivedrab3",
               "mcmc" = "grey30",
               "dc-melding smc2" = "cyan3"),
    labels = c("stage one"  = "stage one",
               "dc-melding" = "dc-melding",
               "mcmc" = "mcmc",
               "dc-melding smc2" = expression("dc-melding " * smc^2))
  ) +
  # Color scale JUST for the true-value line
  scale_colour_manual(
    values = c("true value" = "red"),
    name = NULL
  ) +
  labs(x = NULL, y = NULL, title = expression(phi[2 * "," *3]), fill = 'method') +
  scale_x_discrete(labels = NULL) +
  theme(axis.text = element_text(size = 20),
        plot.title = element_text(size = 40, hjust = .5),
        legend.position = "none", panel.grid.major.y = element_blank())
  # theme(axis.text = element_text(size = 20),
  #       axis.title.y = element_text(size = 40),
  #       plot.title = element_text(size = 40, hjust = .5),
  #       legend.text = element_text(size = 40),
  #       legend.title = element_text(size = 40),
  #       legend.spacing.y = unit(.5, "cm"),
  #       panel.grid.major.y = element_blank())



pdf(file = "C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/boxplot_phi23.pdf",
    width = 18, height = 12)

phi23_boxplot

dev.off()







#-----------------------------------phi9_10---------------------------------#

i <- 300 #10

filename_sub9 <- paste0("RES_n50_submodel9_repN", i, ".RData")
fileaddress_sub9 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub9/",
                           filename_sub9)
load(fileaddress_sub9)
out_stage_one <- out$mu[,6]


filename_sub6 <- paste0("RES_n50_submodel6_repN", i, ".RData")
fileaddress_sub6 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6/",
                           filename_sub6)
filename_sub8 <- paste0("RES_n50_submodel8_repN", i, ".RData")
fileaddress_sub8 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub8/",
                           filename_sub8)
filename_sub10 <- paste0("RES_n50_submodel10_repN", i, ".RData")
fileaddress_sub10 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub10/",
                            filename_sub10)
load(fileaddress_sub6)
A6 <- out$phi_index[,2,21]
load(fileaddress_sub8)
A8 <- out$phi_index[A6,2,6]
load(fileaddress_sub10)
out_dcm <- out$phi[A8,1,6]


filename_mcmc <- paste0("RES_n50_full_mcmc_repN", i, ".RData")
fileaddress_mcmc <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/full_mcmc/",
                           filename_mcmc)
load(fileaddress_mcmc)
out_mcmc <- out_jags[1:5000,"phi9_10"]


filename_sub6_smc2 <- paste0("RES_n50_submodel6_smc2_repN", i, ".RData")
fileaddress_sub6_smc2 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6_smc2/",
                                filename_sub6_smc2)
filename_sub8 <- paste0("RES_n50_submodel8_repN", i, ".RData")
fileaddress_sub8 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub8/",
                           filename_sub8)
filename_sub10 <- paste0("RES_n50_submodel10_repN", i, ".RData")
fileaddress_sub10 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub10/",
                            filename_sub10)
load(fileaddress_sub6_smc2)
A6 <- out$A7[,11]
load(fileaddress_sub8)
A8 <- out$phi_index[A6,2,6]
load(fileaddress_sub10)
out_dcm_smc2 <- out$phi[A8,1,6]


true_value_phi9_10 <- data$phi9_10


out_phi9_10_long <- data.frame(
  value = c(out_stage_one, out_dcm, out_mcmc, out_dcm_smc2),
  method = factor(rep(
    c("stage one", "dc-melding", "mcmc", "dc-melding smc2"),
    times = c(length(out_stage_one),
              length(out_dcm),
              length(out_mcmc),
              length(out_dcm_smc2))
  ),
  levels = c("stage one", "dc-melding", "mcmc", "dc-melding smc2")) # custom the order of the variables
)



phi9_10_boxplot <- ggplot(out_phi9_10_long, aes(x = method, y = value, fill = method)) +
  geom_boxplot(width = .5, alpha = .7) +
  # add the true value
  geom_hline(
    aes(yintercept = true_value_phi9_10, colour = "true value"),
    linetype = "dashed", linewidth = 1) +
  scale_fill_manual(
    values = c("stage one"  = "lightpink",
               "dc-melding" = "olivedrab3",
               "mcmc" = "grey30",
               "dc-melding smc2" = "cyan3"),
    labels = c("stage one"  = "stage one",
               "dc-melding" = "dc-melding",
               "mcmc" = "mcmc",
               "dc-melding smc2" = expression("dc-melding " * smc^2))
  ) +
  # Color scale JUST for the true-value line
  scale_colour_manual(
    values = c("true value" = "red"),
    name = NULL
  ) +
  labs(x = NULL, y = NULL, title = expression(phi[9 * "," * 10]), fill = 'method') +
  scale_x_discrete(labels = NULL) +
  theme(axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 40),
        plot.title = element_text(size = 40, hjust = .5),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 40),
        legend.spacing.y = unit(.5, "cm"),
        panel.grid.major.y = element_blank())



pdf(file = "C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/boxplot_phi89.pdf",
    width = 18, height = 12)

phi9_10_boxplot

dev.off()



# boxplot for phi34 and phi89 together

pdf(file = "C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/boxplot_phi23_9_10.pdf",
    width = 18, height = 12)

(phi23_boxplot | phi9_10_boxplot)

dev.off()








#-----------------------------------phi12--------------------------------#

i <- 10

filename_sub1 <- paste0("RES_n50_submodel1_repN", i, ".RData")
fileaddress_sub1 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub1/",
                           filename_sub1)
load(fileaddress_sub1)
out_stage_one <- out$mu[,6]


filename_sub6 <- paste0("RES_n50_submodel6_repN", i, ".RData")
fileaddress_sub6 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6/",
                           filename_sub6)
filename_sub4 <- paste0("RES_n50_submodel4_repN", i, ".RData")
fileaddress_sub4 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub4/",
                           filename_sub4)
filename_sub2 <- paste0("RES_n50_submodel2_repN", i, ".RData")
fileaddress_sub2 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub2/",
                           filename_sub2)
load(fileaddress_sub6)
A6 <- out$phi_index[,1,21]
load(fileaddress_sub4)
A4 <- out$phi_index[A6,1,6]
load(fileaddress_sub2)
out_dcm <- out$phi[A4,1,6]


filename_mcmc <- paste0("RES_n50_full_mcmc_repN", i, ".RData")
fileaddress_mcmc <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/full_mcmc/",
                           filename_mcmc)
load(fileaddress_mcmc)
out_mcmc <- out_jags[1:5000,"phi12"]


filename_sub6_smc2 <- paste0("RES_n50_submodel6_smc2_repN", i, ".RData")
fileaddress_sub6_smc2 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6_smc2/",
                                filename_sub6_smc2)
filename_sub4 <- paste0("RES_n50_submodel4_repN", i, ".RData")
fileaddress_sub4 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub4/",
                           filename_sub4)
filename_sub2 <- paste0("RES_n50_submodel2_repN", i, ".RData")
fileaddress_sub2 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub2/",
                           filename_sub2)
load(fileaddress_sub6_smc2)
A6 <- out$A5[,11]
load(fileaddress_sub4)
A4 <- out$phi_index[A6,1,6]
load(fileaddress_sub2)
out_dcm_smc2 <- out$phi[A4,1,6]


true_value_phi12 <- data$phi12


out_phi12_long <- data.frame(
  value = c(out_stage_one, out_dcm, out_mcmc, out_dcm_smc2),
  method = factor(rep(
    c("stage one", "dc-melding", "mcmc", "dc-melding smc2"),
    times = c(length(out_stage_one),
              length(out_dcm),
              length(out_mcmc),
              length(out_dcm_smc2))
  ),
  levels = c("stage one", "dc-melding", "mcmc", "dc-melding smc2")) # custom the order of the variables
)



phi12_boxplot <- ggplot(out_phi12_long, aes(x = method, y = value, fill = method)) +
  geom_boxplot(width = .5, alpha = .7) +
  # add the true value
  geom_hline(
    aes(yintercept = true_value_phi12, colour = "true value"),
    linetype = "dashed", linewidth = 1) +
  scale_fill_manual(
    values = c("stage one"  = "lightpink",
               "dc-melding" = "olivedrab3",
               "mcmc" = "grey30",
               "dc-melding smc2" = "cyan3"),
    labels = c("stage one"  = "stage one",
               "dc-melding" = "dc-melding",
               "mcmc" = "mcmc",
               "dc-melding smc2" = expression("dc-melding " * smc^2))
  ) +
  # Color scale JUST for the true-value line
  scale_colour_manual(
    values = c("true value" = "red"),
    name = NULL
  ) +
  labs(x = NULL, y = NULL, title = expression(phi[1 * "," *2]), fill = 'method') +
  scale_x_discrete(labels = NULL) +
  theme(axis.text = element_text(size = 20),
        plot.title = element_text(size = 40, hjust = .5),
        legend.position = "none", panel.grid.major.y = element_blank())
  # theme(axis.text = element_text(size = 20),
  #       axis.title.y = element_text(size = 40),
  #       plot.title = element_text(size = 40, hjust = .5),
  #       legend.text = element_text(size = 40),
  #       legend.title = element_text(size = 40),
  #       legend.spacing.y = unit(.5, "cm"),
  #       panel.grid.major.y = element_blank())



pdf(file = "C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/boxplot_phi12.pdf",
    width = 18, height = 12)

phi12_boxplot

dev.off()







#-----------------------------------phi10_11---------------------------------#

i <- 10

filename_sub11 <- paste0("RES_n50_submodel11_repN", i, ".RData")
fileaddress_sub11 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub11/",
                            filename_sub11)
load(fileaddress_sub11)
out_stage_one <- out$sigma[,6]


filename_sub6 <- paste0("RES_n50_submodel6_repN", i, ".RData")
fileaddress_sub6 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6/",
                           filename_sub6)
filename_sub8 <- paste0("RES_n50_submodel8_repN", i, ".RData")
fileaddress_sub8 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub8/",
                           filename_sub8)
filename_sub10 <- paste0("RES_n50_submodel10_repN", i, ".RData")
fileaddress_sub10 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub10/",
                            filename_sub10)
load(fileaddress_sub6)
A6 <- out$phi_index[,2,21]
load(fileaddress_sub8)
A8 <- out$phi_index[A6,2,6]
load(fileaddress_sub10)
out_dcm <- out$phi[A8,2,6]


filename_mcmc <- paste0("RES_n50_full_mcmc_repN", i, ".RData")
fileaddress_mcmc <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/full_mcmc/",
                           filename_mcmc)
load(fileaddress_mcmc)
out_mcmc <- out_jags[1:5000,"phi10_11"]


filename_sub6_smc2 <- paste0("RES_n50_submodel6_smc2_repN", i, ".RData")
fileaddress_sub6_smc2 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub6_smc2/",
                                filename_sub6_smc2)
filename_sub8 <- paste0("RES_n50_submodel8_repN", i, ".RData")
fileaddress_sub8 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub8/",
                           filename_sub8)
filename_sub10 <- paste0("RES_n50_submodel10_repN", i, ".RData")
fileaddress_sub10 <- paste0("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/sub10/",
                            filename_sub10)
load(fileaddress_sub6_smc2)
A6 <- out$A7[,11]
load(fileaddress_sub8)
A8 <- out$phi_index[A6,2,6]
load(fileaddress_sub10)
out_dcm_smc2 <- out$phi[A8,2,6]


true_value_phi10_11 <- data$phi10_11


out_phi10_11_long <- data.frame(
  value = c(out_stage_one, out_dcm, out_mcmc, out_dcm_smc2),
  method = factor(rep(
    c("stage one", "dc-melding", "mcmc", "dc-melding smc2"),
    times = c(length(out_stage_one),
              length(out_dcm),
              length(out_mcmc),
              length(out_dcm_smc2))
  ),
  levels = c("stage one", "dc-melding", "mcmc", "dc-melding smc2")) # custom the order of the variables
)



phi10_11_boxplot <- ggplot(out_phi10_11_long, aes(x = method, y = value, fill = method)) +
  geom_boxplot(width = .5, alpha = .7) +
  # add the true value
  geom_hline(
    aes(yintercept = true_value_phi10_11, colour = "true value"),
    linetype = "dashed", linewidth = 1) +
  scale_fill_manual(
    values = c("stage one"  = "lightpink",
               "dc-melding" = "olivedrab3",
               "mcmc" = "grey30",
               "dc-melding smc2" = "cyan3"),
    labels = c("stage one"  = "stage one",
               "dc-melding" = "dc-melding",
               "mcmc" = "mcmc",
               "dc-melding smc2" = expression("dc-melding " * smc^2))
  ) +
  # Color scale JUST for the true-value line
  scale_colour_manual(
    values = c("true value" = "red"),
    name = NULL
  ) +
  labs(x = NULL, y = NULL, title = expression(phi[10 * "," * 11]), fill = 'method') +
  scale_x_discrete(labels = NULL) +
  theme(axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 40),
        plot.title = element_text(size = 40, hjust = .5),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 40),
        legend.spacing.y = unit(.5, "cm"),
        panel.grid.major.y = element_blank())



pdf(file = "C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/boxplot_phi10_11.pdf",
    width = 18, height = 12)

phi10_11_boxplot

dev.off()



# boxplot for phi34 and phi89 together

pdf(file = "C:/Users/Yixuan/Documents/codes/SMC/dc_melding/11_submodel/simulation/results/boxplot_phi12_10_11.pdf",
    width = 18, height = 12)

(phi12_boxplot | phi10_11_boxplot)

dev.off()



















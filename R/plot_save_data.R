library(tidyverse)

source("R/plot_funs.R")

file_path <- "C:/Users/guen.grosklos/Google Drive/Missoula_postdoc/Images/"
fig_colors <- c("#2ca25f", "#fc8d59", "#67a9cf", "#f768a1", "#bae4b3", "#fed98e")
nind <- 100
#########################################
# Individual plots
#########################################
# for (sim_num in 1) {
# sim_vars <- data.frame(
#   sim_names = c("Original", "Slow_landscape", "Medium_landscape", "Fast_landscape", "Slow_cams", "Medium_cams", "Fast_cams"),
#   lscape_tag = c("Random", rep("Homogeneous", 3), rep("Random", 3)),
#   all_speed = c(1, 1, 2, 3, rep(1, 3)),
#   cam.dist.set = c(rep(1, 4), 2, 3, 4)
# )
# 
# 
# D.all <- readRDS(paste0("Sim_results/Sim_", sim_vars$sim_names[sim_num], "_slower_speeds.rds"))
# 
# nind <- 100
# 
# lscape_tag <- sim_vars$lscape_tag[sim_num]
# all_speed <- sim_vars$all_speed[sim_num] # 1: Slow, 2: Medium, 3: Fast
# 
# # Cam sample designs (1: random, 2-4: biased slow, medium, fast)
# cam.dist.set <- sim_vars$cam.dist.set[sim_num]
# 
# cam.props.label <- c("Camera Bias: Random",
#                      "Camera Bias: Slow",
#                      "Camera Bias: Medium",
#                      "Camera Bias: Fast")
# 
# 
# D.all <- D.all %>% 
#   dplyr::bind_rows(tibble::tibble(
#     Model = c("TDST", "STE"),
#     Covariate = c("Non-Covariate", "Covariate"),
#     Est = c(nind, nind),
#     SD = c(mean(D.all$SD), mean(D.all$SD))
#   ))
# D.all$Model <- factor(D.all$Model, levels = c("TDST", "REST", "TTE", "PR", "STE"))
# 
# # Plot boxplots of means
# plot_multirun_means()
# 
# # ggsave(
# #   paste0("Abundance_", sim_vars$sim_names[sim_num], "_", sim_vars$lscape_tag[sim_num], ".png"),
# #   plot = last_plot(),
# #   path = file_path,
# #   scale = 1,
# #   width = NA,
# #   height = NA,
# #   # units = c("in", "cm", "mm", "px"),
# #   dpi = 600,
# #   limitsize = TRUE,
# #   bg = NULL
# # )
# 
# # Plot histograms of means?
# 
# # Plot boxplots of coefficient of variation
# plot_multirun_CV()
# 
# # ggsave(
# #   paste0("CV_", sim_vars$sim_names[sim_num], "_", sim_vars$lscape_tag[sim_num], ".png"),
# #   plot = last_plot(),
# #   path = file_path,
# #   scale = 1,
# #   width = NA,
# #   height = NA,
# #   # units = c("in", "cm", "mm", "px"),
# #   dpi = 600,
# #   limitsize = TRUE,
# #   bg = NULL
# # )
# 
# }

#########################################
# Grouped plots
#########################################
D.all <- data.frame(Model = NA,
                    Covariate = NA,
                    Est = NA,
                    SD = NA,
                    Prop_speeds = NA,
                    Run = NA
)
for (sim_num in 2:4) {
# for (sim_num in c(1, 5:7)) {
    sim_vars <- data.frame(
      sim_names = c("Original", "Slow_landscape", "Medium_landscape", "Fast_landscape", "Slow_cams", "Medium_cams", "Fast_cams"),
      lscape_tag = c("Random", rep("Homogeneous", 3), rep("Random", 3)),
      all_speed = c(1, 1, 2, 3, rep(1, 3)),
      cam.dist.set = c(rep(1, 4), 2, 3, 4)
    )
  
  D.in <- readRDS(paste0("Sim_results/Sim_", sim_vars$sim_names[sim_num], ".rds"))
  
  D.all <- D.all %>% 
    dplyr::bind_rows(c(D.in, Run = sim_vars$sim_names[sim_num])) %>% 
    dplyr::mutate(Model = as.character(Model))
  D.all$Model <- replace(D.all$Model, D.all$Model == "MCT", "PR") # Change old PR model name
  
  D.all <- D.all %>% 
    dplyr::bind_rows(tibble::tibble(
      Model = c("TDST", "STE"),
      Covariate = c("Non-Covariate", "Covariate"),
      Est = c(nind, nind),
      SD = c(mean(D.all$SD), mean(D.all$SD)),
      Run = sim_vars$sim_names[sim_num])
    )
}
# Remove NA row
D.all <- D.all[-1,]

D.all$Run[D.all$Run == "Original"] <- "Random"
  
D.all$Model <- factor(D.all$Model, levels = c("TDST", "REST", "TTE", "PR", "STE"))

if (sim_num %in% 2:4) {
  D.all$Run[D.all$Run == "Slow_landscape"] <- "Slow"
  D.all$Run[D.all$Run == "Medium_landscape"] <- "Medium"
  D.all$Run[D.all$Run == "Fast_landscape"] <- "Fast"
  D.all$Run <- factor(D.all$Run, levels = c("Slow", "Medium", "Fast"))
  
  plot_grouped_multirun_means(Unused_cov = "None",
                              subfig_label = "a")
  ggsave(
    "Abundance_homogeneous_speeds.png",
    plot = last_plot(),
    path = file_path,
    scale = 1,
    width = 5,
    height = 3,
    # units = c("in", "cm", "mm", "px"),
    dpi = 600,
    limitsize = TRUE,
    bg = NULL
  )
  
  plot_grouped_multirun_CV(Unused_cov = "None",
                           subfig_label = "b",
                           Title = "Habitat Type")
  ggsave(
    "SD_homogeneous_speeds.png",
    plot = last_plot(),
    path = file_path,
    scale = 1,
    width = 5,
    height = 3,
    # units = c("in", "cm", "mm", "px"),
    dpi = 600,
    limitsize = TRUE,
    bg = NULL
  )
  
} else if (sim_num %in% 5:7){
  D.all$Run[D.all$Run == "Slow_cams"] <- "80% Slow Habitat"
  D.all$Run[D.all$Run == "Medium_cams"] <- "80% Medium Habitat"
  D.all$Run[D.all$Run == "Fast_cams"] <- "80% Fast Habitat"
  D.all$Run <- factor(D.all$Run, levels = c("Random", "80% Slow Habitat", "80% Medium Habitat", "80% Fast Habitat"))
  
  plot_grouped_multirun_means(Unused_cov = "Non-Covariate", 
                              Filter_model = "STE",
                              subfig_label = "c",
                              Cov = "Covariate")
  ggsave(
    "Abundance_cam_speeds_covariate.png",
    plot = last_plot(),
    path = file_path,
    scale = 1,
    width = 5,
    height = 3,
    # units = c("in", "cm", "mm", "px"),
    dpi = 600,
    limitsize = TRUE,
    bg = NULL
  )
  
  plot_grouped_multirun_CV(Unused_cov = "Non-Covariate", 
                           Filter_model = "STE",
                           subfig_label = "d",
                           Title = "Sample Design")
  ggsave(
    "SD_cam_speeds_covariate.png",
    plot = last_plot(),
    path = file_path,
    scale = 1,
    width = 5,
    height = 3,
    # units = c("in", "cm", "mm", "px"),
    dpi = 600,
    limitsize = TRUE,
    bg = NULL
  )
  
  plot_grouped_multirun_means(Unused_cov = "Covariate", 
                              Filter_model = "TDST",
                              subfig_label = "a",
                              Cov = "Non-Covariate")
  ggsave(
    "Abundance_cam_speeds_no_covariate.png",
    plot = last_plot(),
    path = file_path,
    scale = 1,
    width = 5,
    height = 3,
    # units = c("in", "cm", "mm", "px"),
    dpi = 600,
    limitsize = TRUE,
    bg = NULL
  )
  
  plot_grouped_multirun_CV(Unused_cov = "Covariate",
                           Filter_model = "TDST",
                           subfig_label = "b",
                           Title = "Sample Design")
  ggsave(
    "SD_cam_speeds_no_covariate.png",
    plot = last_plot(),
    path = file_path,
    scale = 1,
    width = 5,
    height = 3,
    # units = c("in", "cm", "mm", "px"),
    dpi = 600,
    limitsize = TRUE,
    bg = NULL
  )
}

#########################################
# Camera effort
#########################################
D.all <- data.frame(Model = NA,
                    Covariate = NA,
                    Est = NA,
                    SD = NA,
                    Prop_speeds = NA,
                    ncam = NA
)

ncam_all <- c(10, 25, 50, 100, 150, 250)

for (ncams in ncam_all) {
  if (ncams == 250) {
    D.in <- readRDS(paste0("Sim_results/Sim_Original.rds"))
  } else {
    D.in <- readRDS(paste0("Sim_results/Sim_Original_", ncams, "cam.rds"))
  }
  
  D.all <- D.all %>% 
    dplyr::bind_rows(c(D.in, ncam = ncams)) %>% 
    dplyr::mutate(Model = as.character(Model))
  D.all$Model <- replace(D.all$Model, D.all$Model == "MCT", "PR") # Change old PR model name
  
}
D.all <- D.all[-1,]

D.all.summary <- D.all %>% 
  group_by(Model, Covariate, ncam) %>% 
  summarise(Means = mean(Est, na.rm = T),
            Mean_sd = sd(Est, na.rm = T),
            SDs = mean(SD, na.rm = T),
            CV = SDs/Means,
            MPE = 1/n() * sum((Means - nind)/nind),
            MSE = 1/n() * sum((Means - nind)^2),
            .groups = 'drop')

D.all.summary$Model <- factor(D.all.summary$Model, 
                              levels = c("TDST", "REST", "TTE", "PR", "STE"))

pd <- position_dodge(width = 2)

D.all.summary %>% 
  dplyr::filter(Covariate == "Covariate") %>% 
  # dplyr::filter(Model != "TTE") %>% 
  ggplot(aes(x = ncam, y = Means, shape = Model)) +
  geom_hline(yintercept=nind, linetype="dashed", size=1) +
  geom_line(size = 1) +
  geom_point(size = 2.5, stroke = 1.5) +
  scale_color_manual(values = fig_colors[2:5]) +
  labs(x = "Number of Cameras",
       y = paste0("Covariate Models \n", "Mean Estimates")) +
  scale_x_continuous(breaks = ncam_all) +
  scale_shape_manual(values=c(1, 2, 3, 4)) +
  annotate("text", x = 250, y = 115, label = "c", size = 5) +
  theme(text = element_text(size = 16),
        legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none",
        legend.background = element_blank(),
        legend.spacing.y = unit(0, "mm"),
        legend.box.background = element_rect(colour = "black"))

ggsave(
  "Multi_cam_covariate_Abundance.png",
  plot = last_plot(),
  path = file_path,
  scale = 1,
  width = 6,
  height = 3,
  # units = c("in", "cm", "mm", "px"),
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

D.all.summary %>% 
  dplyr::filter(Covariate == "Covariate") %>% 
  ggplot(aes(x = ncam, y = SDs, shape = Model)) +
  geom_line(size = 1) +
  geom_point(size = 2.5, stroke = 1.5) +
  labs(x = "Number of Cameras",
       y = "Mean SDs") +
  scale_x_continuous(breaks = ncam_all) +
  scale_shape_manual(values=c(1, 2, 3, 4)) +
  annotate("text", x = 250, y = 53, label = "c", size = 5) +
  theme(text = element_text(size = 16),
        legend.title=element_text(size=12),
        legend.text=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = c(0.21, 0.74),
        legend.background = element_blank(),
        legend.spacing.y = unit(0, "mm"))
ggsave(
  "Multi_cam_covariate_SD.png",
  plot = last_plot(),
  path = file_path,
  scale = 1,
  width = 6,
  height = 3,
  # units = c("in", "cm", "mm", "px"),
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

D.all.summary %>% 
  dplyr::filter(Covariate == "Non-Covariate") %>% 
  # dplyr::filter(Model != "TTE") %>% 
  ggplot(aes(x = ncam, y = Means, shape = Model)) +
  geom_hline(yintercept=nind, linetype="dashed", size=1) +
  geom_line(size = 1) +
  geom_point(size = 2.5, stroke = 1.5) +
  scale_color_manual(values = fig_colors[2:5]) +
  labs(x = "Number of Cameras",
       y = paste0("Non-Covariate Models \n", "Mean Estimates")) +
  scale_x_continuous(breaks = ncam_all) +
  scale_shape_manual(values=c(2, 3, 4, 5)) +
  annotate("text", x = 250, y = 105, label = "a", size = 5) +
  theme(text = element_text(size = 16),
        legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none",
        legend.background = element_blank(),
        legend.spacing.y = unit(0, "mm"))
ggsave(
  "Multi_cam_noncovariate_Abundance.png",
  plot = last_plot(),
  path = file_path,
  scale = 1,
  width = 6,
  height = 3,
  # units = c("in", "cm", "mm", "px"),
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

D.all.summary %>% 
  dplyr::filter(Covariate == "Non-Covariate") %>% 
  ggplot(aes(x = ncam, y = SDs, shape = Model)) +
  geom_line(size = 1) +
  geom_point(size = 2.5, stroke = 1.5) +
  labs(x = "Number of Cameras",
       y = "Mean SDs") +
  scale_x_continuous(breaks = ncam_all) +
  scale_shape_manual(values=c(2, 3, 4, 5)) +
  annotate("text", x = 250, y = 25, label = "b", size = 5) +
  theme(text = element_text(size = 16),
        legend.title=element_text(size=12),
        legend.text=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = c(0.21, 0.74),
        legend.background = element_blank(),
        legend.spacing.y = unit(0, "mm"))
ggsave(
  "Multi_cam_noncovariate_SD.png",
  plot = last_plot(),
  path = file_path,
  scale = 1,
  width = 6,
  height = 3,
  # units = c("in", "cm", "mm", "px"),
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)




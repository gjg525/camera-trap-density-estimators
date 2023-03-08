library(tidyverse)

source("R/plot_funs.R")

file_path <- "C:/Users/guen.grosklos/Google Drive/Missoula_postdoc/Images/"

sim_num <- 3
for (sim_num in 3:7) {
sim_vars <- data.frame(
  sim_names = c("Original", "Slow_landscape", "Medium_landscape", "Fast_landscape", "Slow_cams", "Medium_cams", "Fast_cams"),
  lscape_tag = c("Random", rep("Homogeneous", 3), rep("Random", 3)),
  all_speed = c(1, 1, 2, 3, rep(1, 3)),
  cam.dist.set = c(rep(1, 4), 2, 3, 4)
)

D.all <- readRDS(paste0("Sim_results/Sim_", sim_vars$sim_names[sim_num], ".rds"))

nind <- 100

lscape_tag <- sim_vars$lscape_tag[sim_num]
all_speed <- sim_vars$all_speed[sim_num] # 1: Slow, 2: Medium, 3: Fast

# Cam sample designs (1: random, 2-4: biased slow, medium, fast)
cam.dist.set <- sim_vars$cam.dist.set[sim_num]

cam.props.label <- c("Camera Bias: Random",
                     "Camera Bias: Slow",
                     "Camera Bias: Medium",
                     "Camera Bias: Fast")


D.all <- D.all %>% 
  dplyr::bind_rows(tibble::tibble(
    Model = c("TDST", "STE"),
    Covariate = c("Non-Covariate", "Covariate"),
    Est = c(nind, nind),
    SD = c(mean(D.all$SD), mean(D.all$SD))
  ))
D.all$Model <- factor(D.all$Model, levels = c("TDST", "REST", "TTE", "MCT", "STE"))

# Plot boxplots of means
plot_multirun_means()

ggsave(
  paste0("Abundance_", sim_vars$sim_names[sim_num], "_", sim_vars$lscape_tag[sim_num], ".png"),
  plot = last_plot(),
  path = file_path,
  scale = 1,
  width = NA,
  height = NA,
  # units = c("in", "cm", "mm", "px"),
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

# Plot histograms of means?

# Plot boxplots of coefficient of variation
plot_multirun_CV()

ggsave(
  paste0("CV_", sim_vars$sim_names[sim_num], "_", sim_vars$lscape_tag[sim_num], ".png"),
  plot = last_plot(),
  path = file_path,
  scale = 1,
  width = NA,
  height = NA,
  # units = c("in", "cm", "mm", "px"),
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

}

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
nind <- 100
for (sim_num in 2:4) {
  sim_vars <- data.frame(
    sim_names = c("Original", "Slow_landscape", "Medium_landscape", "Fast_landscape", "Slow_cams", "Medium_cams", "Fast_cams"),
    lscape_tag = c("Random", rep("Homogeneous", 3), rep("Random", 3)),
    all_speed = c(1, 1, 2, 3, rep(1, 3)),
    cam.dist.set = c(rep(1, 4), 2, 3, 4)
  )
  
  D.in <- readRDS(paste0("Sim_results/Sim_", sim_vars$sim_names[sim_num], ".rds"))
  
  D.all <- D.all %>% 
    dplyr::bind_rows(c(D.in, Run = sim_vars$sim_names[sim_num]))
  
  D.all <- D.all %>% 
    dplyr::bind_rows(tibble::tibble(
      Model = c("TDST", "STE"),
      Covariate = c("Non-Covariate", "Covariate"),
      Est = c(nind, nind),
      SD = c(mean(D.all$SD), mean(D.all$SD)),
      Run = sim_vars$sim_names[sim_num])
    )
}
D.all <- D.all[-1,]

  
lscape_tag <- sim_vars$lscape_tag[sim_num]
all_speed <- sim_vars$all_speed[sim_num] # 1: Slow, 2: Medium, 3: Fast

# Cam sample designs (1: random, 2-4: biased slow, medium, fast)
cam.dist.set <- sim_vars$cam.dist.set[sim_num]

cam.props.label <- c("Camera Bias: Random",
                     "Camera Bias: Slow",
                     "Camera Bias: Medium",
                     "Camera Bias: Fast")


D.all$Model <- factor(D.all$Model, levels = c("TDST", "REST", "TTE", "MCT", "STE"))
if (sim_num %in% 2:4) {
  D.all$Run <- factor(D.all$Run, levels = c("Slow_landscape", "Medium_landscape", "Fast_landscape"))
} else if (sim_num %in% 5:7){
  D.all$Run <- factor(D.all$Run, levels = c("Slow_cams", "Medium_cams", "Fast_cams"))
}

D.all %>% 
  dplyr::filter(Covariate == "Covariate") %>%
  ggplot(aes(x = Model, y = Est, fill = Run)) +
  geom_boxplot(lwd = .1, fatten = .1) +
  geom_hline(yintercept=100, linetype="dashed", size=1) +
  labs(x = "Model",
       # y = paste(cam.props.label[lscape_var+1], " Landscape \n Abundance Estimates")) +
       y = "Abundance Estimates") +
  scale_fill_manual(values=c('grey20', 'grey40','Grey')) +
  theme(text = element_text(size = 20),
        legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = c(0.17, 0.84),
        legend.background = element_blank(),
        legend.spacing.y = unit(0, "mm"),
        legend.box.background = element_rect(colour = "black"))

D.all %>% 
  dplyr::filter(Covariate == "Non-Covariate") %>%
  ggplot(aes(x = Model, y = SD/Est, fill = Run)) +
  geom_boxplot(lwd = .1, fatten = .1) +
  labs(x = "Model",
       # y = paste(cam.props.label[lscape_var+1], " Landscape \n Abundance Estimates")) +
       y = "Abundance Estimates") +
  scale_fill_manual(values=c('grey20', 'grey40','Grey')) +
  theme(text = element_text(size = 20),
        legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = c(0.17, 0.84),
        legend.background = element_blank(),
        legend.spacing.y = unit(0, "mm"),
        legend.box.background = element_rect(colour = "black"))



  # Plot boxplots of means
  # plot_multirun_means()
  
  # ggsave(
  #   paste0("Abundance_", sim_vars$sim_names[sim_num], "_", sim_vars$lscape_tag[sim_num], ".png"),
  #   plot = last_plot(),
  #   path = file_path,
  #   scale = 1,
  #   width = NA,
  #   height = NA,
  #   # units = c("in", "cm", "mm", "px"),
  #   dpi = 600,
  #   limitsize = TRUE,
  #   bg = NULL
  # )
  
  # Plot histograms of means?
  
  # Plot boxplots of coefficient of variation
  # plot_multirun_CV()
  
  # ggsave(
  #   paste0("CV_", sim_vars$sim_names[sim_num], "_", sim_vars$lscape_tag[sim_num], ".png"),
  #   plot = last_plot(),
  #   path = file_path,
  #   scale = 1,
  #   width = NA,
  #   height = NA,
  #   # units = c("in", "cm", "mm", "px"),
  #   dpi = 600,
  #   limitsize = TRUE,
  #   bg = NULL
  # )
  


# Plot different camera efforts
####################
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
    dplyr::bind_rows(c(D.in, ncam = ncams))
}
D.all <- D.all[-1,]

D.all.summary <- D.all %>% 
  group_by(Model, Covariate, ncam) %>% 
  summarise(Means = mean(Est, na.rm = T),
            SDs = mean(SD, na.rm = T),
            CV = SDs/Means,
            MPE = 1/n() * sum((Means - nind)/nind),
            MSE = 1/n() * sum((Means - nind)^2),
            .groups = 'drop')

pd <- position_dodge(width = 2)

D.all.summary %>% 
  dplyr::filter(Covariate == "Covariate") %>% 
  dplyr::filter(Model != "TTE") %>% 
  ggplot(aes(x = ncam, y = Means, color = Model)) +
  labs(x = "Number of Cameras",
       y = "Mean") +
  geom_line() +
  geom_point() +
  # geom_errorbar(aes(ymin=Means-SDs, ymax=Means+SDs), width=.2,
  #               position=pd) +
  theme(text = element_text(size = 20),
        legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = c(0.17, 0.84),
        legend.background = element_blank(),
        legend.spacing.y = unit(0, "mm"),
        legend.box.background = element_rect(colour = "black"))

D.all.summary %>% 
  dplyr::filter(Covariate == "Non-Covariate") %>% 
  dplyr::filter(Model != "TTE") %>% 
  ggplot(aes(x = ncam, y = Means, color = Model)) +
  labs(x = "Number of Cameras",
       y = "Mean") +
  geom_line() +
  geom_point() +
  theme(text = element_text(size = 20),
        legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = c(0.17, 0.84),
        legend.background = element_blank(),
        legend.spacing.y = unit(0, "mm"),
        legend.box.background = element_rect(colour = "black"))

D.all.summary %>% 
  dplyr::filter(Covariate == "Covariate") %>% 
  ggplot(aes(x = ncam, y = SDs, color = Model)) +
  labs(x = "Number of Cameras",
       y = "CV") +
  geom_line() +
  geom_point() +
  theme(text = element_text(size = 20),
        legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = c(0.17, 0.84),
        legend.background = element_blank(),
        legend.spacing.y = unit(0, "mm"),
        legend.box.background = element_rect(colour = "black"))

D.all.summary %>% 
  dplyr::filter(Covariate == "Non-Covariate") %>% 
  ggplot(aes(x = ncam, y = SDs, color = Model)) +
  labs(x = "Number of Cameras",
       y = "CV") +
  geom_line() +
  geom_point() +
  theme(text = element_text(size = 20),
        legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = c(0.17, 0.84),
        legend.background = element_blank(),
        legend.spacing.y = unit(0, "mm"),
        legend.box.background = element_rect(colour = "black"))




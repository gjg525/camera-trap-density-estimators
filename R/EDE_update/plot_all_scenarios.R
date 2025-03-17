library(dplyr)
library(ggplot2)
tot_animals <- 100

fig_colors <- c("#1B5E20", "#00A8C6", "#FBC02D", "#E65100", "#8E44AD", "#4B6FAD", "#D81B60")

file_names <- c(
  "random",
  "slow",
  "med",
  "fast",
  "all_slow",
  "all_med",
  "all_fast"
)
################################################################################
# loadRData <- function(fileName){
#   #loads an RData file, and returns it
#   load(fileName)
#   get(ls()[ls() != "fileName"])
# }
# 
# source("./R/EDE_update/utils.R")
# source("./R/EDE_update/plot_funs.R")
# 
# 
# D_all <- tibble::tibble()
# for (ii in 1:length(file_names)) {
#   print(ii)
#   all_results <- loadRData(paste0("G:/My Drive/Missoula_postdoc/TDST_Model/Data/results_",
#                                   file_names[ii], 
#                                   "_cam_alt.RData"))
# 
#   D_all <- dplyr::bind_rows(
#     D_all,
#     all_results[[5]] %>%
#       dplyr::mutate(
#         SampDesign = paste0(file_names[ii], "_cam")
#       )
#   )
# }
# 
# save(D_all, file = "G:/My Drive/Missoula_postdoc/TDST_Model/Data/D_all.RData")

# load("G:/My Drive/Missoula_postdoc/TDST_Model/Data/D_all.RData")
load("C:/Users/guen.grosklos/Google Drive/Missoula_postdoc/TDST_Model/Data/D_all.RData")

D_all <- D_all %>% 
  dplyr::filter(Model %in% c("PR", "PR Habitat")) %>% 
  dplyr::mutate(
    Model = paste(Model, Covariate),
    Model = dplyr::case_when(
      Model == "PR Non-Covariate" ~ "PD",
      Model == "PR Covariate" ~ "PR",
      Model == "PR Habitat Non-Covariate" ~ "HAM" # or PDH?
    )
  )

D_all$Model <- factor(
  D_all$Model,
  levels = c("PD", "HAM", "PR")
)

# Random only
for (ff in 1){#:length(file_names)) {
  D_all %>% 
    dplyr::filter(Est < 250) %>%
    dplyr::filter(SampDesign == paste0(file_names[ff], "_cam")) %>% 
    ggplot2::ggplot(ggplot2::aes(x = Model, y = Est, fill = Model)) +
    ggplot2::geom_boxplot(lwd = 0.5, fatten = .5, outlier.shape = NA) +
    ggplot2::geom_hline(yintercept=tot_animals, linetype="dashed", size = 0.7) +
    ggplot2::labs(x = "Model",
                  y = "Posterior Means") +
    ggplot2::scale_fill_manual(values= fig_colors[1:5]) +
    ggplot2::scale_color_manual(values = c('grey0', 'grey40', 'grey60')) +
    ggplot2::theme(text = ggplot2::element_text(size = 16),
                   legend.title = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(colour = "black"),
                   panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
                   legend.position = "none",
                   legend.background = ggplot2::element_blank(),
                   legend.spacing.y = ggplot2::unit(0, "mm"),
                   legend.box.background = ggplot2::element_rect(colour = "black"))
  
  ggplot2::ggsave(
    paste0("C:/Users/guen.grosklos/Google Drive/Missoula_postdoc/TDST_Model/imgs/", 
           file_names[ff],
           "_cam.png"),
    plot = ggplot2::last_plot(),
    # path = file_path,
    # scale = 1,
    width = 5,
    height = 3,
    # units = c("in", "cm", "mm", "px"),
    dpi = 600,
    limitsize = TRUE,
    bg = NULL
  )
  
  D_all %>% 
    dplyr::filter(Est < 250) %>%
    dplyr::filter(SampDesign == paste0(file_names[ff], "_cam")) %>% 
    ggplot2::ggplot(ggplot2::aes(x = Model, y = SD, fill = Model)) +
    ggplot2::geom_boxplot(lwd = 0.5, fatten = .5, outlier.shape = NA) +
    ggplot2::labs(x = "Model",
                  y = "Posterior SDs") +
    ggplot2::scale_fill_manual(values= fig_colors[1:5]) +
    ggplot2::scale_color_manual(values = c('grey0', 'grey40', 'grey60')) +
    ggplot2::theme(text = ggplot2::element_text(size = 16),
                   legend.title = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(colour = "black"),
                   panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
                   legend.position = "none",
                   legend.background = ggplot2::element_blank(),
                   legend.spacing.y = ggplot2::unit(0, "mm"),
                   legend.box.background = ggplot2::element_rect(colour = "black"))
  
  ggplot2::ggsave(
    paste0("C:/Users/guen.grosklos/Google Drive/Missoula_postdoc/TDST_Model/imgs/", 
           file_names[ff],
           "_cam_SD.png"),
    plot = ggplot2::last_plot(),
    # path = file_path,
    # scale = 1,
    width = 5,
    height = 3,
    # units = c("in", "cm", "mm", "px"),
    dpi = 600,
    limitsize = TRUE,
    bg = NULL
  )
  
}

################################################################################
# All other sample designs
D_all %>% 
  dplyr::filter(Est < 250) %>%
  dplyr::filter(SampDesign %in% c("slow_cam", "med_cam", "fast_cam")) %>% 
  dplyr::mutate(
    SampDesign = dplyr::case_when(
      SampDesign == "slow_cam" ~ "80% Slow",
      SampDesign == "med_cam" ~ "80% Medium",
      SampDesign == "fast_cam" ~ "80% Fast"
    )
  ) %>% 
  ggplot2::ggplot(ggplot2::aes(x = SampDesign, y = Est, fill = Model)) +
  ggplot2::geom_boxplot(lwd = 0.5, fatten = .5, outlier.shape = NA) +
  ggplot2::geom_hline(yintercept=tot_animals, linetype="dashed", size = 0.7) +
  ggplot2::labs(x = "Camera Sample Design",
                y = "Posterior Means") +
  ggplot2::scale_fill_manual(values= fig_colors[1:5]) +
  ggplot2::scale_color_manual(values = c('grey0', 'grey40', 'grey60')) +
  ggplot2::theme(text = ggplot2::element_text(size = 16),
                 panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(),
                 axis.line = ggplot2::element_line(colour = "black"),
                 panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
                 legend.position = "none",
                 legend.background = ggplot2::element_blank(),
                 legend.spacing.y = ggplot2::unit(0, "mm"),
                 legend.box.background = ggplot2::element_rect(colour = "black"))

ggplot2::ggsave(
  paste0("C:/Users/guen.grosklos/Google Drive/Missoula_postdoc/TDST_Model/imgs/bias_cam.png"),
  plot = ggplot2::last_plot(),
  # path = file_path,
  # scale = 1,
  width = 5,
  height = 3,
  # units = c("in", "cm", "mm", "px"),
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

D_all %>% 
  dplyr::filter(Est < 250) %>%
  dplyr::filter(SampDesign %in% c("slow_cam", "med_cam", "fast_cam")) %>% 
  dplyr::mutate(
    SampDesign = dplyr::case_when(
      SampDesign == "slow_cam" ~ "80% Slow",
      SampDesign == "med_cam" ~ "80% Medium",
      SampDesign == "fast_cam" ~ "80% Fast"
    )
  ) %>% 
  ggplot2::ggplot(ggplot2::aes(x = SampDesign, y = SD, fill = Model)) +
  ggplot2::geom_boxplot(lwd = 0.5, fatten = .5, outlier.shape = NA) +
  ggplot2::labs(x = "Camera Sample Design",
                y = "Posterior SDs") +
  ggplot2::scale_fill_manual(values= fig_colors[1:5]) +
  ggplot2::scale_color_manual(values = c('grey0', 'grey40', 'grey60')) +
  ggplot2::theme(text = ggplot2::element_text(size = 16),
                 legend.title=element_text(size=12),
                 legend.text=element_text(size=10),
                 legend.position = c(0.915, 0.69),
                 panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(),
                 axis.line = ggplot2::element_line(colour = "black"),
                 panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
                 legend.background = ggplot2::element_rect(color = "white"),
                 legend.spacing.y = ggplot2::unit(0, "mm"),
                 legend.box.background = ggplot2::element_rect(colour = "black"))

ggplot2::ggsave(
  paste0("C:/Users/guen.grosklos/Google Drive/Missoula_postdoc/TDST_Model/imgs/bias_cam_SD.png"),
  plot = ggplot2::last_plot(),
  # path = file_path,
  # scale = 1,
  width = 5,
  height = 3,
  # units = c("in", "cm", "mm", "px"),
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

################################################################################
D_all %>% 
  dplyr::filter(Est < 250) %>%
  dplyr::filter(SampDesign %in% c("all_slow_cam", "all_med_cam", "all_fast_cam")) %>% 
  dplyr::mutate(
    SampDesign = dplyr::case_when(
      SampDesign == "all_slow_cam" ~ "100% Slow",
      SampDesign == "all_med_cam" ~ "100% Medium",
      SampDesign == "all_fast_cam" ~ "100% Fast"
    )
  ) %>% 
  ggplot2::ggplot(ggplot2::aes(x = SampDesign, y = Est, fill = Model)) +
  ggplot2::geom_boxplot(lwd = 0.5, fatten = .5, outlier.shape = NA) +
  ggplot2::geom_hline(yintercept=tot_animals, linetype="dashed", size = 0.7) +
  ggplot2::labs(x = "Camera Sample Design",
                y = "Posterior Means") +
  ggplot2::scale_fill_manual(values= fig_colors[1:5]) +
  ggplot2::scale_color_manual(values = c('grey0', 'grey40', 'grey60')) +
  ggplot2::theme(text = ggplot2::element_text(size = 16),
                 panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(),
                 axis.line = ggplot2::element_line(colour = "black"),
                 panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
                 legend.position = "none",
                 legend.background = ggplot2::element_blank(),
                 legend.spacing.y = ggplot2::unit(0, "mm"),
                 legend.box.background = ggplot2::element_rect(colour = "black"))

ggplot2::ggsave(
  paste0("C:/Users/guen.grosklos/Google Drive/Missoula_postdoc/TDST_Model/imgs/all_bias_cam.png"),
  plot = ggplot2::last_plot(),
  # path = file_path,
  # scale = 1,
  width = 5,
  height = 3,
  # units = c("in", "cm", "mm", "px"),
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)


D_all %>% 
  dplyr::filter(Est < 250 & SD < 40) %>%
  dplyr::filter(SampDesign %in% c("all_slow_cam", "all_med_cam", "all_fast_cam")) %>% 
  dplyr::mutate(
    SampDesign = dplyr::case_when(
      SampDesign == "all_slow_cam" ~ "100% Slow",
      SampDesign == "all_med_cam" ~ "100% Medium",
      SampDesign == "all_fast_cam" ~ "100% Fast"
    )
  ) %>% 
  ggplot2::ggplot(ggplot2::aes(x = SampDesign, y = SD, fill = Model)) +
  ggplot2::geom_boxplot(lwd = 0.5, fatten = .5, outlier.shape = NA) +
  ggplot2::labs(x = "Camera Sample Design",
                y = "Posterior SDs") +
  ggplot2::scale_fill_manual(values= fig_colors[1:5]) +
  ggplot2::scale_color_manual(values = c('grey0', 'grey40', 'grey60')) +
  ggplot2::theme(text = ggplot2::element_text(size = 16),
                 legend.title=element_text(size=12),
                 legend.text=element_text(size=10),
                 legend.position = c(0.915, 0.69),
                 panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(),
                 axis.line = ggplot2::element_line(colour = "black"),
                 panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
                 legend.background = ggplot2::element_rect(color = "white"),
                 legend.spacing.y = ggplot2::unit(0, "mm"),
                 legend.box.background = ggplot2::element_rect(colour = "black"))

ggplot2::ggsave(
  paste0("C:/Users/guen.grosklos/Google Drive/Missoula_postdoc/TDST_Model/imgs/all_bias_cam_SD.png"),
  plot = ggplot2::last_plot(),
  # path = file_path,
  # scale = 1,
  width = 5,
  height = 3,
  # units = c("in", "cm", "mm", "px"),
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

# Try all in one plot?
D_all %>% 
  dplyr::filter(Est < 250) %>%
  dplyr::mutate(
    SampDesign = dplyr::case_when(
      SampDesign == "slow_cam" ~ "80% Slow",
      SampDesign == "med_cam" ~ "80% Medium",
      SampDesign == "fast_cam" ~ "80% Fast",      
      SampDesign == "all_slow_cam" ~ "100% Slow",
      SampDesign == "all_med_cam" ~ "100% Medium",
      SampDesign == "all_fast_cam" ~ "100% Fast",
      .default = "Random"
      
    )
  ) %>% 
  ggplot2::ggplot(ggplot2::aes(x = Model, y = Est, fill = SampDesign)) +
  ggplot2::geom_boxplot(lwd = 0.5, fatten = .5, outlier.shape = NA) +
  ggplot2::geom_hline(yintercept=tot_animals, linetype="dashed", size = 0.7) +
  ggplot2::labs(x = "Camera Sample Design",
                y = "Posterior Means") +
  ggplot2::scale_fill_manual(values= fig_colors[1:7]) +
  ggplot2::scale_color_manual(values = c('grey0', 'grey40', 'grey60')) +
  ggplot2::theme(text = ggplot2::element_text(size = 16),
                 panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(),
                 axis.line = ggplot2::element_line(colour = "black"),
                 panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
                 # legend.position = "none",
                 # legend.background = ggplot2::element_blank(),
                 # legend.spacing.y = ggplot2::unit(0, "mm"),
                 # legend.box.background = ggplot2::element_rect(colour = "black")
  )

D_all %>% 
  dplyr::filter(Est < 250) %>%
  dplyr::mutate(
    SampDesign = dplyr::case_when(
      SampDesign == "slow_cam" ~ "80% Slow",
      SampDesign == "med_cam" ~ "80% Medium",
      SampDesign == "fast_cam" ~ "80% Fast",      
      SampDesign == "all_slow_cam" ~ "100% Slow",
      SampDesign == "all_med_cam" ~ "100% Medium",
      SampDesign == "all_fast_cam" ~ "100% Fast",
      .default = "Random"
      
    )
  ) %>% 
  ggplot2::ggplot(ggplot2::aes(x = Model, y = SD / Est, fill = SampDesign)) +
  ggplot2::geom_boxplot(lwd = 0.5, fatten = .5, outlier.shape = NA) +
  ggplot2::labs(x = "Camera Sample Design",
                y = "Posterior CVs") +
  ggplot2::scale_fill_manual(values= fig_colors[1:7]) +
  ggplot2::scale_color_manual(values = c('grey0', 'grey40', 'grey60')) +
  ggplot2::theme(text = ggplot2::element_text(size = 16),
                 panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(),
                 axis.line = ggplot2::element_line(colour = "black"),
                 panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
                 # legend.position = "none",
                 # legend.background = ggplot2::element_blank(),
                 # legend.spacing.y = ggplot2::unit(0, "mm"),
                 # legend.box.background = ggplot2::element_rect(colour = "black")
  )

# D_all %>%
#   dplyr::filter(Est < 200) %>%
#   # dplyr::filter(SampDesign == "random_cam") %>%
#   ggplot2::ggplot(ggplot2::aes(x = Model, y = Est, fill = Covariate)) +
#   ggplot2::geom_boxplot(position = ggplot2::position_dodge2(preserve = "single"),
#                         outlier.shape = NA) +
#   # ggplot2::geom_violin() +
#   ggplot2::geom_hline(yintercept = tot_animals, linetype = "dashed", size = 1) +
#   ggplot2::labs(
#     x = "Model",
#     y = "Mean Abundance"
#   ) +
#   # ggplot2::coord_cartesian(ylim = quantile(D_all$Est, c(0.01, 0.99), na.rm = T)) +
#   # scale_y_continuous(limits = c(lower, upper)) +
#   ggplot2::theme(
#     text = ggplot2::element_text(size = 20),
#     legend.title = ggplot2::element_blank(),
#     panel.grid.major = ggplot2::element_blank(),
#     panel.grid.minor = ggplot2::element_blank(),
#     panel.background = ggplot2::element_blank(),
#     axis.line = ggplot2::element_line(colour = "black"),
#     panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 1),
#     # legend.position = c(0.85, 0.84),
#     legend.position = c(0.8, 0.05),
#     legend.background = ggplot2::element_blank(),
#     legend.spacing.y = ggplot2::unit(0, "mm"),
#     legend.box.background = ggplot2::element_rect(colour = "black")
#   ) +
#   ggplot2::facet_wrap(~ SampDesign, scales = "free", ncol = 2)
# 
# D_all %>% 
#   dplyr::filter(Est < 200) %>% 
#   dplyr::mutate(CV = SD / Est) %>% 
#   # dplyr::filter(SampDesign == "random_cam") %>%
#   ggplot2::ggplot(ggplot2::aes(x = Model, y = CV, fill = Covariate)) + 
#   ggplot2::geom_boxplot(position = ggplot2::position_dodge2(preserve = "single"),
#                         outlier.shape = NA) +
#   ggplot2::labs(
#     x = "Model",
#     y = "Mean Abundance"
#   ) +
#   # ggplot2::coord_cartesian(ylim = quantile(D_all$Est, c(0.01, 0.99), na.rm = T)) +
#   # scale_y_continuous(limits = c(lower, upper)) +
#   ggplot2::theme(
#     text = ggplot2::element_text(size = 20),
#     legend.title = ggplot2::element_blank(),
#     panel.grid.major = ggplot2::element_blank(),
#     panel.grid.minor = ggplot2::element_blank(),
#     panel.background = ggplot2::element_blank(),
#     axis.line = ggplot2::element_line(colour = "black"),
#     panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 1),
#     # legend.position = c(0.85, 0.84),
#     legend.position = c(0.8, 0.05),
#     legend.background = ggplot2::element_blank(),
#     legend.spacing.y = ggplot2::unit(0, "mm"),
#     legend.box.background = ggplot2::element_rect(colour = "black")
#   ) + 
#   ggplot2::facet_wrap(~ SampDesign, scales = "free", ncol = 2)

# D_all %>% 
#   # dplyr::filter(Est < 250) %>%
#   dplyr::filter(SampDesign == "random_cam") %>% 
#   ggplot2::ggplot(ggplot2::aes(x = Model, y = Est, fill = Model)) +
#   ggplot2::geom_boxplot(lwd = 0.5, fatten = .5, outlier.size = 1) +
#   ggplot2::geom_hline(yintercept=tot_animals, linetype="dashed", size = 0.7) +
#   ggplot2::labs(x = "Model",
#        y = "Posterior Means") +
#   # scale_fill_manual(values= c(fig_colors_dark[1:5], 
#   #                             fig_colors_med[1:5], 
#   #                             fig_colors_light[1:5])) +
#   ggplot2::scale_fill_manual(values= fig_colors[1:5]) +
#   ggplot2::scale_color_manual(values = c('grey0', 'grey40', 'grey60')) +
#   # ggplot2::annotate("text", x = 5.4, y = max(D.all$Est), 
#   #          label = "a", 
#   #          size = 5) +
#   ggplot2::theme(text = ggplot2::element_text(size = 16),
#         legend.title = ggplot2::element_blank(),
#         panel.grid.major = ggplot2::element_blank(),
#         panel.grid.minor = ggplot2::element_blank(),
#         panel.background = ggplot2::element_blank(),
#         axis.line = ggplot2::element_line(colour = "black"),
#         panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
#         legend.position = "none",
#         legend.background = ggplot2::element_blank(),
#         legend.spacing.y = ggplot2::unit(0, "mm"),
#         legend.box.background = ggplot2::element_rect(colour = "black"))
# 
# ggsave(
#   "Abundance_homogeneous_speeds.png",
#   plot = last_plot(),
#   # path = file_path,
#   scale = 1,
#   width = 5,
#   height = 3,
#   # units = c("in", "cm", "mm", "px"),
#   dpi = 600,
#   limitsize = TRUE,
#   bg = NULL
# )
# 
# D_all %>% 
#   dplyr::filter(Est < 250) %>%
#   dplyr::filter(SampDesign == "random_cam") %>% 
#   ggplot2::ggplot(ggplot2::aes(x = Model, y = SD, fill = Model)) +
#   ggplot2::geom_boxplot(lwd = 0.5, fatten = .5, outlier.size = 1) +
#   ggplot2::labs(x = "Model",
#                 y = "Posterior SDs") +
#   # scale_fill_manual(values= c(fig_colors_dark[1:5], 
#   #                             fig_colors_med[1:5], 
#   #                             fig_colors_light[1:5])) +
#   ggplot2::scale_fill_manual(values= fig_colors[1:5]) +
#   ggplot2::scale_color_manual(values = c('grey0', 'grey40', 'grey60')) +
#   # ggplot2::annotate("text", x = 5.4, y = max(D.all$Est), 
#   #          label = "a", 
#   #          size = 5) +
#   ggplot2::theme(text = ggplot2::element_text(size = 16),
#                  legend.title = ggplot2::element_blank(),
#                  panel.grid.major = ggplot2::element_blank(),
#                  panel.grid.minor = ggplot2::element_blank(),
#                  panel.background = ggplot2::element_blank(),
#                  axis.line = ggplot2::element_line(colour = "black"),
#                  panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
#                  legend.position = "none",
#                  legend.background = ggplot2::element_blank(),
#                  legend.spacing.y = ggplot2::unit(0, "mm"),
#                  legend.box.background = ggplot2::element_rect(colour = "black"))
# 
# ggsave(
#   "SD_homogeneous_speeds.png",
#   plot = last_plot(),
#   # path = file_path,
#   scale = 1,
#   width = 5,
#   height = 3,
#   # units = c("in", "cm", "mm", "px"),
#   dpi = 600,
#   limitsize = TRUE,
#   bg = NULL
# )


library(dplyr)
tot_animals <- 100

fig_colors <- c("#2ca25f", "#fc8d59", "#67a9cf", "#f768a1", "#bae4b3", "#fed98e")

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

load("G:/My Drive/Missoula_postdoc/TDST_Model/Data/D_all.RData")

D_all <- D_all %>% 
  dplyr::mutate(
    Model = dplyr::case_when(
      Model == "PR" & Covariate == "Covariate" ~ "PR Covariate",
      Model == "PR" & Covariate == "Non-Covariate" ~ "PR Non-Covariate",
      # Model == "PR Habitat" ~ "Habitat PR",
      Model == "TDST" ~ "TDST Camera",
      Model == "TDST w/ Priors" ~ "TDST Tele Priors",
      .default = Model
    )
  )

for (ff in 1:length(file_names)) {
  D_all %>% 
    dplyr::filter(Est < 250) %>%
    dplyr::filter(SampDesign == paste0(file_names[ff], "_cam")) %>% 
    ggplot2::ggplot(ggplot2::aes(x = Model, y = Est, fill = Model)) +
    ggplot2::geom_boxplot(lwd = 0.5, fatten = .5, outlier.size = 1) +
    ggplot2::geom_hline(yintercept=tot_animals, linetype="dashed", size = 0.7) +
    ggplot2::labs(x = "Model",
                  y = "Posterior Means") +
    ggplot2::scale_fill_manual(values= fig_colors[1:5]) +
    ggplot2::scale_color_manual(values = c('grey0', 'grey40', 'grey60')) +
    ggplot2::theme(text = ggplot2::element_text(size = 20),
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
    paste0("G:/My Drive/Missoula_postdoc/TDST_Model/imgs/", 
           file_names[ff],
           "_cam.png"),
    plot = ggplot2::last_plot(),
    # path = file_path,
    # scale = 1,
    width = 8,
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
    ggplot2::geom_boxplot(lwd = 0.5, fatten = .5, outlier.size = 1) +
    ggplot2::labs(x = "Model",
                  y = "Posterior Means") +
    ggplot2::scale_fill_manual(values= fig_colors[1:5]) +
    ggplot2::scale_color_manual(values = c('grey0', 'grey40', 'grey60')) +
    ggplot2::theme(text = ggplot2::element_text(size = 20),
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
    paste0("G:/My Drive/Missoula_postdoc/TDST_Model/imgs/", 
           file_names[ff],
           "_cam_SD.png"),
    plot = ggplot2::last_plot(),
    # path = file_path,
    # scale = 1,
    width = 8,
    height = 3,
    # units = c("in", "cm", "mm", "px"),
    dpi = 600,
    limitsize = TRUE,
    bg = NULL
  )
  
}

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


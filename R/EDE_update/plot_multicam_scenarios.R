library(dplyr)
library(ggplot2)
tot_animals <- 100

fig_colors <- c("#1B5E20", "#00A8C6", "#FBC02D", "#E65100", "#8E44AD", "#4B6FAD", "#D81B60")

design_names <- c(
  "Random", "Slow_80_bias", "Medium_80_bias", "Fast_80_bias", "Slow_bias", "Medium_bias", "Fast_bias"
)

ncams <- c(10, 25, 50, 75, 100)

D_all_in <- c()
for (dn in design_names) {
  for (nc in ncams) {
    load(
      file = paste0(
        # "G:/My Drive/Missoula_postdoc/TDST_Model/Data/results_",
        "C:/Users/guen.grosklos/Google Drive/Missoula_postdoc/TDST_Model/Data/results_",
        dn,
        "_",
        nc,
        "_cam.RData")
    )
    
    D_all_in <- D_all_in %>% 
      dplyr::bind_rows(
        save_results[[5]] %>% 
          dplyr::mutate(
            Design = dn
          )
      )
  }
}

D_all <- D_all_in %>% 
  dplyr::mutate(
    Model = paste(Model, Covariate),
    Model = dplyr::case_when(
      Model == "PR Non-Covariate" ~ "PD",
      Model == "PR Covariate" ~ "PR",
      Model == "PR Habitat Non-Covariate" ~ "HAM" # or PDH?
    ),
    Design = dplyr::case_when(
      Design %in% c("Slow_80_bias") ~ "80% Slow",
      Design %in% c("Medium_80_bias") ~ "80% Medium",
      Design %in% c("Fast_80_bias") ~ "80% Fast",
      Design %in% c("Slow_bias") ~ "100% Slow",
      Design %in% c("Medium_bias") ~ "100% Medium",
      Design %in% c("Fast_bias") ~ "100% Fast",
      Design %in% c("Random") ~ "Random"
    )
  )
  

D_all$Model <- factor(
  D_all$Model,
  levels = c("PD", "HAM", "PR")
)
D_all$Design <- factor(
  D_all$Design,
  levels = c("Random", "80% Slow", "80% Medium", "80% Fast", "100% Slow", "100% Medium", "100% Fast")
)

# Maybe put in a table - PR requires more data than PD and HAM
D_counts <- D_all %>% 
  dplyr::group_by(Model, cams, Design) %>% 
  dplyr::summarise(
    count = sum(Est > 250, na.rm = T) / 1000,
    .groups = 'drop'
  ) %>% 
  dplyr::mutate(cams = as.character(cams)) 

D_counts$cams <- factor(
  D_counts$cams,
  levels = c("10", "25", "50", "75", "100")
)  

D_counts %>% 
  ggplot2::ggplot(ggplot2::aes(x = Design, y = count, fill = cams)) +
  ggplot2::geom_col(position = "dodge") +
  ggplot2::facet_wrap(~ Model, ncol = 1) +
  ggplot2::labs(
    x = "Sampling Design",
    y = "Proportion of Simulations Omitted"
  ) +
  ggplot2::scale_fill_manual(values= fig_colors, name = "Number of Cameras") +
  ggplot2::theme(text = ggplot2::element_text(size = 16),
                 legend.title = ggplot2::element_blank(),
                 panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(),
                 axis.line = ggplot2::element_line(colour = "black"),
                 panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1)
                 # legend.position = "none",
                 # legend.background = ggplot2::element_blank(),
                 # legend.spacing.y = ggplot2::unit(0, "mm"),
                 # legend.box.background = ggplot2::element_rect(colour = "black")
  )

# D_counts %>% 
#   ggplot2::ggplot(ggplot2::aes(x = Model, y = count, fill = Design)) +
#   ggplot2::geom_col(position = "dodge") +
#   ggplot2::facet_wrap(~ cams, ncol = 1) +
#   ggplot2::labs(
#     x = "Sampling Design",
#     y = "Proportion of Simulations Omitted"
#   ) +
#   ggplot2::scale_fill_manual(values= fig_colors, name = "Number of Cameras") +
#   ggplot2::theme(text = ggplot2::element_text(size = 16),
#                  legend.title = ggplot2::element_blank(),
#                  panel.grid.major = ggplot2::element_blank(),
#                  panel.grid.minor = ggplot2::element_blank(),
#                  panel.background = ggplot2::element_blank(),
#                  axis.line = ggplot2::element_line(colour = "black"),
#                  panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1)
#                  # legend.position = "none",
#                  # legend.background = ggplot2::element_blank(),
#                  # legend.spacing.y = ggplot2::unit(0, "mm"),
#                  # legend.box.background = ggplot2::element_rect(colour = "black")
#   )

ggplot2::ggsave(
  paste0("C:/Users/guen.grosklos/Google Drive/Missoula_postdoc/TDST_Model/imgs/Omitted_counts.png"),
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



D_cam_means <- D_all %>% 
  dplyr::filter(Est < 250) %>%
  dplyr::group_by(cams, Model, Design) %>% #, Speed, Bias) %>% 
  dplyr::summarise(
    mm = mean(Est, na.rm = T), 
    ss = sd(Est, na.rm = T), 
    hh = mean(SD, na.rm = T),
    .groups = 'drop'
  ) %>% 
  dplyr::mutate(
    Speed = dplyr::case_when(
      Design %in% c("80% Slow", "100% Slow") ~ "Slow",
      Design %in% c("80% Medium", "100% Medium") ~ "Medium",
      Design %in% c("80% Fast", "100% Fast") ~ "Fast",
      Design %in% c("Random") ~ "Random"
    ),
    Bias = dplyr::case_when(
      Design %in% c("80% Slow", "80% Medium", "80% Fast") ~ "80%",
      Design %in% c("100% Slow", "100% Medium", "100% Fast") ~ "100%",
      Design %in% c("Random") ~ "Random"
    )
  )

# D_cam_means %>% 
#   # dplyr::filter(Design == "Slow_bias") %>% 
#   dplyr::filter(Model == "HAM") %>% 
#   # dplyr::filter(SampDesign == paste0(file_names[ff], "_cam")) %>%
#   ggplot2::ggplot(ggplot2::aes(x = cams, y = mm, color = Design)) +
#   ggplot2::geom_line() +
#   ggplot2::geom_line(
#     data = D_cam_means %>% 
#       dplyr::filter(Model == "PR Non-Covariate",
#                     Design == "Random"),
#     ggplot2::aes(x = cams, y = mm), 
#     color = "black", linetype = "dashed")

# GOOD ONE
D_cam_means %>% 
  dplyr::filter(Model == "PR") %>%
  ggplot2::ggplot(ggplot2::aes(x = cams, y = hh, color = Design)) + #, shape = Speed)) +
  ggplot2::geom_line(size = 1) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_line(
    data = D_cam_means %>% 
      dplyr::filter(Model == "PD",
                    Design == "Random"),
    ggplot2::aes(x = cams, y = hh), 
    color = "black", linetype = "dashed", size = 1) +
  ggplot2::labs(x = "Number of Cameras",
                y = "Posterior SDs") +
  ggplot2::scale_fill_manual(values= fig_colors) +
  ggplot2::scale_color_manual(values = fig_colors) +
  ggplot2::theme(text = ggplot2::element_text(size = 16),
                 # legend.title = ggplot2::element_blank(),
                 panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(),
                 axis.line = ggplot2::element_line(colour = "black"),
                 panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1)
                 # legend.position = "none",
                 # legend.background = ggplot2::element_blank(),
                 # legend.spacing.y = ggplot2::unit(0, "mm"),
                 # legend.box.background = ggplot2::element_rect(colour = "black")
  )

ggplot2::ggsave(
  paste0("C:/Users/guen.grosklos/Google Drive/Missoula_postdoc/TDST_Model/imgs/PR_multicam_SD.png"),
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

D_cam_means %>% 
  dplyr::filter(Model == "HAM") %>%
  ggplot2::ggplot(ggplot2::aes(x = cams, y = hh, color = Design, shape = Speed)) +
  ggplot2::geom_line(size = 1) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_line(
    data = D_cam_means %>% 
      dplyr::filter(Model == "PD",
                    Design == "Random"),
    ggplot2::aes(x = cams, y = hh), 
    color = "black", linetype = "dashed", size = 1) +
  ggplot2::labs(x = "Number of Cameras",
                y = "Posterior SDs") +
  ggplot2::scale_fill_manual(values= fig_colors) +
  ggplot2::scale_color_manual(values = fig_colors) +
  ggplot2::theme(text = ggplot2::element_text(size = 16),
                 # legend.title = ggplot2::element_blank(),
                 panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(),
                 axis.line = ggplot2::element_line(colour = "black"),
                 panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1)
                 # legend.position = "none",
                 # legend.background = ggplot2::element_blank(),
                 # legend.spacing.y = ggplot2::unit(0, "mm"),
                 # legend.box.background = ggplot2::element_rect(colour = "black")
  )

ggplot2::ggsave(
  "C:/Users/guen.grosklos/Google Drive/Missoula_postdoc/TDST_Model/imgs/HAM_multicam_SD.png",
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

D_cam_means %>% 
  # dplyr::filter(Design == "Slow_bias") %>% 
  dplyr::filter(Model == "HAM") %>%
  # dplyr::filter(Model == "PR") %>%
  ggplot2::ggplot(ggplot2::aes(x = cams, y = mm, color = Design, shape = Speed)) +
  ggplot2::geom_line(size = 1) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_hline(yintercept=tot_animals, linetype="dashed", size = 1) +
  ggplot2::scale_fill_manual(values= fig_colors) +
  ggplot2::scale_color_manual(values = fig_colors) +
  ggplot2::labs(x = "Number of Cameras",
                y = "Posterior Means") +
  ggplot2::theme(text = ggplot2::element_text(size = 16),
                 # legend.title = ggplot2::element_blank(),
                 panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(),
                 axis.line = ggplot2::element_line(colour = "black"),
                 panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1)
                 # legend.position = "none",
                 # legend.background = ggplot2::element_blank(),
                 # legend.spacing.y = ggplot2::unit(0, "mm"),
                 # legend.box.background = ggplot2::element_rect(colour = "black")
  )

ggplot2::ggsave(
  "C:/Users/guen.grosklos/Google Drive/Missoula_postdoc/TDST_Model/imgs/HAM_multicam_Mean.png",
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

D_cam_means %>% 
  dplyr::filter(Model == "PR") %>%
  ggplot2::ggplot(ggplot2::aes(x = cams, y = mm, color = Design)) +
  ggplot2::geom_line(size = 1) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_hline(yintercept=tot_animals, linetype="dashed", size = 1) +
  ggplot2::scale_fill_manual(values= fig_colors) +
  ggplot2::scale_color_manual(values = fig_colors) +
  ggplot2::labs(x = "Number of Cameras",
                y = "Posterior Means") +
  ggplot2::theme(text = ggplot2::element_text(size = 16),
                 # legend.title = ggplot2::element_blank(),
                 panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(),
                 axis.line = ggplot2::element_line(colour = "black"),
                 panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1)
                 # legend.position = "none",
                 # legend.background = ggplot2::element_blank(),
                 # legend.spacing.y = ggplot2::unit(0, "mm"),
                 # legend.box.background = ggplot2::element_rect(colour = "black")
  )

ggplot2::ggsave(
  "C:/Users/guen.grosklos/Google Drive/Missoula_postdoc/TDST_Model/imgs/PR_multicam_Mean.png",
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
###########################################################################

  # ggplot2::ggsave(
  #   paste0("G:/My Drive/Missoula_postdoc/TDST_Model/imgs/", 
  #          file_names[ff],
  #          "_cam_SD.png"),
  #   plot = ggplot2::last_plot(),
  #   # path = file_path,
  #   # scale = 1,
  #   width = 5,
  #   height = 3,
  #   # units = c("in", "cm", "mm", "px"),
  #   dpi = 600,
  #   limitsize = TRUE,
  #   bg = NULL
  # )



###############################################################################


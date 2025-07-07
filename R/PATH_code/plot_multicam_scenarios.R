library(dplyr)
library(ggplot2)
tot_animals <- 100

sim_dir <- "G:/My Drive/Missoula_postdoc/PATH_model/sim_results/"

img_dir <- "G:/My Drive/Missoula_postdoc/PATH_model/imgs"
  
fig_colors <- c("#1B5E20", "#00A8C6", "#FBC02D", "#E65100", "#8E44AD", "#4B6FAD", "#D81B60")

design_names <- c(
  "Random", "Slow_80_bias", "Medium_80_bias", "Fast_80_bias", "Slow_bias", "Medium_bias", "Fast_bias"
)

ncams <- c(10, 25, 50, 75, 100)

D_all_in <- c()
D_all_counts <- c()
for (dn in design_names) {
  for (nc in ncams) {
    load(
      file = paste0(
        sim_dir,
        dn,
        "_",
        nc,
        "_cam.RData")
    )
    
    D_all_in <- D_all_in %>% 
      dplyr::bind_rows(
        save_results[[5]]
      )
    
  }
}

D_all <- D_all_in %>% 
  dplyr::mutate(
    Design = dplyr::case_when(
      cam_design %in% c("Slow_80_bias") ~ "80% Slow",
      cam_design %in% c("Medium_80_bias") ~ "80% Medium",
      cam_design %in% c("Fast_80_bias") ~ "80% Fast",
      cam_design %in% c("Slow_bias") ~ "100% Slow",
      cam_design %in% c("Medium_bias") ~ "100% Medium",
      cam_design %in% c("Fast_bias") ~ "100% Fast",
      cam_design %in% c("Random") ~ "Random"
    )
  )
  
# # Outliers defined as 1.5 times the interquartile range
# # So 250 is about the mean upper bound over all results
# outlier <- D_all  %>% 
#   # dplyr::group_by(Model, Covariate, Design) %>% 
#   dplyr::summarise(
#     q1 = quantile(Est, 0.25, na.rm = T),
#     q3 = quantile(Est, 0.75, na.rm = T),
#     lower_bound = q1 - (1.5 * (q3 - q1)),
#     upper_bound = q3 + (1.5 * (q3 - q1))
#   )
  
D_all$Model <- factor(
  D_all$Model,
  levels = c("IS", "PATH")
)
D_all$Design <- factor(
  D_all$Design,
  levels = c("Random", "80% Slow", "80% Medium", "80% Fast", "100% Slow", "100% Medium", "100% Fast")
)
# 
################################################################################
# Summarize results
D_cam_means <- D_all %>% 
  dplyr::filter(Est < 250) %>%
  dplyr::mutate(MAE = abs(tot_animals - Est)) |> 
  dplyr::group_by(cams, Model, Design) %>% 
  dplyr::summarise(
    Mean_MAE = mean(MAE, na.rm = T),
    Mean_Est = mean(Est, na.rm = T), 
    SD_Est = sd(Est, na.rm = T), 
    Mean_SD = mean(SD, na.rm = T),
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

# Plot all SDs in one plot
D_cam_means_fctr <- D_cam_means %>% 
  dplyr::filter(Model == "PATH") %>%
  dplyr::bind_rows(
    D_cam_means |> 
      dplyr::filter(Model == "IS" & Design == "Random") 
  ) |> 
  dplyr::mutate(Model = paste(Model, Design))

D_cam_means_fctr$Model <- factor(
  D_cam_means_fctr$Model,
  levels = c("IS Random", "PATH Random", 
             "PATH 80% Slow", "PATH 100% Slow", 
             "PATH 80% Medium", "PATH 100% Medium",
             "PATH 80% Fast", "PATH 100% Fast")
)

fctr_pal <- c("#D81B60", "black", "#1B5E20", "#1B5E20", 
              "#00A8C6", "#00A8C6", "#E65100", "#E65100")

names(fctr_pal) <- levels(D_cam_means_fctr$Model)
custom_colors <- ggplot2::scale_colour_manual(name = "Model", values = fctr_pal)

D_cam_means_fctr |> 
  ggplot2::ggplot(ggplot2::aes(x = cams, y = Mean_SD, color = Model, linetype = Model)) + 
  ggplot2::geom_line(size = 1.2) +
  ggplot2::geom_point(size = 2) +
  custom_colors +
  ggplot2::scale_linetype_manual(
    values = c(
      "IS Random" = "solid",
      "PATH Random" = "solid",
      "PATH 80% Slow" = "dashed",
      "PATH 100% Slow" = "dotted",
      "PATH 80% Medium" = "dashed",
      "PATH 100% Medium" = "dotted",
      "PATH 80% Fast" = "dashed",
      "PATH 100% Fast" = "dotted"
    )
  ) +
  ggplot2::labs(x = "Number of Cameras",
                y = "Mean SDs") + 
  ggplot2::scale_x_continuous(breaks = c(10,25,50, 75, 100), 
                              labels = c("10","25","50", "75", "100")) +
  ggplot2::annotate("text", x = 100, y = 105,
                    label = "b",
                    size = 5) +
  ggplot2::theme(
    axis.title=element_text(size = 16),
    axis.text = ggplot2::element_text(size = 16),
    legend.position = "none",# c(0.85, 0.72),
    legend.text = ggplot2::element_text(size = 10),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(colour = "black"),
    panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
    legend.background = ggplot2::element_rect(color = "white"),
    legend.spacing.y = ggplot2::unit(0, "mm"),
    legend.box.background = ggplot2::element_rect(colour = "black")
  )

ggplot2::ggsave(
  paste0(img_dir, "/cam_vs_SD.png"),
  plot = ggplot2::last_plot(),
  # path = file_path,
  # scale = 1,
  width = 7,
  height = 4,
  # units = c("in", "cm", "mm", "px"),
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

# Plot all MAEs
D_cam_means_fctr |> 
  ggplot2::ggplot(ggplot2::aes(x = cams, y = Mean_MAE, color = Model, linetype = Model)) + 
  ggplot2::geom_line(size = 1) +
  ggplot2::geom_point(size = 2) +
  custom_colors +
  ggplot2::scale_linetype_manual(
    values = c(
      "IS Random" = "solid",
      "PATH Random" = "solid",
      "PATH 80% Slow" = "dashed",
      "PATH 100% Slow" = "dotted",
      "PATH 80% Medium" = "dashed",
      "PATH 100% Medium" = "dotted",
      "PATH 80% Fast" = "dashed",
      "PATH 100% Fast" = "dotted"
    )
  ) +
  ggplot2::labs(x = "Number of Cameras",
                y = "Mean Errors") + 
  ggplot2::scale_x_continuous(breaks = c(10,25,50, 75, 100), 
                              labels = c("10","25","50", "75", "100")) +
  ggplot2::annotate("text", x = 100, y = 48,
                    label = "a",
                    size = 5) +
  ggplot2::theme(
    axis.title=element_text(size = 16),
    axis.text = ggplot2::element_text(size = 16),
    legend.position = "none",# c(0.85, 0.72),
    legend.title = ggplot2::element_text(size=13),
    legend.text = ggplot2::element_text(size = 13),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(colour = "black"),
    panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
    legend.background = ggplot2::element_blank(), # ggplot2::element_rect(color = "white"), #
    legend.spacing.y = ggplot2::unit(0, "mm"),
    legend.box.background = ggplot2::element_rect(colour = "black")
  )

ggplot2::ggsave(
  paste0(img_dir, "/cam_vs_MAE.png"),
  plot = ggplot2::last_plot(),
  # path = file_path,
  # scale = 1,
  width = 7,
  height = 4,
  # units = c("in", "cm", "mm", "px"),
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)


################################################################################
# tmp_pal <- c("black", "#1B5E20", "#00A8C6", "#4B6FAD", "#8E44AD")
# # fixed_tag_palette <- setNames(tmp_pal[seq_along(c("IS", "PATH"))], c("IS", "PATH"))
# 
# # Plot sds
# D_cam_means %>% 
#   dplyr::filter(Design == "Random") %>%
#   dplyr::mutate(Model = paste(Model, Design)) |> 
#   ggplot2::ggplot(
#     ggplot2::aes(x = cams, y = Mean_SD, color = Model, linetype = Model)
#   ) + 
#   ggplot2::geom_line(size = 1) +
#   ggplot2::scale_linetype_manual(
#     values = c("IS Random" = "dashed", "PATH Random" = "solid")
#   ) +
#   ggplot2::geom_point(size = 2) +
#   ggplot2::labs(x = "Number of Cameras",
#                 y = "Posterior SDs") +
#   ggplot2::scale_fill_manual(values= tmp_pal) +
#   ggplot2::scale_color_manual(values = tmp_pal) +
#   ggplot2::annotate("text", x = 8, y = 80,# max(D_cam_means$Mean_SD, na.rm = T), 
#                     label = "a",
#                     size = 5) +
#   ggplot2::ggtitle("Random Design") +
#   ggplot2::theme(
#     axis.title=element_text(size = 16),
#     axis.text = ggplot2::element_text(size = 16),
#     legend.position = c(0.85, 0.72),
#     legend.text = ggplot2::element_text(size = 8),
#     panel.grid.major = ggplot2::element_blank(),
#     panel.grid.minor = ggplot2::element_blank(),
#     panel.background = ggplot2::element_blank(),
#     axis.line = ggplot2::element_line(colour = "black"),
#     panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
#     legend.background = ggplot2::element_rect(color = "white"),
#     legend.spacing.y = ggplot2::unit(0, "mm"),
#     legend.box.background = ggplot2::element_rect(colour = "black")
#   )
# 
# D_cam_means %>% 
#   dplyr::filter(Design %in% c("80% Slow", "100% Slow") & Model %in% "PATH") %>%
#   dplyr::bind_rows(
#     D_cam_means |> 
#       dplyr::filter(Model == "IS" & Design == "Random")
#   ) |> 
#   dplyr::mutate(Model = paste(Model, Design)) |> 
#   ggplot2::ggplot(
#     ggplot2::aes(x = cams, y = Mean_SD, color = Model, linetype = Model)
#   ) + 
#   ggplot2::geom_line(size = 1) +
#   ggplot2::scale_linetype_manual(
#     values = c(
#       "IS Random" = "dashed", 
#       "PATH 80% Slow" = "solid", 
#       "PATH 100% Slow" = "solid"
#     )
#   ) +
#   ggplot2::geom_point(size = 2) +
#   ggplot2::labs(x = "Number of Cameras",
#                 y = "Posterior SDs") +
#   ggplot2::scale_fill_manual(values= tmp_pal) +
#   ggplot2::scale_color_manual(values = tmp_pal) +
#   ggplot2::annotate("text", x = 8, y = 80,# max(D_cam_means$Mean_SD, na.rm = T), 
#                     label = "d",
#                     size = 5) +
#   ggplot2::ggtitle("Slow Habitat Design") +
#   ggplot2::theme(
#     axis.title=element_text(size = 16),
#     axis.text = ggplot2::element_text(size = 16),
#     legend.position = c(0.85, 0.72),
#     legend.text = ggplot2::element_text(size = 8),
#     panel.grid.major = ggplot2::element_blank(),
#     panel.grid.minor = ggplot2::element_blank(),
#     panel.background = ggplot2::element_blank(),
#     axis.line = ggplot2::element_line(colour = "black"),
#     panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
#     legend.background = ggplot2::element_rect(color = "white"),
#     legend.spacing.y = ggplot2::unit(0, "mm"),
#     legend.box.background = ggplot2::element_rect(colour = "black")
#   )
# 
# D_cam_means %>% 
#   dplyr::filter(Design %in% c("80% Medium", "100% Medium") & Model %in% "PATH") %>%
#   dplyr::bind_rows(
#     D_cam_means |> 
#       dplyr::filter(Model == "IS" & Design == "Random")
#   ) |> 
#   dplyr::mutate(Model = paste(Model, Design)) |> 
#   ggplot2::ggplot(
#     ggplot2::aes(x = cams, y = Mean_SD, color = Model, linetype = Model)
#   ) + 
#   ggplot2::geom_line(size = 1) +
#   ggplot2::scale_linetype_manual(
#     values = c(
#       "IS Random" = "dashed", 
#       "PATH 80% Medium" = "solid", 
#       "PATH 100% Medium" = "solid"
#     )
#   ) +
#   ggplot2::geom_point(size = 2) +
#   ggplot2::labs(x = "Number of Cameras",
#                 y = "Posterior SDs") +
#   ggplot2::scale_fill_manual(values= tmp_pal) +
#   ggplot2::scale_color_manual(values = tmp_pal) +
#   ggplot2::annotate("text", x = 8, y = 80,# max(D_cam_means$Mean_SD, na.rm = T), 
#                     label = "c",
#                     size = 5) +
#   ggplot2::ggtitle("Medium Habitat Design") +
#   ggplot2::theme(
#     axis.title=element_text(size = 16),
#     axis.text = ggplot2::element_text(size = 16),
#     legend.position = c(0.85, 0.72),
#     legend.text = ggplot2::element_text(size = 8),
#     panel.grid.major = ggplot2::element_blank(),
#     panel.grid.minor = ggplot2::element_blank(),
#     panel.background = ggplot2::element_blank(),
#     axis.line = ggplot2::element_line(colour = "black"),
#     panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
#     legend.background = ggplot2::element_rect(color = "white"),
#     legend.spacing.y = ggplot2::unit(0, "mm"),
#     legend.box.background = ggplot2::element_rect(colour = "black")
#   )
# 
# D_cam_means %>% 
#   dplyr::filter(Design %in% c("80% Fast", "100% Fast") & Model %in% "PATH") %>%
#   dplyr::bind_rows(
#     D_cam_means |> 
#       dplyr::filter(Model == "IS" & Design == "Random")
#   ) |> 
#   dplyr::mutate(Model = paste(Model, Design)) |> 
#   ggplot2::ggplot(
#     ggplot2::aes(x = cams, y = Mean_SD, color = Model, linetype = Model)
#   ) + 
#   ggplot2::geom_line(size = 1) +
#   ggplot2::scale_linetype_manual(
#     values = c(
#       "IS Random" = "dashed", 
#       "PATH 80% Fast" = "solid", 
#       "PATH 100% Fast" = "solid"
#     )
#   ) +
#   ggplot2::geom_point(size = 2) +
#   ggplot2::labs(x = "Number of Cameras",
#                 y = "Posterior SDs") +
#   ggplot2::scale_fill_manual(values= tmp_pal) +
#   ggplot2::scale_color_manual(values = tmp_pal) +
#   ggplot2::annotate("text", x = 8, y = 80,# max(D_cam_means$Mean_SD, na.rm = T), 
#                     label = "b",
#                     size = 5) +
#   ggplot2::ggtitle("Fast Habitat Design") +
#   ggplot2::theme(
#     axis.title=element_text(size = 16),
#     axis.text = ggplot2::element_text(size = 16),
#     legend.position = c(0.85, 0.72),
#     legend.text = ggplot2::element_text(size = 8),
#     panel.grid.major = ggplot2::element_blank(),
#     panel.grid.minor = ggplot2::element_blank(),
#     panel.background = ggplot2::element_blank(),
#     axis.line = ggplot2::element_line(colour = "black"),
#     panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
#     legend.background = ggplot2::element_rect(color = "white"),
#     legend.spacing.y = ggplot2::unit(0, "mm"),
#     legend.box.background = ggplot2::element_rect(colour = "black")
#   )

################################################################################
################################################################################
# # Try boxplots
# D_all_alt <- D_all %>%
#   dplyr::filter(Model == "PATH") |> 
#   dplyr::bind_rows(
#     D_all |>
#       dplyr::filter(Model == "IS" & Design == "Random")
#   ) |>
#   dplyr::mutate(MAE = abs(tot_animals - Est)) |> 
#   dplyr::mutate(Model = paste(Model, Design))
# 
# D_all_alt$Model <- factor(
#   D_all_alt$Model,
#   levels = c("IS Random", "PATH Random", 
#              "PATH 80% Slow", "PATH 100% Slow", 
#              "PATH 80% Medium", "PATH 100% Medium",
#              "PATH 80% Fast", "PATH 100% Fast")
# )
# 
# D_all_alt |> 
#   ggplot2::ggplot(ggplot2::aes(x = as.factor(cams), y = SD, fill = Model)) +
#   ggplot2::geom_boxplot(lwd = 0.5, fatten = .5, outlier.shape = NA) +
#   ggplot2::labs(x = "Number of Cameras",
#                 y = "Posterior SDs") +
#   # ggplot2::scale_fill_manual(values= fig_colors[1:5]) +
#   # ggplot2::scale_color_manual(values = c('grey0', 'grey40', 'grey60')) +
#   # ggplot2::geom_hline(
#   #   yintercept = D_cam_means |>
#   #     dplyr::filter(Model == "IS" & Design == "Random") |>
#   #     dplyr::pull(Mean_SD),
#   #   linetype="dashed",
#   #   size = 0.7
#   # ) |>
#   ggplot2::annotate("text", x = 0.5, y = 32,
#                     label = "b",
#                     size = 5) +
#   ggplot2::theme(text = ggplot2::element_text(size = 16),
#                  legend.title=element_text(size=12),
#                  legend.text=element_text(size=10),
#                  legend.position = c(0.9, 0.77),
#                  panel.grid.major = ggplot2::element_blank(),
#                  panel.grid.minor = ggplot2::element_blank(),
#                  panel.background = ggplot2::element_blank(),
#                  axis.line = ggplot2::element_line(colour = "black"),
#                  panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
#                  legend.background = ggplot2::element_rect(color = "white"),
#                  legend.spacing.y = ggplot2::unit(0, "mm"),
#                  legend.box.background = ggplot2::element_rect(colour = "black"))
# 
# ################################################################################
# # Plot means
# D_cam_means %>% 
#   dplyr::filter(Design == "Random") %>%
#   dplyr::mutate(Model = paste(Model, Design)) |> 
#   ggplot2::ggplot(ggplot2::aes(x = cams, y = Mean_Est, color = Model)) + #, shape = Speed)) +
#   ggplot2::geom_line(size = 1) +
#   ggplot2::geom_point(size = 2) +
#   ggplot2::labs(x = "Number of Cameras",
#                 y = "Posterior Means") +
#   ggplot2::scale_fill_manual(values= fig_colors) +
#   ggplot2::scale_color_manual(values = fig_colors) +
#   ggplot2::annotate("text", x = 8, y = 80,# max(D_cam_means$Mean_SD, na.rm = T), 
#                     label = "b",
#                     size = 5) +
#   ggplot2::ggtitle("Random Design") +
#   ggplot2::theme(
#     axis.title=element_text(size = 16),
#     axis.text = ggplot2::element_text(size = 16),
#     legend.position = c(0.85, 0.72),
#     legend.text = ggplot2::element_text(size = 8),
#     panel.grid.major = ggplot2::element_blank(),
#     panel.grid.minor = ggplot2::element_blank(),
#     panel.background = ggplot2::element_blank(),
#     axis.line = ggplot2::element_line(colour = "black"),
#     panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
#     legend.background = ggplot2::element_rect(color = "white"),
#     legend.spacing.y = ggplot2::unit(0, "mm"),
#     legend.box.background = ggplot2::element_rect(colour = "black")
#   )
# 
# D_cam_means %>% 
#   dplyr::filter(Design %in% c("80% Slow", "100% Slow") & Model %in% "PATH") %>%
#   dplyr::bind_rows(
#     D_cam_means |> 
#       dplyr::filter(Model == "IS" & Design == "Random")
#   ) |> 
#   dplyr::mutate(Model = paste(Model, Design)) |> 
#   ggplot2::ggplot(ggplot2::aes(x = cams, y = Mean_Est, color = Model)) + #, shape = Speed)) +
#   ggplot2::geom_line(size = 1) +
#   ggplot2::geom_point(size = 2) +
#   ggplot2::labs(x = "Number of Cameras",
#                 y = "Posterior Means") +
#   ggplot2::scale_fill_manual(values= fig_colors) +
#   ggplot2::scale_color_manual(values = fig_colors) +
#   ggplot2::annotate("text", x = 8, y = 80,# max(D_cam_means$Mean_SD, na.rm = T), 
#                     label = "b",
#                     size = 5) +
#   ggplot2::ggtitle("Slow Habitat Design") +
#   ggplot2::theme(
#     axis.title=element_text(size = 16),
#     axis.text = ggplot2::element_text(size = 16),
#     legend.position = c(0.85, 0.72),
#     legend.text = ggplot2::element_text(size = 8),
#     panel.grid.major = ggplot2::element_blank(),
#     panel.grid.minor = ggplot2::element_blank(),
#     panel.background = ggplot2::element_blank(),
#     axis.line = ggplot2::element_line(colour = "black"),
#     panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
#     legend.background = ggplot2::element_rect(color = "white"),
#     legend.spacing.y = ggplot2::unit(0, "mm"),
#     legend.box.background = ggplot2::element_rect(colour = "black")
#   )
# 
# D_cam_means %>% 
#   dplyr::filter(Design %in% c("80% Medium", "100% Medium") & Model %in% "PATH") %>%
#   dplyr::bind_rows(
#     D_cam_means |> 
#       dplyr::filter(Model == "IS" & Design == "Random")
#   ) |> 
#   dplyr::mutate(Model = paste(Model, Design)) |> 
#   ggplot2::ggplot(ggplot2::aes(x = cams, y = Mean_Est, color = Model)) + #, shape = Speed)) +
#   ggplot2::geom_line(size = 1) +
#   ggplot2::geom_point(size = 2) +
#   ggplot2::labs(x = "Number of Cameras",
#                 y = "Posterior Means") +
#   ggplot2::scale_fill_manual(values= fig_colors) +
#   ggplot2::scale_color_manual(values = fig_colors) +
#   ggplot2::annotate("text", x = 8, y = 80,# max(D_cam_means$Mean_SD, na.rm = T), 
#                     label = "b",
#                     size = 5) +
#   ggplot2::ggtitle("Medium Habitat Design") +
#   ggplot2::theme(
#     axis.title=element_text(size = 16),
#     axis.text = ggplot2::element_text(size = 16),
#     legend.position = c(0.85, 0.72),
#     legend.text = ggplot2::element_text(size = 8),
#     panel.grid.major = ggplot2::element_blank(),
#     panel.grid.minor = ggplot2::element_blank(),
#     panel.background = ggplot2::element_blank(),
#     axis.line = ggplot2::element_line(colour = "black"),
#     panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
#     legend.background = ggplot2::element_rect(color = "white"),
#     legend.spacing.y = ggplot2::unit(0, "mm"),
#     legend.box.background = ggplot2::element_rect(colour = "black")
#   )
# 
# D_cam_means %>% 
#   dplyr::filter(Design %in% c("80% Fast", "100% Fast") & Model %in% "PATH") %>%
#   dplyr::bind_rows(
#     D_cam_means |> 
#       dplyr::filter(Model == "IS" & Design == "Random")
#   ) |> 
#   dplyr::mutate(Model = paste(Model, Design)) |> 
#   ggplot2::ggplot(ggplot2::aes(x = cams, y = Mean_Est, color = Model)) + #, shape = Speed)) +
#   ggplot2::geom_line(size = 1) +
#   ggplot2::geom_point(size = 2) +
#   ggplot2::labs(x = "Number of Cameras",
#                 y = "Posterior Means") +
#   ggplot2::scale_fill_manual(values= fig_colors) +
#   ggplot2::scale_color_manual(values = fig_colors) +
#   ggplot2::annotate("text", x = 8, y = 80,# max(D_cam_means$Mean_SD, na.rm = T), 
#                     label = "b",
#                     size = 5) +
#   ggplot2::ggtitle("Fast Habitat Design") +
#   ggplot2::theme(
#     axis.title=element_text(size = 16),
#     axis.text = ggplot2::element_text(size = 16),
#     legend.position = c(0.85, 0.72),
#     legend.text = ggplot2::element_text(size = 8),
#     panel.grid.major = ggplot2::element_blank(),
#     panel.grid.minor = ggplot2::element_blank(),
#     panel.background = ggplot2::element_blank(),
#     axis.line = ggplot2::element_line(colour = "black"),
#     panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
#     legend.background = ggplot2::element_rect(color = "white"),
#     legend.spacing.y = ggplot2::unit(0, "mm"),
#     legend.box.background = ggplot2::element_rect(colour = "black")
#   )
# 
# 
################################################################################
# # Maybe put in a table - PR requires more data than PD and HAM
# D_omit <- D_all %>% 
#   dplyr::group_by(Model, cams, Design) %>% 
#   dplyr::summarise(
#     count = sum(Est > 250, na.rm = T) / 1000,
#     .groups = 'drop'
#   ) %>% 
#   dplyr::mutate(cams = as.character(cams)) 
# 
# D_omit$cams <- factor(
#   D_omit$cams,
#   levels = c("10", "25", "50", "75", "100")
# )  
# 
# D_omit %>% 
#   ggplot2::ggplot(ggplot2::aes(x = Design, y = count, fill = cams)) +
#   ggplot2::geom_col(position = "dodge") +
#   ggplot2::facet_wrap(~ Model, ncol = 1) +
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
# 
# D_omit %>%
#   # dplyr::filter(!(Design %in% c("100% Fast", "100% Medium", "100% Slow"))) %>% 
#   ggplot2::ggplot(ggplot2::aes(x = Design, y = count, fill = Model)) +
#   ggplot2::geom_col(position = "dodge", color = 'black') +
#   ggplot2::facet_wrap(~ cams, ncol = 1) +
#   ggplot2::labs(
#     x = "Sampling Design",
#     y = "Proportion of Simulations Omitted"
#   ) +
#   ggplot2::scale_fill_manual(values= fig_colors, name = "Number of Cameras") +
#   ggplot2::scale_color_manual(values= fig_colors, name = "Number of Cameras") +
#   # ggplot2::ylim(0,0.2) +
#   ggplot2::scale_y_continuous(limits = c(0,0.2), expand = c(0,0)) +
#   ggplot2::theme(axis.title=element_text(size = 16),
#                  axis.text = ggplot2::element_text(size = 8),
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

# ggplot2::ggsave(
#   paste0(img_dir, "Omitted_counts.png"),
#   plot = ggplot2::last_plot(),
#   # path = file_path,
#   # scale = 1,
#   width = 7,
#   height = 5,
#   # units = c("in", "cm", "mm", "px"),
#   dpi = 600,
#   limitsize = TRUE,
#   bg = NULL
# )


library(dplyr)
library(ggplot2)
tot_animals <- 100

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
    
    D_all_counts <- D_all_counts %>% 
      dplyr::bind_rows(
        save_results[[4]]$count_data %>% 
          purrr::map2_dfr(
            .y = seq_along(save_results[[4]]$count_data),
            .f = function(tbl, idx) {
              tbl %>% mutate(iteration = idx)
            }
          ) %>% 
          dplyr::mutate(
            Design = dn,
            cams = nc
          )
      )
  }
}

D_counts_mean <- D_all_counts %>% 
  dplyr::group_by(iteration, Design, cams) %>% 
  dplyr::summarise(tot_count = sum(count),
                   .groups = 'drop')

D_all <- D_all_in %>% 
  dplyr::left_join(
    D_counts_mean,
    by = join_by(iteration, cams, Design)
  ) %>% 
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
  levels = c("PD", "HAM", "PR")
)
D_all$Design <- factor(
  D_all$Design,
  levels = c("Random", "80% Slow", "80% Medium", "80% Fast", "100% Slow", "100% Medium", "100% Fast")
)

# Maybe put in a table - PR requires more data than PD and HAM
D_omit <- D_all %>% 
  dplyr::group_by(Model, cams, Design) %>% 
  dplyr::summarise(
    count = sum(Est > 250, na.rm = T) / 1000,
    .groups = 'drop'
  ) %>% 
  dplyr::mutate(cams = as.character(cams)) 

D_omit$cams <- factor(
  D_omit$cams,
  levels = c("10", "25", "50", "75", "100")
)  

D_omit %>% 
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

D_omit %>%
  # dplyr::filter(!(Design %in% c("100% Fast", "100% Medium", "100% Slow"))) %>% 
  ggplot2::ggplot(ggplot2::aes(x = Design, y = count, fill = Model)) +
  ggplot2::geom_col(position = "dodge", color = 'black') +
  ggplot2::facet_wrap(~ cams, ncol = 1) +
  ggplot2::labs(
    x = "Sampling Design",
    y = "Proportion of Simulations Omitted"
  ) +
  ggplot2::scale_fill_manual(values= fig_colors, name = "Number of Cameras") +
  ggplot2::scale_color_manual(values= fig_colors, name = "Number of Cameras") +
  # ggplot2::ylim(0,0.2) +
  ggplot2::scale_y_continuous(limits = c(0,0.2), expand = c(0,0)) +
  ggplot2::theme(axis.title=element_text(size = 16),
                 axis.text = ggplot2::element_text(size = 8),
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

ggplot2::ggsave(
  paste0("C:/Users/guen.grosklos/Google Drive/Missoula_postdoc/TDST_Model/imgs/Omitted_counts.png"),
  plot = ggplot2::last_plot(),
  # path = file_path,
  # scale = 1,
  width = 7,
  height = 5,
  # units = c("in", "cm", "mm", "px"),
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

# # Plot counts vs precision
# nested_data <- D_all %>%
#   dplyr::filter(Est < 250, !is.na(Est)) %>% 
#   group_by(cams, Model, Design) %>%
#   tidyr::nest() %>%
#   mutate(
#     model = purrr::map(data, ~lm(tot_count ~ SD, data = .x)),
#     pred_data = purrr::map(data, ~tibble(
#       SD = seq(min(.x$SD), max(.x$SD), length.out = 100)
#     )),
#     predictions = purrr::map2(model, pred_data, ~predict(.x, newdata = .y, 
#                                                   interval = "confidence")),
#     pred_data = purrr::map2(pred_data, predictions, ~bind_cols(.x, as_tibble(.y)))
#   ) %>%
#   select(cams, Model, Design, pred_data) %>%
#   tidyr::unnest(pred_data)
# 
# D_all %>%
#   dplyr::filter(Est < 250, !is.na(Est)) %>% 
#   dplyr::filter(cams == 10) %>% 
#   ggplot() +
#   # Add the actual data points
#   geom_point(aes(x = SD, y = tot_count, color = Model)) +
#   # Add prediction lines
#   geom_line(data = nested_data %>% 
#               dplyr::filter(cams == 10), 
#             aes(x = SD, y = fit, color = Model)) +
#   # Add confidence ribbons
#   geom_ribbon(data = nested_data %>% 
#                 dplyr::filter(cams == 10), 
#               aes(x = SD, y = fit, ymin = lwr, ymax = upr, 
#                   fill = Model), alpha = 0.2) +
#   # Facet by Design if there are multiple designs
#   facet_wrap(~ Design, scales = "free") +
#   # Add labels and theme
#   labs(
#     title = "Prediction of tot_count vs SD",
#     subtitle = "Grouped by cams, Model, and Design",
#     x = "Standard Deviation (SD)",
#     y = "Total Count",
#     color = "Model",
#     fill = "Model",
#     linetype = "Cams",
#     shape = "Cams"
#   ) +
#   theme_minimal()


D_cam_means <- D_all %>% 
  dplyr::filter(Est < 250) %>%
  dplyr::group_by(cams, Model, Design) %>% #, Speed, Bias) %>% 
  dplyr::summarise(
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

# D_cam_means %>% 
#   # dplyr::filter(Design == "Slow_bias") %>% 
#   dplyr::filter(Model == "HAM") %>% 
#   # dplyr::filter(SampDesign == paste0(file_names[ff], "_cam")) %>%
#   ggplot2::ggplot(ggplot2::aes(x = cams, y = Mean_Est, color = Design)) +
#   ggplot2::geom_line() +
#   ggplot2::geom_line(
#     data = D_cam_means %>% 
#       dplyr::filter(Model == "PR Non-Covariate",
#                     Design == "Random"),
#     ggplot2::aes(x = cams, y = Mean_Est), 
#     color = "black", linetype = "dashed")

# GOOD ONE
D_cam_means %>% 
  dplyr::filter(Model == "PR") %>%
  ggplot2::ggplot(ggplot2::aes(x = cams, y = Mean_SD, color = Design)) + #, shape = Speed)) +
  ggplot2::geom_line(size = 1) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_line(
    data = D_cam_means %>% 
      dplyr::filter(Model == "PD",
                    Design == "Random"),
    ggplot2::aes(x = cams, y = Mean_SD), 
    color = "black", size = 1) +
  ggplot2::geom_point(
    data = D_cam_means %>% 
      dplyr::filter(Model == "PD",
                    Design == "Random"),
    ggplot2::aes(x = cams, y = Mean_SD), 
    color = "black", shape = 4, size = 2, stroke = 1.5
  ) +
  ggplot2::labs(x = "Number of Cameras",
                y = "Posterior SDs") +
  ggplot2::scale_fill_manual(values= fig_colors) +
  ggplot2::scale_color_manual(values = fig_colors) +
  ggplot2::annotate("text", x = 8, y = 80,# max(D_cam_means$Mean_SD, na.rm = T), 
                               label = "b",
                               size = 5) +
  ggplot2::theme(
    axis.title=element_text(size = 16),
    axis.text = ggplot2::element_text(size = 16),
    legend.position = c(0.85, 0.72),
    legend.text = ggplot2::element_text(size = 8),
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
  ggplot2::ggplot(ggplot2::aes(x = cams, y = Mean_SD, color = Design)) + #, shape = Speed)) +
  ggplot2::geom_line(size = 1) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_line(
    data = D_cam_means %>% 
      dplyr::filter(Model == "PD",
                    Design == "Random"),
    ggplot2::aes(x = cams, y = Mean_SD), 
    color = "black", size = 1) +
  ggplot2::geom_point(
    data = D_cam_means %>% 
      dplyr::filter(Model == "PD",
                    Design == "Random"),
    ggplot2::aes(x = cams, y = Mean_SD), 
    color = "black", shape = 4, size = 2, stroke = 1.5
  ) +
  ggplot2::labs(x = "Number of Cameras",
                y = "Posterior SDs") +
  ggplot2::scale_fill_manual(values= fig_colors) +
  ggplot2::scale_color_manual(values = fig_colors) +
  ggplot2::annotate("text", x = 8, y = 105,# max(D_cam_means$Mean_SD, na.rm = T), 
                    label = "d",
                    size = 5) +
  ggplot2::guides(color = ggplot2::guide_legend(ncol = 2)) +
  ggplot2::theme(    
    axis.title=element_text(size = 16),
    axis.text = ggplot2::element_text(size = 16),
    legend.position = c(0.72, 0.72),
    legend.text = ggplot2::element_text(size = 8),
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
  ggplot2::ggplot(ggplot2::aes(x = cams, y = Mean_Est, color = Design)) + #, shape = Speed)) +
  ggplot2::geom_line(size = 1) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_line(
    data = D_cam_means %>% 
      dplyr::filter(Model == "PD",
                    Design == "Random"),
    ggplot2::aes(x = cams, y = Mean_Est), 
    color = "black", size = 1) +
  ggplot2::geom_point(
    data = D_cam_means %>% 
      dplyr::filter(Model == "PD",
                    Design == "Random"),
    ggplot2::aes(x = cams, y = Mean_Est), 
    color = "black", shape = 4, size = 2, stroke = 1.5
  ) +
  ggplot2::geom_hline(yintercept=tot_animals, linetype="dashed", size = 1) +
  ggplot2::scale_fill_manual(values= fig_colors) +
  ggplot2::scale_color_manual(values = fig_colors) +
  ggplot2::labs(x = "Number of Cameras",
                y = paste0("HAM \n", "Posterior Means")) +
  ggplot2::annotate("text", x = 8, y = 130,# max(D_cam_means$Mean_SD, na.rm = T), 
                    label = "c",
                    size = 5) +
  ggplot2::theme(    
    axis.title=element_text(size = 16),
    axis.text = ggplot2::element_text(size = 16),
    legend.title = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(colour = "black"),
    panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
    legend.position = "none",
    legend.background = ggplot2::element_blank(),
    legend.spacing.y = ggplot2::unit(0, "mm"),
    legend.box.background = ggplot2::element_rect(colour = "black")
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
  ggplot2::ggplot(ggplot2::aes(x = cams, y = Mean_Est, color = Design)) +
  ggplot2::geom_line(size = 1) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_line(
    data = D_cam_means %>% 
      dplyr::filter(Model == "PD",
                    Design == "Random"),
    ggplot2::aes(x = cams, y = Mean_Est), 
    color = "black", size = 1) +
  ggplot2::geom_point(
    data = D_cam_means %>% 
      dplyr::filter(Model == "PD",
                    Design == "Random"),
    ggplot2::aes(x = cams, y = Mean_Est), 
    color = "black", shape = 4, size = 2, stroke = 1.5
  ) +
  ggplot2::geom_hline(yintercept=tot_animals, linetype="dashed", size = 1) +
  ggplot2::scale_fill_manual(values= fig_colors) +
  ggplot2::scale_color_manual(values = fig_colors) +
  ggplot2::labs(x = "Number of Cameras",
                y = paste0("PR \n", "Posterior Means")) +
  ggplot2::annotate("text", x = 8, y = 102,# max(D_cam_means$Mean_SD, na.rm = T), 
                    label = "a",
                    size = 5) +
  ggplot2::theme(    
    axis.title=element_text(size = 16),
    axis.text = ggplot2::element_text(size = 16),
    legend.title = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(colour = "black"),
    panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
    legend.position = "none",
    legend.background = ggplot2::element_blank(),
    legend.spacing.y = ggplot2::unit(0, "mm"),
    legend.box.background = ggplot2::element_rect(colour = "black")
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




################################
## Plot ABM simulations
################################
################################################################################
#' @export
#'
plot_ABM <- function(study_design,
                     cam_design,
                     cam_locs,
                     animalxy.all
                     ) {
  bounds <- unlist(study_design$bounds)
  ncam <- cam_design$ncam

  # Sample area border
  b.df <- data.frame(c(bounds, study_design$dx))

  # Camera locations
  dt.cam <- data.frame(
    group = rep(1:ncam, each = 3),
    polygon.x = unlist(cam_locs$x),
    polygon.y = unlist(cam_locs$y)
  )


  g <- ggplot() +
    # geom_tile(aes(X-.5*dx, Y-.5*dy, fill = Speed)) +
    scale_fill_manual(values = c("grey20", "grey50", "grey80")) +
    geom_path(
      data = animalxy.all,
      aes(X, Y, col = as.factor(Animal_ID))
    ) +
    geom_rect(
      data = b.df,
      aes(xmin = bounds[1],
          xmax = bounds[2],
          ymin = bounds[1],
          ymax = bounds[2]),
      color = "black",
      fill = NA
    ) +
    geom_polygon(
      data = dt.cam,
      aes(x = polygon.x,
          y = polygon.y,
          group = group),
      color = "black",
      fill = NA
    ) +
    labs(x = "X", y = "Y", fill = "Landscape Type") +
    guides(
      colour = "none",
      alpha = "none"
    ) +
    theme(panel.background = element_blank())
  g + coord_fixed()
}

################################################################################
#' @export
#'
plot_ABM_single <- function(study_design,
                            cam_design,
                            cam_locs,
                            animalxy.all,
                            ID_in = 1
                            ) {
  # Plot movement of one animal only
  bounds <- unlist(study_design$bounds)

  # Sample area border
  b.df <- data.frame(c(bounds, study_design$dx))

  g <- ggplot() +
    scale_fill_manual(values = c("grey20", "grey50", "grey80")) +
    geom_path(
      data = animalxy.all |>
        dplyr::filter(Animal_ID == ID_in),
      aes(X, Y, col = as.factor(Animal_ID))
    ) +
    labs(x = "X", y = "Y", fill = "Landscape Type") +
    guides(
      colour = "none",
      alpha = "none"
    ) +
    theme(panel.background = element_blank())
  g + coord_fixed()
}

################################################################################
#' @export
#'
plot_heat_custom_cells <- function(animal_locs, bounds, q) {
  # NEED TO UPDATE
  # q is the total number of grid cells (must be perfect square)

  n.all <- raster(,
    nrows = q^0.5, ncols = q^0.5, xmn = 0,
    xmx = q^0.5, ymn = 0, ymx = q^0.5, crs = NA
  )
  n.all[] <- 0

  n.temp <- matrix(0, nrow = q^0.5, ncol = q^0.5)

  for (xx in 1:nrow(animal_locs)) {
    x.round <- ceiling(animal_locs$x[xx] * (q^0.5 / max(bounds)))
    y.round <- ceiling(animal_locs$y[xx] * (q^0.5 / max(bounds)))

    # Transpose for converting matrix to raster
    n.temp[x.round, y.round] <- n.temp[x.round, y.round] + 1
  }

  n_plot <- lattice::levelplot(n.temp,
    layers = 1,
    main = list(expression("Animal space-use"), cex = 1),
    cuts = 254,
    margin = FALSE,
    scales = list(draw = FALSE),
    xlab = "x",
    ylab = "y",
    subset(n.temp > 0),
    col.regions = colorRampPalette(rev(RColorBrewer::brewer.pal(3, "Blues")), bias = 3)
  )

  gridExtra::grid.arrange(n_plot, ncol = 1, nrow = 1)
}

################################################################################
#' @export
#'
plot_lscape_types <- function(subset_lscape) {
  # NEED TO UPDATE
  # Sample area border
  b.df <- data.frame(c(bounds, dx))

  # lscape type locations
  lscape_locs <- subset_lscape |>
    dplyr::mutate(
      xmin = X - dx,
      xmax = X,
      ymin = Y - dx,
      ymax = Y,
      a = (ymax - ymin) * (xmax - xmin)
    )

  # Triangle viewshed cameras
  dt.triangle <- data.frame(
    group = rep(1:ncam, each = 3),
    polygon.x = tri_cam_samps$x,
    polygon.y = tri_cam_samps$y
  )


  g <- ggplot(lscape_speeds) +
    # geom_tile(aes(X-.5*dx, Y-.5*dy, fill = Speed)) +
    scale_fill_manual(values = c("grey20", "grey50", "grey80")) +
    geom_path(
      data = animalxy.all,
      aes(x, y, col = as.factor(Animal_ID))
    ) +
    geom_rect(
      data = b.df,
      aes(xmin = bounds[1], xmax = bounds[2], ymin = bounds[1], ymax = bounds[2]),
      color = "black",
      fill = NA
    ) +
    geom_rect(
      data = lscape_locs,
      aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax),
      color = "black",
      fill = NA
    ) +
    geom_polygon(
      data = dt.triangle,
      aes(x = polygon.x, y = polygon.y, group = group),
      color = "black",
      fill = NA
    ) +
    labs(x = "X", y = "Y", fill = "Landscape Type") +
    guides(
      colour = "none",
      alpha = "none"
    ) +
    theme(panel.background = element_blank())
  g + coord_fixed()
}

################################################################################
#' @export
#'
plot_space_use <- function(study_design,
                           animalxy.all
                           ) {
  q <- study_design$q
  bounds <- unlist(study_design$bounds)

  u.abm.all <- matrix(0, nrow = q^0.5, ncol = q^0.5)

  for (xx in 1:nrow(animalxy.all)) {
    x.round <- ceiling(animalxy.all$X[xx] * (q^0.5 / max(bounds)))
    y.round <- ceiling(animalxy.all$Y[xx] * (q^0.5 / max(bounds)))

    # Transpose for converting matrix to raster
    u.abm.all[x.round, y.round] <- u.abm.all[x.round, y.round] + 1
  }

  u_plot <- levelplot(u.abm.all,
                      layers = 1,
                      main = list(expression("Animal space-use"), cex = 1),
                      cuts = 254,
                      margin = FALSE,
                      scales = list(draw = FALSE),
                      xlab = "x",
                      ylab = "y",
                      col.regions = colorRampPalette(rev(brewer.pal(11, "Spectral")), bias = 1))

  grid.arrange(u_plot, ncol = 1, nrow = 1)
}

################################################################################
#' @export
#'
plot_space_use_zoom <- function(animalxy.all, bounds, q) {
  # NEED TO UPDATE
  u.abm.all <- matrix(0, nrow = q^0.5, ncol = q^0.5)

  for (xx in 1:nrow(animalxy.all)) {
    x.round <- ceiling(animalxy.all$x[xx] * (q^0.5 / max(bounds)))
    y.round <- ceiling(animalxy.all$y[xx] * (q^0.5 / max(bounds)))

    # Transpose for converting matrix to raster
    u.abm.all[x.round, y.round] <- u.abm.all[x.round, y.round] + 1
  }

  new_bounds <- which(u.abm.all != 0, arr.ind = T)

  u.abm.zoom <- u.abm.all[
    (min(new_bounds[, 1]) - 1):(max(new_bounds[, 1]) + 1),
    (min(new_bounds[, 2]) - 1):(max(new_bounds[, 2]) + 1)
  ]

  u_plot <- levelplot(u.abm.zoom,
    layers = 1,
    main = list(expression("Animal space-use"), cex = 1),
    cuts = 254,
    margin = FALSE,
    scales = list(
      x = list(
        at = seq(1, (max(new_bounds[, 1]) - min(new_bounds[, 1]) + 2), length = 5),
        labels = as.character(seq(round(min(animalxy.all$x)) - 1, round(max(animalxy.all$x)) + 1, length = 5))
      ),
      y = list(
        at = seq(1, (max(new_bounds[, 2]) - min(new_bounds[, 2]) + 2), length = 5),
        labels = as.character(seq(round(min(animalxy.all$y)) - 1, round(max(animalxy.all$y)) + 1, length = 5))
      )
    ),
    xlab = "x",
    ylab = "y",
    col.regions = colorRampPalette(rev(brewer.pal(11, "Spectral")), bias = 1)
  )

  grid.arrange(u_plot, ncol = 1, nrow = 1)
}


################################
## Plot data density plots
################################
# # Density plot for count data
################################################################################
#' @export
#'
plot_count_data <- function(count_data_in,
                            fill = "Speed") {
  count_data_in$Speed <- factor(count_data_in$Speed, c("Slow", "Medium", "Fast"))

  ggplot(count_data_in, aes(x = count, linetype = .data[[fill]])) +
    geom_density(fill = "grey", position = "identity", alpha = 0.4, adjust = 2, size = 1.5) +
    labs(x = "Count", y = "Frequency", fill = "Landscape Type") +
    scale_x_continuous(
      limits = c(0, max(count_data_in$count)),
      expand = expansion(mult = 0, add = 0)
    ) +
    scale_y_continuous(
      limits = c(0, NA),
      expand = expansion(mult = 0, add = c(0, 0.01))
    ) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    guides(linetype = guide_legend(title = "Habitat Type")) +
    theme(
      text = element_text(size = 20),
      legend.title = element_text(size = 14),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      legend.position = c(0.85, 0.75),
      legend.background = element_blank(),
      legend.spacing.y = unit(0, "mm")
    )

  # ggplot(count_data_in, aes(x=count, fill = .data[[fill]])) +
  #   geom_histogram(position = "identity", alpha = 0.4) +
  #   labs(x = "Count", y = "Frequency", fill="Landscape Type") +
  #   theme(axis.text=element_text(size=20),
  #         axis.title=element_text(size=22),
  #         legend.text=element_text(size=20),
  #         legend.title = element_text(size=20),
  #         panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(),
  #         panel.background = element_blank(),
  #         panel.border = element_rect(colour = "black", fill=NA, size=1))
}

# # Density plot for encounter data
################################################################################
#' @export
#'
plot_encounter_data <- function(encounter_data_in,
                                fill = "Speed") {
  encounter_data_in$Speed <- factor(encounter_data_in$Speed, c("Slow", "Medium", "Fast"))

  ggplot(encounter_data_in, aes(x = encounter, linetype = .data[[fill]])) +
    geom_density(fill = "grey", position = "identity", alpha = 0.4, adjust = 2, size = 1.5) +
    labs(x = "Number of Encounters", y = "Frequency", fill = "Landscape Type") +
    scale_x_continuous(
      limits = c(0, max(encounter_data_in$encounter)),
      expand = expansion(mult = 0, add = c(0, 0.3))
    ) +
    scale_y_continuous(
      limits = c(0, NA),
      expand = expansion(mult = 0, add = c(0, 0.01))
    ) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    guides(linetype = guide_legend(title = "Habitat Type")) +
    theme(
      text = element_text(size = 20),
      legend.title = element_text(size = 14),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      legend.position = c(0.85, 0.75),
      legend.background = element_blank(),
      legend.spacing.y = unit(0, "mm")
    )

  # ggplot(encounter_data_in, aes(x=encounter, fill = .data[[fill]])) +
  #   geom_histogram(position = "identity", alpha = 0.4) +
  #   labs(x = "Count", y = "Frequency", fill="Landscape Type") +
  #   theme(axis.text=element_text(size=20),
  #         axis.title=element_text(size=22),
  #         legend.text=element_text(size=20),
  #         legend.title = element_text(size=20),
  #         panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(),
  #         panel.background = element_blank(),
  #         panel.border = element_rect(colour = "black", fill=NA, size=1))
}

# # Density plot for staying time data
################################################################################
#' @export
#'
plot_staytime_data <- function(stay_time_raw_in,
                               fill = "Speed") {
  stay_time_raw_in$Speed <- factor(stay_time_raw_in$Speed, c("Slow", "Medium", "Fast"))

  ggplot(stay_time_raw_in, aes(x = t_stay, linetype = .data[[fill]])) +
    geom_density(fill = "grey", position = "identity", alpha = 0.4, adjust = 2, size = 1.5) +
    labs(x = "Staying Time", y = "Frequency", fill = "Landscape Type") +
    scale_x_continuous(
      limits = c(0, max(stay_time_raw_in$t_stay)),
      expand = expansion(mult = 0, add = 0)
    ) +
    scale_y_continuous(
      limits = c(0, NA),
      expand = expansion(mult = 0, add = c(0, 0.01))
    ) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    guides(linetype = guide_legend(title = "Habitat Type")) +
    theme(
      text = element_text(size = 20),
      legend.title = element_text(size = 14),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      legend.position = c(0.85, 0.75),
      legend.background = element_blank(),
      legend.spacing.y = unit(0, "mm")
    )
}

# # Density plot for TTE data
################################################################################
#' @export
#'
plot_TTE_data <- function(TTE_data_raw_in,
                          fill = "Speed") {
  TTE_data_raw_in$Speed <- factor(TTE_data_raw_in$Speed, c("Slow", "Medium", "Fast"))

  ggplot(TTE_data_raw_in, aes(x = TTE, linetype = .data[[fill]])) +
    geom_density(fill = "grey", position = "identity", alpha = 0.4, adjust = 2, size = 1.5) +
    labs(x = "Time to Encounter", y = "Frequency", fill = "Landscape Type") +
    scale_x_continuous(
      limits = c(0, max(TTE_data_raw_in$TTE)),
      expand = expansion(mult = 0, add = 0)
    ) +
    scale_y_continuous(
      limits = c(0, NA),
      expand = expansion(mult = 0, add = c(0, 0.01))
    ) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    guides(linetype = guide_legend(title = "Habitat Type")) +
    theme(
      text = element_text(size = 20),
      legend.title = element_text(size = 14),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      legend.position = c(0.85, 0.75),
      legend.background = element_blank(),
      legend.spacing.y = unit(0, "mm")
    )
}

################################
## Plot single-run results
################################

################################################################################
#' @export
#'
plot_onerun_results <- function(D.all) {
  ggplot(D.all, aes(x = Model)) +
    geom_point(aes(y = Est)) +
    geom_hline(yintercept = nind, linetype = "dashed") +
    geom_errorbar(aes(ymin = Est - SD, ymax = Est + SD), width = 0.2)
}


################################
## Plot multi-run results
################################

################################################################################
#' @export
#'
plot_multirun_means <- function(study_design,
                                D.all) {

  ggplot(D.all, aes(x = Model, y = Est, fill = Covariate)) +
    geom_boxplot(position = position_dodge2(preserve = "single"),
                 outlier.shape = NA) +
    geom_hline(yintercept = study_design$tot_animals, linetype = "dashed", size = 1) +
    labs(
      x = "Model",
      y = "Mean Abundance"
    ) +
    coord_cartesian(ylim = quantile(D.all$Est, c(0.01, 0.99), na.rm = T)) +
    # scale_y_continuous(limits = c(lower, upper)) +
    theme(
      text = element_text(size = 20),
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      legend.position = c(0.17, 0.84),
      legend.background = element_blank(),
      legend.spacing.y = unit(0, "mm"),
      legend.box.background = element_rect(colour = "black")
    )
}

################################################################################
#' @export
#'
plot_grouped_multirun_means <- function(Unused_cov = "None",
                                        Filter_model = "None",
                                        Cov = "None") {
  # NEED TO UPDATE

  # Filter with covariate labels not used
  D.all |>
    dplyr::filter(Model != Filter_model) |>
    dplyr::filter(Covariate != Unused_cov) |>
    ggplot(aes(x = Model, y = Est, fill = Run)) +
    geom_boxplot(lwd = 0.5, fatten = .5, outlier.size = 1) +
    geom_hline(yintercept = nind, linetype = "dashed", size = 0.7) +
    labs(
      x = "Model",
      y = if_else(Unused_cov == "None",
        "Posterior Means",
        paste0(Cov, " Models \n", "Posterior Means")
      )
    ) +
    theme(
      text = element_text(size = 16),
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      legend.position = "none",
      legend.background = element_blank(),
      legend.spacing.y = unit(0, "mm"),
      legend.box.background = element_rect(colour = "black")
    )
}

################################################################################
#' @export
#'
plot_multirun_sds <- function(D.all) {
  ggplot(D.all, aes(x = Model, y = SD, fill = Covariate)) +
    geom_boxplot(position = position_dodge2(preserve = "single"),
                 outlier.shape = NA) +
    labs(
      x = "Model",
      y = "SD"
    ) +
    coord_cartesian(ylim = quantile(D.all$SD, c(0.01, 0.99), na.rm = T)) +
    theme(
      text = element_text(size = 20),
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      legend.position = c(0.17, 0.84),
      legend.background = element_blank(),
      legend.spacing.y = unit(0, "mm"),
      legend.box.background = element_rect(colour = "black")
    )
}

################################################################################
#' @export
#'
plot_multirun_CV <- function(D.all) {
  ggplot(D.all, aes(x = Model, y = SD / Est, fill = Covariate)) +
    geom_boxplot(position = position_dodge2(preserve = "single"),
                 outlier.shape = NA) +
    labs(
      x = "Model",
      y = "CV"
    ) +
    coord_cartesian(ylim = quantile(D.all$SD / D.all$Est, c(0.01, 0.99), na.rm = T)) +
    theme(
      text = element_text(size = 20),
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      legend.position = c(0.17, 0.84),
      legend.background = element_blank(),
      legend.spacing.y = unit(0, "mm"),
      legend.box.background = element_rect(colour = "black")
    )
}

################################################################################
#' @export
#'
plot_grouped_multirun_CV <- function(Unused_cov = "None",
                                     Filter_model = "None",
                                     Title) {
  # NEED TO UPDATE
  D.all |>
    dplyr::filter(Model != Filter_model) |>
    dplyr::filter(Covariate != Unused_cov) |>
    ggplot(aes(x = Model, y = SD, fill = Run)) +
    geom_boxplot(lwd = 0.5, fatten = .5, outlier.size = 1) +
    labs(
      x = "Model",
      y = "Posterior SDs"
    ) +
    guides(fill = guide_legend(title = Title)) +
    ylim(c(0, max(D.all$SD, na.rm = T) + 0.1)) +
    scale_fill_manual(values = if (sim_num == 4) {
      c("grey40", "grey60", "grey80")
    } else {
      c("grey20", "grey40", "grey60", "grey80")
    }) +
    theme(
      text = element_text(size = 16),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      legend.position = if (Unused_cov == "None") {
        c(0.12, 0.77)
      } else {
        c(0.22, 0.75)
      },
      legend.background = element_blank(),
      legend.spacing.y = unit(0, "mm")
    )
}

################################################################################
#' @export
#'
plot_multirun_hist <- function(D.all) {
  mean_data <- D.all |>
    group_by(Model) |>
    summarise(
      Means = mean(Est),
      SDs = sd(Est)
    )

  ggplot(D.all, aes(x = Est, fill = Model)) +
    geom_density(position = "identity", alpha = 0.3, adjust = 1) +
    # geom_vline(data = mean_data, aes(xintercept = Means, color = Model)) +
    labs(x = "Estimate", y = "Density", fill = "Model") +
    theme(
      text = element_text(size = 20),
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      legend.background = element_blank(),
      legend.spacing.y = unit(0, "mm"),
      legend.box.background = element_rect(colour = "black")
    )
}

################################################################################
#' @export
#'
plot_model_proportions <- function() {
  # NEED TO UPDATE

  # number of counts across whole landscape for each covariate type
  slow.counts <- animal_data |>
    filter(lscape_type == "Slow") |>
    pull(t_spent)
  med.counts <- animal_data |>
    filter(lscape_type == "Medium") |>
    pull(t_spent)
  fast.counts <- animal_data |>
    filter(lscape_type == "Fast") |>
    pull(t_spent)

  Prop_all <- D.all |>
    dplyr::filter(Covariate == "Covariate") |>
    unnest_wider(Prop_speeds, names_sep = "_") |>
    rename(
      Slow = Prop_speeds_1,
      Medium = Prop_speeds_2,
      Fast = Prop_speeds_3
    ) |>
    select(Model, Slow, Medium, Fast) |>
    pivot_longer(!Model, names_to = "Speed", values_to = "Proportions") |>
    group_by(Model, Speed) |>
    summarise(
      Means = mean(Proportions),
      SDs = sd(Proportions),
      .groups = "drop"
    ) |>
    add_row(
      Model = "ABM",
      Speed = c("Slow", "Medium", "Fast"),
      Means = c(slow.counts, med.counts, fast.counts)
    )


  ggplot(Prop_all, aes(x = Speed, y = Means, fill = Model)) +
    geom_bar(stat = "identity", color = "black", position = position_dodge()) +
    geom_errorbar(aes(ymin = Means - SDs, ymax = Means + SDs),
      width = .2,
      position = position_dodge(.9), size = .7
    ) +
    # scale_y_continuous(limits=c(0, max(Prop_all$Means) +.01), expand = c(0, 0)) +
    labs(
      x = "Landscape Type",
      y = "Relative Distributions"
    ) +
    scale_fill_manual(values = c("grey40", fig_colors[1:4])) +
    theme(
      text = element_text(size = 20),
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      legend.position = c(0.87, 0.7),
      legend.background = element_blank(),
      legend.spacing.y = unit(0, "mm"),
      legend.box.background = element_rect(colour = "black")
    )
}

################################################################################
#' @export
#'
plot_ABM_stay_proportions <- function(stay_time_raw_in) {
  # NEED TO UPDATE
  # number of counts across whole landscape for each covariate type
  slow.counts <- animal_data |>
    filter(lscape_type == "Slow") |>
    pull(t_spent)
  med.counts <- animal_data |>
    filter(lscape_type == "Medium") |>
    pull(t_spent)
  fast.counts <- animal_data |>
    filter(lscape_type == "Fast") |>
    pull(t_spent)


  Prop_all <- stay_time_raw_in |>
    dplyr::group_by(speed) |>
    dplyr::summarise(
      Model = "Stay Time Data",
      Means = mean(t_stay),
      SDs = sd(t_stay),
      .groups = "drop"
    ) |>
    dplyr::add_row(
      Model = "Total Habitat Use",
      speed = c("Slow", "Medium", "Fast"),
      Means = c(slow.counts, med.counts, fast.counts)
    ) |>
    dplyr::group_by(Model) |>
    dplyr::mutate(Means = Means / sum(Means))

  Prop_all$speed <- factor(Prop_all$speed, c("Slow", "Medium", "Fast"))
  Prop_all$Model <- factor(Prop_all$Model, c("Total Habitat Use", "Stay Time Data"))


  ggplot(Prop_all, aes(x = speed, y = Means, fill = Model)) +
    geom_bar(stat = "identity", color = "black", position = position_dodge()) +
    # geom_errorbar(aes(ymin=Means-SDs, ymax=Means+SDs), width=.2,
    #               position=position_dodge(.9), size = .7) +
    labs(
      x = "Habitat Type",
      y = "Proportion Occupied"
    ) +
    scale_fill_manual(values = c("grey", "grey40")) +
    scale_y_continuous(
      limits = c(0, NA),
      expand = expansion(mult = 0, add = c(0, 0.05))
    ) +
    theme(
      text = element_text(size = 20),
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      legend.position = c(0.85, 0.87),
      legend.background = element_blank(),
      legend.spacing.y = unit(0, "mm")
    )
}

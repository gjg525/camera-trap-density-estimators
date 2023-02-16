

################################
## Plot ABM simulations
################################
plot_ABM <- function() {
  # Sample area border
  b.df <- data.frame(c(bounds,dx))

  # camera cell locations
  cam.df <- tibble(cam = 1:ncam,
                   xmin = dx*(cam.samps - 1) %% q^0.5,
                   ymin = dx*ceiling(cam.samps / q^0.5) - dx,
                   xmax = xmin + dx,
                   ymax = ymin + dy,
                   a = (ymax - ymin) * (xmax - xmin)
  )
  # Triangle viewshed cameras
  dt.triangle <- data.frame(group = rep(1:ncam, each = 3),
                            polygon.x = tri_cam_samps$x,
                            polygon.y = tri_cam_samps$y)


  g <- ggplot(lscape_speeds) +
    # geom_tile(aes(X-.5*dx, Y-.5*dy, fill = Speed)) +
    scale_fill_manual(values= c("grey20", "grey50", "grey80")) +
    geom_path(data = animalxy.all,
              aes(x, y, col = as.factor(Animal_ID))) +
    geom_rect(data = b.df,
              aes(xmin = bounds[1], xmax = bounds[2], ymin = bounds[1], ymax = bounds[2]),
              color="black",
              fill= NA) +
    # geom_rect(data = cam.df,
    #           aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax),
    #           color="black",
    #           fill= NA) +
    geom_polygon(data = dt.triangle,
                 aes(x=polygon.x, y=polygon.y, group=group),
                 color = "black",
                 fill = NA) +
    labs(x = "X", y = "Y", fill="Landscape Type") +
    guides(colour = "none",
           alpha = "none") +
    theme(panel.background = element_blank())
  g + coord_fixed()
}

plot_space_use <- function() {
  u.abm.all <- raster(, nrows = q^0.5, ncols = q^0.5, xmn = 0,
                      xmx = q^0.5, ymn = 0, ymx = q^0.5, crs = NA)
  u.abm.all[] <- 0
  u.abm.all <- stack(replicate(t.steps,u.abm.all))

  u.temp = matrix(0, nrow = q^0.5, ncol = q^0.5)

  animalxy.discrete <- animalxy.all |>
    filter(t %in% 1:t.steps)
  for(xx in 1:nrow(animalxy.discrete)) {
    x.round <- ceiling(animalxy.discrete$x[xx]*(q^0.5/max(bounds)))
    y.round <- ceiling(animalxy.discrete$y[xx]*(q^0.5/max(bounds)))

    # Transpose for converting matrix to raster
    u.temp[x.round,y.round] <- u.temp[x.round,y.round]+1
  }

  u_plot <- levelplot(u.temp,
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


################################
## Plot data density plots
################################
# # Density plot for count data
plot_count_data <- function(fill = "speed") {
  ggplot(count_data, aes(x=count, fill = .data[[fill]])) +
    geom_density(position = "identity", alpha = 0.4, adjust = 3) +
    labs(x = "Count", y = "Frequency", fill="Landscape Type") +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=22),
          legend.text=element_text(size=20),
          legend.title = element_text(size=20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1))

  # ggplot(count_data, aes(x=count, fill = .data[[fill]])) +
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
plot_encounter_data <- function(fill = "speed") {
  ggplot(encounter_data, aes(x=encounter, fill = .data[[fill]])) +
    geom_density(position = "identity", alpha = 0.4, adjust = 2) +
    labs(x = "Number of Encounters", y = "Frequency", fill="Landscape Type") +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=22),
          legend.text=element_text(size=20),
          legend.title = element_text(size=20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1))

  # ggplot(encounter_data, aes(x=encounter, fill = .data[[fill]])) +
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
plot_staytime_data <- function(fill = "speed") {
  ggplot(stay_time_raw, aes(x=t_stay, fill = .data[[fill]])) +
    geom_density(position = "identity", alpha = 0.4, adjust = 5) +
    labs(x = "Staying Time", y = "Frequency", fill="Landscape Type") +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=22),
          legend.text=element_text(size=20),
          legend.title = element_text(size=20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1))
}

# # Density plot for TTE data
plot_TTE_data <- function(fill = "speed") {
  ggplot(TTE_data_raw, aes(x=TTE, fill = .data[[fill]])) +
    geom_density(position = "identity", alpha = 0.4, adjust = 2) +
    labs(x = "Time to Encounter", y = "Frequency", fill = "Landscape Type") +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=22),
          legend.text=element_text(size=20),
          legend.title = element_text(size=20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1))
}

################################
## Plot single-run results
################################

plot_onerun_results <- function() {
  ggplot(D.all, aes(x = Model)) +
    geom_point(aes(y = Est)) +
    geom_hline(yintercept=nind, linetype="dashed") +
    geom_errorbar(aes(ymin = Est - SD, ymax = Est + SD), width = 0.2)
}


################################
## Plot multi-run results
################################

plot_multirun_means <- function() {
  ggplot(D.all, aes(x = Model, y = Est, fill=Covariate)) +
    geom_boxplot(lwd = .1, fatten = .1) +
    geom_hline(yintercept=nind, linetype="dashed", size=1) +
    labs(x = "Model",
         y = "Abundance Estimates") +
    # y = paste(cam.props.label[cam.dist.set],"\n Abundance Estimates")) +
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
}

plot_multirun_sds <- function() {
  ggplot(D.all, aes(x = Model, y = SD)) +
    geom_boxplot(lwd = .1, fatten = .1, fill = "Gray") +
    labs(x = "Model",
         y = paste(cam.props.label[cam.dist.set],"\n SD Results")) +
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
}

plot_multirun_hist <- function() {
  mean_data <- D.all |>
    group_by(Model) |>
    summarise(Means = mean(Est),
              SDs = sd(Est))

  ggplot(D.all, aes(x = Est, fill = Model)) +
    geom_density(position = "identity", alpha = 0.3, adjust = 1) +
    # geom_vline(data = mean_data, aes(xintercept = Means, color = Model)) +
    labs(x = "Estimate", y = "Density", fill="Model") +
    theme(text = element_text(size = 20),
          legend.title=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.background = element_blank(),
          legend.spacing.y = unit(0, "mm"),
          legend.box.background = element_rect(colour = "black"))
}


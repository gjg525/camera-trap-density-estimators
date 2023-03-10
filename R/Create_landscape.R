# Create landscapes for ABMs

lscape_creator <- function(tag = "Random", speed_index = 1) {
  # Speed bounds for each landscape type
  speed_bounds <- data.frame(Speed = character(),
                             Min = double(),
                             Max = double()
  )
  speed_bounds <- dplyr::add_row(speed_bounds,
                                 Speed = "Slow",
                                 Min = 0.02,
                                 Max = 0.04)
  speed_bounds <- dplyr::add_row(speed_bounds,
                                 Speed = "Medium",
                                 Min = 0.09,
                                 Max = 0.11)
  speed_bounds <- dplyr::add_row(speed_bounds,
                                 Speed = "Fast",
                                 Min = 0.4,
                                 Max = 0.6)
  # speed_bounds <- dplyr::add_row(speed_bounds,
  #                                Speed = "Slow",
  #                                Min = 0.09,
  #                                Max = 0.11)
  # speed_bounds <- dplyr::add_row(speed_bounds,
  #                                Speed = "Medium",
  #                                Min = 0.3,
  #                                Max = 0.5)
  # speed_bounds <- dplyr::add_row(speed_bounds,
  #                                Speed = "Fast",
  #                                Min = 0.9,
  #                                Max = 1.1)
  # speed_bounds <- dplyr::add_row(speed_bounds,
  #                                Speed = "Slow",
  #                                Min = 0.2,
  #                                Max = 0.3)
  # speed_bounds <- dplyr::add_row(speed_bounds,
  #                                Speed = "Medium",
  #                                Min = 0.2,
  #                                Max = 0.3)
  # speed_bounds <- dplyr::add_row(speed_bounds,
  #                                Speed = "Fast",
  #                                Min = 0.2,
  #                                Max = 0.3)
  speed_bounds$Min <- speed_bounds$Min*dx
  speed_bounds$Max <- speed_bounds$Max*dx


  #############################################################################
  # Create dataframe for randomly-distributed speeds
  #############################################################################
  if (tag == "Random") {
    lscape_speeds <- data.frame(Index = 1:q,
                                Speed_index = sample(nrow(speed_bounds),
                                                     q,
                                                     replace = T)) |>
      dplyr::mutate(Speed = speed_bounds$Speed[Speed_index],
             Value = runif(q,
                           speed_bounds$Min[Speed_index],
                           speed_bounds$Max[Speed_index]),
             X = (Index - 1) %% q^0.5 + 1,
             Y = ceiling(Index / q^0.5),
             Direction = NA,
             Kappa = default_kappa,
             Road = "Off Trail")
  }

  #############################################################################
  # # Create dataframe for homogeneous landscape
  #############################################################################
  # All spaces have the same speed
  if (tag == "Homogeneous") {
    all_inds <- c(1:q)
  
    lscape_speeds <- data.frame(Index = 1:q,
                                Speed_index = NA)
  
    lscape_speeds$Speed_index[all_inds] <- speed_index
  
  
    lscape_speeds <- lscape_speeds |>
      dplyr::mutate(Speed = speed_bounds$Speed[Speed_index],
             Value = runif(q,
                           speed_bounds$Min[Speed_index],
                           speed_bounds$Max[Speed_index]),
             X = (Index - 1) %% q^0.5 + 1,
             Y = ceiling(Index / q^0.5),
             Direction = NA,
             Kappa = default_kappa,
             Road = "Off Trail")
  }
  return(lscape_speeds)
}

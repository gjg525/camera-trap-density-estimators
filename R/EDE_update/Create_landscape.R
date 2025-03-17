# Create landscapes for ABMs

################################################################################
#' @export
#'
lscape_creator <- function(study_design, lscape_design) {
  q <- study_design$q
  tag <- lscape_design$lscape_tag
  num_roads <- lscape_design$num_roads
  
  speed_bounds <- tibble::tibble(
    Speed = unlist(lscape_design$Speed_ID),
    Min = unlist(lscape_design$Speed_mins),
    Max = unlist(lscape_design$Speed_maxes)
  )

  lscape_defs <- tibble::tibble(
    Index = 1:q,
    X = (Index - 1) %% q^0.5 + 1,
    Y = ceiling(Index / q^0.5),
    Direction = NA,
    Kappa = lscape_design$default_kappa,
    Road = "Off Trail"
  )

  road_speeds <- tibble::tibble(
    Road = unlist(lscape_design$Trail_ID),
    Speed = unlist(lscape_design$Trail_speed)
  )|>
    dplyr::left_join(speed_bounds, by = "Speed")

  #############################################################################
  # Create dataframe for randomly-distributed speeds (default)
  #############################################################################
  if (tag == "Random") {
    lscape_defs <- lscape_defs |>
      dplyr::bind_cols(dplyr::sample_n(speed_bounds, q, replace = T)) |>
      dplyr::mutate(
        Value = runif(
          q,
          Min,
          Max
        )
      )

    # lscape_defs <- data.frame(Index = 1:q,
    #                             Speed_index = sample(nrow(speed_bounds),
    #                                                  q,
    #                                                  replace = T)) |>
    #   mutate(Speed = speed_bounds$Speed[Speed_index],
    #          Value = runif(q,
    #                        speed_bounds$Min[Speed_index],
    #                        speed_bounds$Max[Speed_index]),
    #          X = (Index - 1) %% q^0.5 + 1,
    #          Y = ceiling(Index / q^0.5),
    #          Direction = NA,
    #          Kappa = default_kappa,
    #          Road = "Off Trail")

    # #############################################################################
    # # # Create custom landscape
    # #############################################################################
  } else if (tag == "Random road") {
    trail_length <- 15
    x_lims_vert <- seq(min(lscape_defs$X) + 1, max(lscape_defs$X) - trail_length - 1)
    y_lims_vert <- seq(min(lscape_defs$Y) + 1, max(lscape_defs$Y) - trail_length - 1)
    x_lims_hor <- seq(min(lscape_defs$X) + 1, max(lscape_defs$X) - trail_length - 1)
    y_lims_hor <- seq(min(lscape_defs$Y) + 1, max(lscape_defs$Y) - trail_length - 1)
    
    for (road in 1:num_roads) {
      # Vertical trails
      sample_x <- sample(x_lims_vert, 1)
      sample_y <- sample(y_lims_vert, 1)
      y_dir <- sample(c(pi/2, 3 * pi / 2), 1)
      
      lscape_defs <- Cell_bias_xy(
        lscape_defs,
        sample_x,
        sample_y:(sample_y + trail_length - 1),
        y_dir,
        10,
        "On Trail"
      )
      
      # remove x coords to avoid collision
      x_lims_vert <- x_lims_vert[!(x_lims_vert %in% (sample_x - 1):(sample_x + 1))]
      
      # horizontal trails
      sample_x <- sample(x_lims_hor, 1)
      sample_y <- sample(y_lims_hor, 1)
      x_dir <- sample(c(0, pi), 1)
      
      lscape_defs <- Cell_bias_xy(
        lscape_defs,
        sample_x:(sample_x + trail_length - 1),
        sample_y,
        x_dir,
        10,
        "On Trail"
      )
      
      # remove y coord to avoid collision
      y_lims_hor <- y_lims_hor[!(y_lims_hor %in% (sample_y - 1):(sample_y + 1))]
      
      if(length(x_lims_vert) == 0 | length(y_lims_hor) == 0) {
        break
      }
    }
  } else if (tag == "grid") {
    road_kappa <- 3
    start_inds <- round(seq(0, 29, length.out = num_roads + 2))
    road_inds_all <- c()
    
    for (road in 1:(num_roads / 2)) {
      lscape_defs <- Cell_bias(
        lscape_defs,
        (start_inds[2 * road] * 30 + 2):(start_inds[2 * road] * 30 + 29),
        0,
        road_kappa,
        "On Trail"
      )
      lscape_defs <- Cell_bias(
        lscape_defs,
        (start_inds[2 * road + 1] * 30 + 2):(start_inds[2 * road + 1] * 30 + 29),
        pi,
        road_kappa,
        "On Trail"
      )
      lscape_defs <- Cell_bias(
        lscape_defs,
        seq(start_inds[2 * road + 1] + 30, start_inds[2 * road + 1] + 840, by = 30),
        pi / 2,
        road_kappa,
        "On Trail"
      )
      lscape_defs <- Cell_bias(
        lscape_defs,
        seq(start_inds[2 * road] + 30, start_inds[2 * road] + 840, by = 30),
        3 * pi / 2,
        road_kappa,
        "On Trail"
      )
      
      road_inds_all <- c(
        road_inds_all, 
        (start_inds[2 * road] * 30 + 2):(start_inds[2 * road] * 30 + 29),
        (start_inds[2 * road + 1] * 30 + 2):(start_inds[2 * road + 1] * 30 + 29),
        seq(start_inds[2 * road + 1] + 30, start_inds[2 * road + 1] + 870, by = 30),
        seq(start_inds[2 * road] + 30, start_inds[2 * road] + 870, by = 30)
      )
    }
    
    road_overlap <- road_inds_all[duplicated(road_inds_all)]
    lscape_defs <- Cell_bias(lscape_defs, road_overlap, NA, 10, "On Trail")
    
    # Border movement
    lscape_defs <- Cell_bias(lscape_defs, seq(30, 900, by = 30), 3 * pi / 2, road_kappa, "On Trail")
    lscape_defs <- Cell_bias(lscape_defs, 2:30, pi, road_kappa, "On Trail")
    lscape_defs <- Cell_bias(lscape_defs, seq(1, 871, by = 30), pi / 2, road_kappa, "On Trail")
    lscape_defs <- Cell_bias(lscape_defs, 871:899, 0, road_kappa, "On Trail")
  } else if (tag == "cross") {
    start_inds <- 15
    road_inds_all <- c()
    
      lscape_defs <- Cell_bias(
        lscape_defs,
        (start_inds * 30 + 2):(start_inds * 30 + 29),
        0,
        10,
        "On Trail"
      )
      # lscape_defs <- Cell_bias(
      #   lscape_defs,
      #   (start_inds[2 + 1] * 30 + 2):(start_inds[2 + 1] * 30 + 29),
      #   pi,
      #   10,
      #   "On Trail"
      # )
      # lscape_defs <- Cell_bias(
      #   lscape_defs,
      #   seq(start_inds[2 + 1] + 30, start_inds[2 + 1] + 840, by = 30),
      #   pi / 2,
      #   10,
      #   "On Trail"
      # )
      lscape_defs <- Cell_bias(
        lscape_defs,
        seq(start_inds + 30, start_inds + 840, by = 30),
        3 * pi / 2,
        10,
        "On Trail"
      )
      
      road_inds_all <- c(
        road_inds_all, 
        (start_inds * 30 + 2):(start_inds * 30 + 29),
        # (start_inds[2 + 1] * 30 + 2):(start_inds[2 + 1] * 30 + 29),
        # seq(start_inds[2 + 1] + 30, start_inds[2 + 1] + 870, by = 30),
        seq(start_inds + 30, start_inds + 870, by = 30)
      )
    
    road_overlap <- road_inds_all[duplicated(road_inds_all)]
    lscape_defs <- Cell_bias(lscape_defs, road_overlap, NA, 10, "On Trail")
    
    # Border movement
    lscape_defs <- Cell_bias(lscape_defs, seq(30, 900, by = 30), 3 * pi / 2, 5, "On Trail")
    lscape_defs <- Cell_bias(lscape_defs, 2:30, pi, 5, "On Trail")
    lscape_defs <- Cell_bias(lscape_defs, seq(1, 871, by = 30), pi / 2, 5, "On Trail")
    lscape_defs <- Cell_bias(lscape_defs, 871:899, 0, 5, "On Trail")
  } else if (tag == "circ") {
    road_kappa <- 5
    # Outer trail
    lscape_defs <- Cell_bias_xy(
      lscape_defs,
      c(3:28),
      3,
      0,
      road_kappa,
      "On Trail"
    )
    lscape_defs <- Cell_bias_xy(
      lscape_defs,
      28,
      3:28,
      pi/2,
      road_kappa,
      "On Trail"
    )
    lscape_defs <- Cell_bias_xy(
      lscape_defs,
      c(3:28),
      28,
      pi,
      road_kappa,
      "On Trail"
    )
    lscape_defs <- Cell_bias_xy(
      lscape_defs,
      3,
      4:28,
      3*pi/2,
      road_kappa,
      "On Trail"
    )
    # 
    # Second trail
    lscape_defs <- Cell_bias_xy(
      lscape_defs,
      11,
      11:20,
      pi/2,
      road_kappa,
      "On Trail"
    )
    lscape_defs <- Cell_bias_xy(
      lscape_defs,
      12:20,
      11,
      pi,
      road_kappa,
      "On Trail"
    )
    lscape_defs <- Cell_bias_xy(
      lscape_defs,
      20,
      12:20,
      3*pi/2,
      road_kappa,
      "On Trail"
    )
    lscape_defs <- Cell_bias_xy(
      lscape_defs,
      11:19,
      20,
      0,
      road_kappa,
      "On Trail"
    )
    # # Middle trail
    # lscape_defs <- Cell_bias_xy(
    #   lscape_defs,
    #   8,
    #   8:23,
    #   pi/2,
    #   10,
    #   "On Trail"
    # )
    # lscape_defs <- Cell_bias_xy(
    #   lscape_defs,
    #   9:23,
    #   8,
    #   pi,
    #   10,
    #   "On Trail"
    # )
    # lscape_defs <- Cell_bias_xy(
    #   lscape_defs,
    #   23,
    #   9:23,
    #   3*pi/2,
    #   10,
    #   "On Trail"
    # )
    # lscape_defs <- Cell_bias_xy(
    #   lscape_defs,
    #   8:22,
    #   23,
    #   0,
    #   10,
    #   "On Trail"
    # )
    # 
    # # Inner trail
    # lscape_defs <- Cell_bias_xy(
    #   lscape_defs,
    #   13:18,
    #   13,
    #   0,
    #   10,
    #   "On Trail"
    # )
    # lscape_defs <- Cell_bias_xy(
    #   lscape_defs,
    #   18,
    #   13:18,
    #   pi/2,
    #   10,
    #   "On Trail"
    # )
    # lscape_defs <- Cell_bias_xy(
    #   lscape_defs,
    #   13:18,
    #   18,
    #   pi,
    #   10,
    #   "On Trail"
    # )
    # lscape_defs <- Cell_bias_xy(
    #   lscape_defs,
    #   13,
    #   14:18,
    #   3*pi/2,
    #   10,
    #   "On Trail"
    # )
    # 
    # # bisecting roads
    # # Inner
    # lscape_defs <- Cell_bias_xy(
    #   lscape_defs,
    #   15,
    #   12:15,
    #   3*pi/2,
    #   10,
    #   "On Trail"
    # )
    # lscape_defs <- Cell_bias_xy(
    #   lscape_defs,
    #   16:19,
    #   15,
    #   0,
    #   10,
    #   "On Trail"
    # )
    # lscape_defs <- Cell_bias_xy(
    #   lscape_defs,
    #   15,
    #   16:19,
    #   pi/2,
    #   10,
    #   "On Trail"
    # )
    # lscape_defs <- Cell_bias_xy(
    #   lscape_defs,
    #   12:14,
    #   15,
    #   pi,
    #   10,
    #   "On Trail"
    # )
    
    
  } else if (tag == "squares") {
    # Squares
    lscape_defs <- Cell_bias(
      lscape_defs,
      c(840:869, 272:299),
      0,
      10,
      "On Trail"
    )
    lscape_defs <- Cell_bias(
      lscape_defs,
      c(32:59, 632:659),
      pi,
      10,
      "On Trail"
    )
    lscape_defs <- Cell_bias(
      lscape_defs,
      # c(seq(35, 845, by = 30), seq(270, 630, by = 30)),
      c(seq(32, 842, by = 30), seq(59, 630, by = 30)),
      pi / 2,
      10,
      "On Trail"
    )
    lscape_defs <- Cell_bias(
      lscape_defs,
      c(seq(65, 885, by = 30), seq(271, 631, by = 30)),
      3 * pi / 2,
      10,
      "On Trail"
    )  
    } else if (tag == "metapop") {
    # "Metapopulation"
    # Initialize roads
    lscape_defs <- Cell_bias(lscape_defs, 8531:8570, 0, 10, "On Trail")
    lscape_defs <- Cell_bias(lscape_defs, 1531:1570, pi, 10, "On Trail")
    lscape_defs <- Cell_bias(lscape_defs, seq(3115, 7015, by = 100), pi / 2, 10, "On Trail")
    lscape_defs <- Cell_bias(lscape_defs, seq(3185, 7085, by = 100), 3 * pi / 2, 10, "On Trail")

    # direct individuals toward roads
    lscape_defs <- Cell_bias(
      lscape_defs,
      c(3101:3114, 7070:7084),
      0,
      1,
      "Off Trail"
    )
    lscape_defs <- Cell_bias(
      lscape_defs,
      c(seq(7031, 8431, by = 100), seq(70, 1470, by = 100)),
      pi / 2,
      1,
      "Off Trail"
    )
    lscape_defs <- Cell_bias(
      lscape_defs,
      c(3116:3130, 7086:7100),
      pi,
      1,
      "Off Trail"
    )
    lscape_defs <- Cell_bias(
      lscape_defs,
      c(seq(8631, 9931, by = 100), seq(1670, 3070, by = 100)),
      3 * pi / 2,
      1,
      "Off Trail"
    )

    # Create diagonal biases to make trajectories look more natural
    lscape_defs <- Cell_bias(
      lscape_defs,
      c(3001:3014, 2901:2914, seq(7030, 8430, by = 100), seq(7029, 8429, by = 100)),
      pi / 4,
      1,
      "Off Trail"
    )
    lscape_defs <- Cell_bias(
      lscape_defs,
      c(3016:3030, 2916:2930, seq(71, 1471, by = 100), seq(72, 1472, by = 100)),
      3 * pi / 4,
      1,
      "Off Trail"
    )
    lscape_defs <- Cell_bias(
      lscape_defs,
      c(7186:7200, 7286:7300, seq(1671, 3071, by = 100), seq(1672, 3072, by = 100)),
      5 * pi / 4,
      1,
      "Off Trail"
    )
    lscape_defs <- Cell_bias(
      lscape_defs,
      c(seq(8630, 9930, by = 100), seq(8629, 9929, by = 100), 7170:7184, 7270:7284),
      7 * pi / 4,
      1,
      "Off Trail"
    )


    # Define "Dead Spaces"
    lscape_defs <- Cell_bias(
      lscape_defs,
      c(seq(7070, 8470, by = 100), seq(8670, 9970, by = 100)),
      0,
      10,
      NA
    )
    lscape_defs <- Cell_bias(
      lscape_defs,
      c(seq(31, 1431, by = 100), seq(1631, 3131, by = 100)),
      pi,
      10,
      NA
    )
    lscape_defs <- Cell_bias(
      lscape_defs,
      c(7001:7014, 7016:7031),
      pi / 2,
      10,
      NA
    )
    lscape_defs <- Cell_bias(
      lscape_defs,
      c(3170:3184, 3186:3200),
      3 * pi / 2,
      10,
      NA
    )

    for (road in 32:69) {
      # # up-right
      lscape_defs <- Cell_bias(
        lscape_defs,
        c((road * 100 + 1):(road * 100 + 14), (road * 100 + 31):(road * 100 + 69)),
        pi / 4,
        10,
        NA
      )
      # # up-left
      lscape_defs <- Cell_bias(
        lscape_defs,
        (road * 100 + 16):(road * 100 + 31),
        3 * pi / 4,
        10,
        NA
      )
      # # down-right
      lscape_defs <- Cell_bias(
        lscape_defs,
        (road * 100 + 70):(road * 100 + 84),
        7 * pi / 4,
        10,
        NA
      )
      # # down-left
      lscape_defs <- Cell_bias(
        lscape_defs,
        (road * 100 + 86):(road * 100 + 100),
        5 * pi / 4,
        10,
        NA
      )
    }

    for (road in 70:84) {
      # up-right
      lscape_defs <- Cell_bias(
        lscape_defs,
        (road * 100 + 32):(road * 100 + 69),
        pi / 4,
        10,
        NA
      )
    }
    for (road in 0:14) {
      # up-left
      lscape_defs <- Cell_bias(
        lscape_defs,
        (road * 100 + 32):(road * 100 + 69),
        3 * pi / 4,
        10,
        NA
      )
    }
    for (road in 86:99) {
      # down-right
      lscape_defs <- Cell_bias(
        lscape_defs,
        (road * 100 + 32):(road * 100 + 69),
        7 * pi / 4,
        10,
        NA
      )
    }
    for (road in 16:31) {
      # down-left
      lscape_defs <- Cell_bias(
        lscape_defs,
        (road * 100 + 32):(road * 100 + 69),
        5 * pi / 4,
        10,
        NA
      )
    }
  }

  if (tag != "Random") {
    lscape_defs <- dplyr::left_join(lscape_defs, road_speeds, by = "Road") |>
      dplyr::mutate(
        Min = tidyr::replace_na(Min, 0),
        Max = tidyr::replace_na(Max, 0),
        Value = runif(q, Min, Max)
      )
  }

  return(lscape_defs)
}

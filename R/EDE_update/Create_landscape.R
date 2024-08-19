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
  } else if (tag == "grid") {
    # Gridded road

    # Parms for intersecting roads landscape
    start_inds <- round(seq(0, 99, length.out = num_roads + 2))

    # # # There's probably a better way to do this
    road_inds_all <- c()
    for (road in 1:(num_roads / 2)) {
      lscape_defs <- Cell_bias(
        lscape_defs,
        (start_inds[2 * road] * 100 + 2):(start_inds[2 * road] * 100 + 99),
        0,
        10,
        "On Trail"
      )
      lscape_defs <- Cell_bias(
        lscape_defs,
        (start_inds[2 * road + 1] * 100 + 2):(start_inds[2 * road + 1] * 100 + 99),
        pi,
        10,
        "On Trail"
      )
      lscape_defs <- Cell_bias(
        lscape_defs,
        seq(start_inds[2 * road + 1] + 100, start_inds[2 * road + 1] + 9800, by = 100),
        pi / 2,
        10,
        "On Trail"
      )
      lscape_defs <- Cell_bias(
        lscape_defs,
        seq(start_inds[2 * road] + 100, start_inds[2 * road] + 9800, by = 100),
        3 * pi / 2,
        10,
        "On Trail"
      )
      # Come back to solve for road intersections
      road_inds_all <- c(
        road_inds_all, (start_inds[2 * road] * 100 + 2):(start_inds[2 * road] * 100 + 99),
        (start_inds[2 * road + 1] * 100 + 2):(start_inds[2 * road + 1] * 100 + 99),
        seq(start_inds[2 * road + 1] + 100, start_inds[2 * road + 1] + 9800, by = 100),
        seq(start_inds[2 * road] + 100, start_inds[2 * road] + 9800, by = 100)
      )
    }
    road_overlap <- road_inds_all[duplicated(road_inds_all)]
    lscape_defs <- Cell_bias(lscape_defs, road_overlap, 0, 1, "On Trail")

    # Define movement around border
    lscape_defs <- Cell_bias(
      lscape_defs,
      1:99,
      0,
      10,
      "On Trail"
    )
    lscape_defs <- Cell_bias(
      lscape_defs,
      seq(100, 9900, by = 100),
      pi / 2,
      10,
      "On Trail"
    )
    lscape_defs <- Cell_bias(
      lscape_defs,
      9902:10000,
      pi,
      10,
      "On Trail"
    )
    lscape_defs <- Cell_bias(
      lscape_defs,
      seq(101, 9901, by = 100),
      3 * pi / 2,
      10,
      "On Trail"
    )
  } else if (tag == "circ") {
    # Circular Movement
    # Horizontal roads
    lscape_defs <- Cell_bias(
      lscape_defs,
      c(9536:9565, seq(1064, 9964, by = 100)),
      0,
      10,
      "On Trail"
    )
    lscape_defs <- Cell_bias(
      lscape_defs,
      c(536:565, seq(35, 9035, by = 100)),
      pi,
      10,
      "On Trail"
    )

    # 'Reflective bounds' on road borders
    for (road in c(0:4, 51:94)) {
      lscape_defs <- Cell_bias(
        lscape_defs,
        (road * 100 + 36):(road * 100 + 65),
        pi / 2,
        10,
        "On Trail"
      )
    }
    for (road in c(6:50, 96:99)) {
      lscape_defs <- Cell_bias(
        lscape_defs,
        (road * 100 + 36):(road * 100 + 65),
        3 * pi / 2,
        10,
        "On Trail"
      )
    }
    for (road in 0:50) {
      lscape_defs <- Cell_bias(
        lscape_defs,
        (road * 100 + 2):(road * 100 + 34),
        5 * pi / 8,
        1,
        "On Trail"
      )
    }
    for (road in 51:98) {
      lscape_defs <- Cell_bias(
        lscape_defs,
        (road * 100 + 2):(road * 100 + 34),
        3 * pi / 8,
        1,
        "On Trail"
      )
    }
    for (road in 51:99) {
      lscape_defs <- Cell_bias(
        lscape_defs,
        (road * 100 + 66):(road * 100 + 99),
        13 * pi / 8,
        1,
        "On Trail"
      )
    }
    for (road in 1:50) {
      lscape_defs <- Cell_bias(
        lscape_defs,
        (road * 100 + 66):(road * 100 + 99),
        11 * pi / 8,
        1,
        "On Trail"
      )
    }
  } else if (tag == "squares") {
    # Squares
    lscape_defs <- Cell_bias(
      lscape_defs,
      c(9505:9594, 3030:3069),
      0,
      10,
      "On Trail"
    )
    lscape_defs <- Cell_bias(
      lscape_defs,
      c(506:595, 7031:7070),
      pi,
      10,
      "On Trail"
    )
    lscape_defs <- Cell_bias(
      lscape_defs,
      c(seq(505, 9405, by = 100), seq(3070, 6970, by = 100)),
      pi / 2,
      10,
      "On Trail"
    )
    lscape_defs <- Cell_bias(
      lscape_defs,
      c(seq(695, 9595, by = 100), seq(3130, 7030, by = 100)),
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

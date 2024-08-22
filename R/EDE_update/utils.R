
################################################################################
#' @export
#'
create_cam_samp_design <- function(study_design,
                                   lscape_defs,
                                   cam_design) {
  q <- study_design$q
  dx <- study_design$dx
  dy <- study_design$dy
  cam_length <- cam_design$cam_length

  cam_inds <- sample_speeds(
    cam.dist.prop = cam_design,
    lscape_speeds = lscape_defs
    )
  
  # cam_inds <- place_cams(study_design,
  #                         lscape_defs,
  #                         cam_design)

  cam_locs <- tibble::tibble(
    cam_ID = 1:cam_design$ncam,
    lscape_index = cam_inds,
    x = ((cam_inds - 1) %% q^0.5 + 1) * dx - dx,
    y = ceiling(cam_inds / q^0.5) * dy - dy
  ) |>
    dplyr::group_by(cam_ID) |>
    dplyr::summarise(
      lscape_index = lscape_index,
      x = list(runif(1, x + cam_length * 0.5 + 0.01 * dx, x + dx - cam_length * 0.5 - 0.01 * dx) +
        c(0, -0.5 * cam_length, 0.5 * cam_length)),
      y = list(runif(1, y + .01 * dy, y + dx - cam_length - 0.01 * dy) +
        c(0, cam_length, cam_length)),
      vertex = list(c(1, 2, 3)),
      cam_area = calc_tri_area(unlist(x), unlist(y)),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      Speed = lscape_defs$Speed[lscape_index],
      Road = lscape_defs$Road[lscape_index]
    )


  return(cam_locs)
}

#' ################################################################################
#' #' @export
#' #'
#' place_cams <- function(study_design,
#'                        lscape_defs,
#'                        cam_design) {
#'   # COME BACK TO CLEAN AND ALLOW FOR CUSTOM CAMERA PLACEMENT
#' 
#'   cam_design_name <- cam_design$Design[1]
#' 
#'   # Create camera sample designs
#'   if (cam_design_name == "Random") {
#'     # # Random camera sampling
#'     cam_samps <- sample(which(!is.na(lscape_defs$Road)), cam_design$ncam, replace = F)
#'   } else if (cam_design_name == "Speed") {
#'     # Biased sample design (lscape_defs)
#'     sample_speeds()
#'   } else if (cam_design_name == "Road") {
#'     cam_samps <- c(
#'       sample(which(lscape_defs$Road == "On Trail"),
#'         round(cam_design$ncam * cam_design$Props[1]),
#'         replace = F
#'       ),
#'       sample(which(lscape_defs$Road == "Off Trail"),
#'              cam_design$ncam - round(cam_design$ncam * cam_design$Props[1]),
#'         replace = F
#'       )
#'     )
#'   }
#' 
#'   return(cam_samps)
#' }

# Biased sample design (lscape_defs)
################################################################################
sample_speeds <- function(cam.dist.prop = NULL, lscape_speeds) {
  if (cam.dist.prop$Design == "Random") {
    cam.samps <- sample(lscape_speeds$Index, cam.dist.prop$ncam, replace = F)
  } else{
    ps <- cam.dist.prop$ncam * unlist(cam.dist.prop$Props)
    cam.samps <- c(sample(lscape_speeds |>
                            filter(Speed == "Slow") |>
                            pull(Index),
                          ps[1],
                          replace=F),
                   sample(lscape_speeds |>
                            filter(Speed == "Medium") |>
                            pull(Index),
                          ps[2],
                          replace=F),
                   sample(lscape_speeds |>
                            filter(Speed == "Fast") |>
                            pull(Index),
                          ps[3],
                          replace=F))
  }
}

################################################################################
#' @export
#'
place_animals <- function(study_design, lscape_defs) {
  # Randomly place across landscape
  q <- study_design$q
  num_groups <- study_design$num_groups
  group_sizes <- study_design$group_sizes
  group_spread <- study_design$group_spread
  tot_animals <- study_design$tot_animals
  road_mat <- matrix(lscape_defs$Road, nrow = q^0.5, q^0.5)

  # Note: Allow for custom placements
  animal.inds.0 <- tibble::tibble(
    IC_index = sample(which(!is.na(road_mat)), num_groups)
  ) |>
    dplyr::summarise(
      group_ID = 1:num_groups,
      group_size = unlist(group_sizes),
      X_ind = ((IC_index - 1) %% q^0.5 + 1),
      Y_ind = (ceiling(IC_index / q^0.5))
    )

  # Initial animal placement
  animalxy.0 <- animal.inds.0 |>
    dplyr::group_by(group_ID) |>
    dplyr::summarise(
      X = truncnorm::rtruncnorm(group_size,
        a = 0,
        b = q^0.5,
        X_ind - 0.5,
        group_spread
      ),
      Y = truncnorm::rtruncnorm(group_size,
        a = 0,
        b = q^0.5,
        Y_ind - 0.5,
        group_spread
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(Animal_ID = 1:tot_animals)

  # Make sure no animal is outside of bounds
  while (any(is.na(road_mat[cbind(ceiling(animalxy.0$X), ceiling(animalxy.0$Y))]))) {
    na_inds <- which(is.na(road_mat[cbind(ceiling(animalxy.0$X), ceiling(animalxy.0$Y))]))
    animalxy.0[na_inds, ] <- animalxy.0[na_inds, ] |>
      dplyr::rowwise() |>
      dplyr::mutate(
        X = truncnorm::rtruncnorm(1,
          a = 0,
          b = q^0.5,
          X - 0.5,
          0.5 * group_spread
        ),
        Y = truncnorm::rtruncnorm(1,
          a = 0,
          b = q^0.5,
          Y - 0.5,
          0.5 * group_spread
        )
      )
  }

  # Convert to x-y coords
  animalxy.0 <- animalxy.0 |>
    dplyr::mutate(
      X = X * study_design$dx,
      Y = Y * study_design$dy
    )

  return(animalxy.0)
}

################################################################################
#' @export
#'
create_activity_mat <- function(study_design,
                                animalxy.0) {

  if (study_design$activity_sync == "sync") {
    activity_mat <- rbinom(study_design$t_steps,
                           1,
                           unlist(study_design$activity_prob))

    animalxy.0 <- animalxy.0 |>
      # rowwise() |>
      dplyr::mutate(activity_mat = list(activity_mat))
  } else if(study_design$activity_sync == "non-sync") {
    animalxy.0 <- animalxy.0 |>
      dplyr::rowwise() |>
      dplyr::mutate(activity_mat = list(rbinom(study_design$t_steps,
                                               1,
                                               unlist(study_design$activity_prob))))
  } else if(study_design$activity_sync == "shift") {
    raw_activity_prob <- unlist(study_design$activity_prob)

    animalxy.0 <- animalxy.0 |>
      dplyr::rowwise() |>
      dplyr::mutate(
        shift_num = sample(1:(24/study_design$dt), 1),
        shift_activity_prob = list(c(raw_activity_prob[shift_num:length(raw_activity_prob)],
                                raw_activity_prob[1:(shift_num-1)])),
        activity_mat = list(rbinom(study_design$t_steps,
                                               1,
                                   unlist(shift_activity_prob))))
  }
}

################################################################################
#' @export
#'
get_covariate_info <- function(data_in, cov_name) {
  if (is.na(cov_name)) {
    covariate_labels <- NA
  } else {
    covariate_labels <- unique(data_in |> pull(cov_name))
    covariate_labels <- covariate_labels[!is.na(covariate_labels)]
  }

  return(covariate_labels)
}

################################################################################
#' @export
#'
create_covariate_mat <- function(lscape_defs, study_design, covariate_labels) {

  # Create covariate matrix with 0, 1 values
  if (any(is.na(covariate_labels))){
    Z <- matrix(0,
                nrow = study_design$q,
                ncol = 1)
    Z[!is.na(lscape_defs$Road)] <- 1

  } else {
    Z <- matrix(0,
                nrow = study_design$q,
                ncol = length(covariate_labels))
    for (zz in 1:length(covariate_labels)) {
      cov_vec <- lscape_defs |> 
        dplyr::filter(Speed == covariate_labels[zz]) %>% 
        dplyr::pull(Index)
      Z[cov_vec, zz] <- 1
    }
  }
  return(Z)
}

################################################################################
#' @export
#'
exp_na_covs <- function(Z, coeffs) {
  Z_coeffs <- Z %*% coeffs
  Z_coeffs[rowSums(Z) == 0] <- NA
  out <- exp(Z_coeffs)
}

###################################
# Data collection functions
###################################
################################################################################
#' @export
#'
calc_tri_area <- function(x, y) {
  # Calculate area of triangle using vertices
  tri_area <- 0.5 * abs(x[1] * (y[2] - y[3]) + x[2] * (y[3] - y[1]) + x[3] * (y[1] - y[2]))
}

################################################################################
#' @export
#'
calc_tri_signs <- function(p1, p2, p3) {
  tri_sign <- (p1[1] - p3[1]) * (p2[2] - p3[2]) - (p2[1] - p3[1]) * (p1[2] - p3[2])
}

# Determine if a point lies within a triangle using where the point lies relative to the triangle half-planes
################################################################################
#' @export
#'
calc_in_triangle <- function(p_viewshed, point) {
  A_1 <- calc_tri_signs(point, p_viewshed[1, ], p_viewshed[2, ])
  A_2 <- calc_tri_signs(point, p_viewshed[3, ], p_viewshed[1, ])
  A_3 <- calc_tri_signs(point, p_viewshed[2, ], p_viewshed[3, ])

  neg_vals <- (A_1 < 0) || (A_2 < 0) || (A_3 < 0)
  pos_vals <- (A_1 > 0) || (A_2 > 0) || (A_3 > 0)

  in_triangle <- ifelse(!(neg_vals && pos_vals),
    TRUE,
    FALSE
  )
}

# Calculate where an animal's trajectory intersects with a camera viewshed and the time spent within viewshed
################################################################################
#' @export
#'
calc_intersects <- function(p_viewshed, p_animal, speed, t) {
  X <- c()
  Y <- c()
  t_stay <- c()
  in_cam_all <- c()
  t_in_all <- c()
  t_out_all <- c()
  # Capture encounters as an individual leaves the viewshed
  encounter <- 1
  num_encounters <- 1

  # Get lengths of all edges
  r.length <- rbind(
    p_viewshed[2, ] - p_viewshed[1, ],
    p_viewshed[3, ] - p_viewshed[2, ],
    p_viewshed[1, ] - p_viewshed[3, ]
  )

  for (nn in 1:(nrow(p_animal) - 1)) {
    X.temp <- c()
    Y.temp <- c()
    t_tot <- 0 # time spent in camera
    t_in <- NA
    t_out <- NA

    animal.step <- p_animal[nn + 1, ] - p_animal[nn, ]

    # check for intersections with viewshed
    # Don't use intersections when start point is on the line (a,b = 0),
    # or if a full step lands on a line (a,b = 1)
    a <- c()
    b <- c()
    int.check <- c()
    for (rr in 1:3) {
      a[rr] <- det(cbind(p_animal[nn, ] - p_viewshed[rr, ], r.length[rr, ])) /
        det(cbind(r.length[rr, ], animal.step))
      b[rr] <- det(cbind((p_animal[nn, ] - p_viewshed[rr, ]), animal.step)) /
        det(cbind(r.length[rr, ], animal.step))
    }

    int.check <- a < 1 & a > 0 & b < 1 & b > 0 |
      (a == 0 & b < 1 & b > 0) |
      (a == 1 & b < 1 & b > 0) |
      (b == 0 & a < 1 & a > 0) |
      (b == 1 & a < 1 & a > 0)

    # Check if animal is already inside viewshed
    in_cam <- calc_in_triangle(p_viewshed, p_animal[nn, ])

    if (in_cam) {
      # Individual steps outside of camera
      if (any(int.check == T)) {
        int.ind <- which(int.check == T)

        # X,Y coordinates at intersection
        X.temp <- p_viewshed[int.ind, 1] + b[int.ind] * r.length[int.ind, 1]
        Y.temp <- p_viewshed[int.ind, 2] + b[int.ind] * r.length[int.ind, 2]

        # Find distance from start point to intersection
        dist_1 <- sqrt((X.temp - p_animal[nn, 1])^2 +
          (Y.temp - p_animal[nn, 2])^2)

        # Calculate time spent in camera
        t_tot <- ifelse(speed[nn + 1] > 0,
                        dist_1 / speed[nn + 1],
                        t[nn + 1] - t[nn])

        # Calculate enter/exit times
        t_in <- NA
        t_out <- t[nn] + t_tot

        num_encounters <- num_encounters + 1
      } else {
        # Individual remains in camera

        # Find distance from start point to end point
        dist_1 <- sqrt((p_animal[nn + 1, 1] - p_animal[nn, 1])^2 +
          (p_animal[nn + 1, 2] - p_animal[nn, 2])^2)

        # Calculate time spent in camera
        t_tot <- ifelse(speed[nn + 1] > 0,
                        dist_1 / speed[nn + 1],
                        t[nn + 1] - t[nn])

        # Encounter if individual remains in viewshed on final step
        if (t[nn] == t[length(t) - 1]) {
          num_encounters <- num_encounters + 1
        }
      }
    } else {
      # Individual starts outside of camera
      if (sum(int.check) == 1) {
        # Individual steps into and stays in camera
        int.ind <- which(int.check == T)

        # X,Y coordinates at intersection
        X.temp <- p_viewshed[int.ind, 1] + b[int.ind] * r.length[int.ind, 1]
        Y.temp <- p_viewshed[int.ind, 2] + b[int.ind] * r.length[int.ind, 2]

        # Find distance from intersection to end point
        dist_1 <- sqrt((p_animal[nn + 1, 1] - X.temp)^2 +
          (p_animal[nn + 1, 2] - Y.temp)^2)

        # Calculate time spent in camera
        t_tot <- ifelse(speed[nn + 1] > 0,
                        dist_1 / speed[nn + 1],
                        t[nn + 1] - t[nn])

        # Find distance from start point to intersection
        dist_2 <- sqrt((p_animal[nn, 1] - X.temp)^2 +
          (p_animal[nn, 2] - Y.temp)^2)

        # Calculate enter/exit times
        t_in <- t[nn] + ifelse(speed[nn + 1] > 0,
                               dist_2 / speed[nn + 1],
                               t[nn + 1] - t[nn])
        # t_out <- NA

        # Encounter if individual remains in viewshed on final step
        if (t[nn] == t[length(t) - 1]) {
          num_encounters <- num_encounters + 1
        }
      } else if (sum(int.check) == 2) {
        # Individual passes through camera
        int.ind <- which(int.check == T)

        # X,Y coordinates at intersection
        X.temp <- p_viewshed[int.ind, 1] + b[int.ind] * r.length[int.ind, 1]
        Y.temp <- p_viewshed[int.ind, 2] + b[int.ind] * r.length[int.ind, 2]

        # Find distance between points of intersection
        dist_1 <- sqrt((X.temp[1] - X.temp[2])^2 +
          (Y.temp[1] - Y.temp[2])^2)

        # Calculate time spent in camera
        t_tot <- ifelse(speed[nn + 1] > 0,
                        dist_1 / speed[nn + 1],
                        t[nn + 1] - t[nn])

        # Find distance from start point to first intersection
        dist_2 <- sqrt((p_animal[nn, 1] - X.temp[1])^2 +
          (p_animal[nn, 2] - Y.temp[1])^2)

        # Calculate enter/exit times
        t_in <- t[nn] + ifelse(speed[nn + 1] > 0,
                               dist_2 / speed[nn + 1],
                               t[nn + 1] - t[nn])
        t_out <- t_in + t_tot

        num_encounters <- num_encounters + 1
      }
    }

    # Accumulate all values
    X <- c(X, X.temp)
    Y <- c(Y, Y.temp)
    t_stay <- c(t_stay, t_tot)
    in_cam_all <- c(in_cam_all, in_cam)
    encounter <- c(encounter, num_encounters)
    t_in_all <- c(t_in_all, t_in)
    t_out_all <- c(t_out_all, t_out)
  }

  # encounter <- c(encounter, encounter[length(encounter)])
  return(list(cbind(X, Y), t_stay, in_cam_all, encounter, t_in_all, t_out_all))
}

###################################
# Landscape creator helper
###################################
# Define individual cells
################################################################################
#' @export
#'
Cell_bias <- function(lscape_data, Index, Dir, Kappa, Road_ID) {
  lscape_data$Direction[Index] <- Dir
  lscape_data$Kappa[Index] <- Kappa
  lscape_data$Road[Index] <- Road_ID

  return(lscape_data)
}

########################################
# Collect space-use for each animal
########################################
################################################################################
get_animal_data <- function(animalxy.all, lscape_defs) {
  animal_data <- animalxy.all |>
    group_by(Animal_ID) |>
    mutate(t_diff = c(0,t[2:length(t)] - t[1:(length(t)-1)])) |>
    group_by(Animal_ID, lscape_type) |>
    summarise(t_spent = sum(t_diff),
              .groups = 'drop') |>
    dplyr::group_by(lscape_type) |>
    dplyr::summarise(nn = n(),
      prop_tot_time = sum(t_spent)) |>
  dplyr::left_join(lscape_defs |>
                     dplyr::rename("lscape_type" = "Speed") |>
                     dplyr::group_by(lscape_type) |>
                     dplyr::summarise(n_units = n()),
                   by = dplyr::join_by(lscape_type))

  return(animal_data)
}

################################################################################
calc_group_size <- function(cam_captures) {
  # Mean group size for adjusted STE
  mean_group_size <- cam_captures |>
    dplyr::filter(in_cam == T) |>
    dplyr::group_by(lscape_index, t) |>
    dplyr::summarise(
      group_size = n(),
      .groups = "drop"
    ) |>
    dplyr::summarise(
      mean_group_size = mean(group_size)
    ) |>
    pull(mean_group_size)

  return(mean_group_size)
}

################################################################################
summarise_data <- function(count_data,
                           encounter_data,
                           stay_time_raw
                           ) {
  count_data_summary <- count_data |>
    group_by(Road) |>
    summarise(
      Data_type = "Count",
      total = sum(count) / n()
    )
  count_data_summary$total <- count_data_summary$total / sum(count_data_summary$total)

  encounter_data_summary <- encounter_data |>
    group_by(Road) |>
    summarise(
      Data_type = "Encounter",
      total = sum(encounter) / n()
    )
  encounter_data_summary$total <- encounter_data_summary$total / sum(encounter_data_summary$total)

  staytime_data_summary <- stay_time_raw |>
    group_by(Road) |>
    summarise(
      Data_type = "Stay Time",
      total = sum(t_stay) / n()
    )
  staytime_data_summary$total <- staytime_data_summary$total / sum(staytime_data_summary$total)

  all_summaries <- rbind(count_data_summary, encounter_data_summary, staytime_data_summary)
}

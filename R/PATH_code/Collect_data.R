########################################
# Collect data for triangle viewshed cameras
########################################

################################################################################
    get_cam_captures <- function(animalxy) {
      t_steps <- study_design$t_steps
      dt <- study_design$dt

      # collect data within each grid cell to reduce subsequent calculations
      cell_captures <- animalxy |>
        dplyr::ungroup() |>
        dplyr::filter(lscape_index %in% cam_locs$lscape_index) |>
        dplyr::mutate(pass_i = cumsum(c(1, (abs(ii[-n()] - ii[-1]) > 1) |
          (abs(lscape_index[-n()] - lscape_index[-1]) > 0)))) |>
        dplyr::group_by(pass_i) |>
        dplyr::mutate(next_i = ifelse(max(t) == t_steps * dt,
          max(ii),
          max(ii) + 1
        )) |>
        dplyr::ungroup()

      # # Use animal's next step to capture all trajectories within cell
      cell_captures <- cell_captures |>
        dplyr::add_row(
          animalxy |>
            dplyr::ungroup() |>
            dplyr::filter(ii %in% unique(cell_captures$next_i)) |>
            dplyr::mutate(pass_i = unique(cell_captures$pass_i))
          ) |>
        dplyr::arrange(ii) |>
        dplyr::select(-next_i)
      
      # Remove repeat rows (occurs for max time step)
      cell_captures <- dplyr::distinct(cell_captures)

      cell_check <- cell_captures |>
        dplyr::group_by(pass_i) |>
        dplyr::mutate(
          num_animals = length(unique(Animal_ID)),
          num_lscapes = length(unique(lscape_index))
        )

      # Check that only 1 animal is included in each pass
      if (any(cell_check$num_animals > 1)) {
        stop("More than 1 animal in a pass")
      }
      # Check that only 2 cells are included in each pass
      if (any(cell_check$num_lscapes > 2)) {
        stop("More than 2 cells in a pass")
      }

      if (nrow(cell_captures) > 0) {
        # Collect all camera viewshed captures
        # in_cam determines whether an individual starts within the camera viewshed
        cam_captures <- cell_captures |>
          dplyr::group_by(pass_i) |>
          # dplyr::mutate(pass_i = pass_i[t == min(t)]) |>
          dplyr::summarise(
            xy_index = list(cam_locs[cam_locs$lscape_index %in% lscape_index[1], 3:4]),
            cam_intersects = list(
              calc_intersects(
                matrix(unlist(xy_index),
                  nrow = length(unlist(xy_index)) / 2,
                  ncol = 2
                ),
                cbind(X, Y),
                trav_speeds,
                t
              )
            ),
            t_stay = c(unlist(lapply(cam_intersects, "[[", 2)), NA),
            in_cam = c(unlist(lapply(cam_intersects, "[[", 3)), NA),
            encounter = unlist(lapply(cam_intersects, "[[", 4)),
            t_in = c(unlist(lapply(cam_intersects, "[[", 5)), NA),
            t_out = c(unlist(lapply(cam_intersects, "[[", 6)), NA),
            lscape_index = lscape_index[1],
            Animal_ID = Animal_ID,
            t = t,
            ii = ii,
            X = X,
            Y = Y,
            trav_speeds = trav_speeds,
            Speed = lscape_type[1],
            Road = Road,
            .groups = "drop"
          ) |>
          dplyr::filter(t_stay > 0) |>
          dplyr::group_by(pass_i) |>
          dplyr::mutate(pass_i = pass_i + 0.001 * cumsum(c(1, abs(encounter[-n()] - encounter[-1]) > 0))) |>
          dplyr::select(-cam_intersects) |>
          ungroup()
      } else {
        cam_captures <- cell_captures
      }
      return(cam_captures)
    }

################################################################################
    get_count_data <- function(cam_locs, cam_captures, animalxy) {
      # count_data <- animalxy |>
      #   dplyr::select(Animal_ID, group_ID, X, Y, t, lscape_index,lscape_type) %>%
      #   dplyr::filter(lscape_index %in% cam_locs$lscape_index &
      #                   t %in% 1:study_design$t_steps) |>
      #   dplyr::rowwise() %>%
      #   dplyr::mutate(
      #     xy_index = list(cam_locs[cam_locs$lscape_index %in% lscape_index, 3:4]),
      #     in_cam = calc_in_triangle(
      #       matrix(unlist(xy_index),
      #              nrow = length(unlist(xy_index)) / 2,
      #              ncol = 2
      #       ),
      #       cbind(X, Y)
      #     )
      #   ) %>%
      #   dplyr::filter(in_cam == T) %>%
      #   dplyr::group_by(lscape_index, t) |>
      #   dplyr::summarise(
      #     group_size = n(),
      #     .groups = "drop"
      #   ) |>
      #   dplyr::group_by(lscape_index) |>
      #   dplyr::summarise(
      #     count = sum(group_size),
      #     mean_group = mean(group_size),
      #     .groups = "drop"
      #   ) |>
      #   dplyr::full_join(cam_locs |>
      #                      select(cam_ID, lscape_index, Speed),
      #                    by = c("lscape_index")) |>
      #   dplyr::arrange(cam_ID) |>
      #   dplyr::mutate(count = replace(count, is.na(count), 0))
      
      count_data <- cam_captures |>
        dplyr::filter(in_cam == T) |>
        dplyr::group_by(lscape_index, t) |>
        dplyr::summarise(
          group_size = n(),
          .groups = "drop"
        ) |>
        dplyr::group_by(lscape_index) |>
        dplyr::summarise(
          count = sum(group_size),
          mean_group = mean(group_size),
          .groups = "drop"
        ) |>
        dplyr::full_join(cam_locs |>
                           select(cam_ID, lscape_index, Speed, Road),
                         by = c("lscape_index")) |>
        dplyr::arrange(cam_ID) |>
        dplyr::mutate(count = replace(count, is.na(count), 0))

      return(count_data)
    }

################################################################################
    get_encounter_data <- function(cam_locs, cam_captures) {
      encounter_data <- cam_captures |>
        dplyr::group_by(lscape_index) |>
        dplyr::summarise(
          encounter = length(unique(pass_i)),
          .groups = "drop"
        ) |>
        dplyr::full_join(cam_locs |>
                           select(cam_ID, lscape_index, Speed, Road),
                         by = c("lscape_index")) |>
        dplyr::arrange(cam_ID) |>
        dplyr::mutate(encounter = replace(encounter, is.na(encounter), 0))

      return(encounter_data)
    }


################################################################################
    get_stay_time_data <- function(cam_locs, cam_captures) {
      
      stay_time_raw <- cam_captures |>
        dplyr::group_by(pass_i) |>
        dplyr::summarise(
          t_stay = sum(t_stay),
          lscape_index = lscape_index[1],
          Speed = Speed[1],
          Road = Road[1],
          .groups = "drop"
        )

      stay_time_ref <- tibble::tibble(
        Speed = c("Slow", "Medium", "Fast"),
        mean_stay = 1
      ) %>%
        dplyr::filter(
          !(Speed %in% unique(stay_time_raw$Speed))
        )

      stay_time_normalize <- stay_time_raw %>%
        dplyr::group_by(Speed) %>%
        dplyr::summarise(
          mean_stay = mean(t_stay),
          .groups = 'drop'
        ) %>%
        dplyr::bind_rows(stay_time_ref) %>%
        dplyr::summarise(
          sum_stay = sum(mean_stay)
        ) %>%
        dplyr::pull(sum_stay)
      
      # Format staying time data
      stay_time_data <- stay_time_raw %>% 
        dplyr::mutate(t_stay = t_stay / stay_time_normalize) %>%
        # dplyr::add_row(lscape_index = cam_locs$lscape_index[cam_locs$lscape_index %notin% cam_captures$lscape_index]) |>
        dplyr::full_join(cam_locs |>
                           select(cam_ID, lscape_index),
                         by = "lscape_index") |>
        dplyr::arrange(cam_ID) |>
        dplyr::select(lscape_index, t_stay) |>
        dplyr::group_by(lscape_index) |>
        dplyr::mutate(encounter = 1:n()) |>
        dplyr::ungroup() |>
        pivot_wider(names_from = encounter, values_from = t_stay) |>
        dplyr::select(-lscape_index)

      return(list(stay_time_raw, stay_time_data))
}

################################################################################
    get_TTE_data <- function(study_design, cam_locs, cam_captures) {
      num_occ <- study_design$num_occ
      t_steps<- study_design$t_steps
      TTE_censor <- study_design$TTE_censor

      TTE_data <- data.frame(
        lscape_index = double(),
        TTE = double(),
        occasion = double()
      )

      # Loop through all occasions
      for (jj in 1:num_occ) {
        occ_times <- t_steps * c((jj - 1) / num_occ, jj / num_occ)
        TTE_temp <- cam_captures |>
          dplyr::filter(
            lscape_index %in% cam_locs$lscape_index,
            t >= occ_times[1],
            t < occ_times[2]
          )
        if (nrow(TTE_temp) > 0) {
          TTE_data <- rbind(TTE_data, TTE_temp |>
            dplyr::group_by(lscape_index) |>
            dplyr::summarise(
              TTE = min(t) - (jj - 1) * TTE_censor,
              occasion = jj
            ))
        }
        TTE_data <- dplyr::add_row(TTE_data,
          lscape_index = cam_locs$lscape_index[cam_locs$lscape_index %notin% TTE_temp$lscape_index],
          TTE = NA,
          occasion = jj
        )
      }
      TTE_data <- TTE_data |>
        dplyr::mutate(
          Speed = lscape_defs$Speed[lscape_index],
          Road = lscape_defs$Road[lscape_index]
        )

      TTE_data_raw <- TTE_data

      # Format TTE data
      TTE_data <- TTE_data |>
        dplyr::add_row(lscape_index = cam_locs$lscape_index[cam_locs$lscape_index %notin% TTE_data$lscape_index]) |>
        dplyr::left_join(cam_locs,
          by = "lscape_index"
        ) |>
        dplyr::select(lscape_index, TTE, occasion) |>
        dplyr::distinct() |>
        pivot_wider(names_from = occasion, values_from = TTE)

      TTE_data <- as.matrix(dplyr::left_join(data.frame(lscape_index = cam_locs$lscape_index),
        TTE_data,
        by = "lscape_index"
      ) |>
        dplyr::select(-lscape_index))
      TTE_data[TTE_data == 0] <- NA

      return(list(TTE_data_raw, TTE_data))
    }

################################################################################
    get_STE_data <- function(study_design, cam_locs, cam_captures) {
      cam_A <- cam_locs$cam_area[1]
      ncam <- dim(cam_locs)[1]
      t_steps <- study_design$t_steps
      dt <- study_design$dt

      STE_data <- cam_captures |>
        dplyr::filter(in_cam == T) |>
        dplyr::select(Animal_ID, t, lscape_index) |>
        dplyr::group_by(t) |>
        dplyr::summarise(
          cam_inds = list(sample(c(unique(lscape_index),
                                   rep(0, ncam - length(unique(lscape_index)))))),
          STE = min(which(unlist(cam_inds) > 0)) * cam_A,
          .groups = "drop"
        ) |>
        dplyr::select(t, STE) |>
        dplyr::full_join(tibble::tibble(
          t = seq(dt, t_steps * dt, by = dt)),
          by = "t") |>
      dplyr::arrange(t)

      return(STE_data)
}


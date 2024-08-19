########################################
# Collect data for triangle viewshed cameras
########################################
Collect_tele_data = function(animalxy.all, study_design, cam_locs) {
  # First, collect data within each grid cell to reduce subsequent calculations
  # Creates IDs (pass_i) for each pass through a cell and landscape indices for the cell that the individual passes into (next_i)
  cell_captures_tele <- animalxy.all |>
    ungroup() |> 
    dplyr::mutate(pass_i = cumsum(c(1, (abs(ii[-length(ii)] - ii[-1]) > 1) +
                                      (abs(Animal_ID[-length(Animal_ID)] - Animal_ID[-1]) > 0) +                                     (abs(lscape_index[-length(Animal_ID)] - lscape_index[-1]) > 0)))) |> 
    dplyr::group_by(pass_i) |>
    dplyr::mutate(next_i = ifelse(max(t) == study_design$t_steps,
                                  max(ii),
                                  max(ii) + 1)) |>
    dplyr::ungroup()
  
  # # Use animal's next step to capture all trajectories within cell
  cell_captures_tele <- cell_captures_tele |>
    dplyr::add_row(animalxy.all[unique(cell_captures_tele$next_i), ] |>
                     dplyr::ungroup() |>
                     dplyr::mutate(pass_i = unique(cell_captures_tele$pass_i))) |>
    dplyr::arrange(ii) |>
    dplyr::select(-next_i)
  
  # Remove repeat rows
  cell_captures_tele <- dplyr::distinct(cell_captures_tele)
  
  cell_check <- cell_captures_tele |>
    dplyr::group_by(pass_i) |>
    dplyr::mutate(num_animals = length(unique(Animal_ID)),
                  num_lscapes = length(unique(lscape_index)))
  # Check that only 1 animal is included in each pass
  if (any(cell_check$num_animals > 1)) {
    stop("More than 1 animal in a pass")
  }
  # Check that only 2 cells are included in each pass
  if (any(cell_check$num_lscapes > 2)) {
    stop("More than 2 cells in a pass")
  }
  
  ########################################
  # Collect stay time data for each encounter
  ########################################
  stay_time_raw_tele <- cell_captures_tele |>
    dplyr::group_by(pass_i) |>
    dplyr::summarise(t_stay = max(t) - min(t),
                     lscape_index = lscape_index[1],
                     speed = lscape_type[1],
                     .groups = 'drop')
  
  
  # # Format staying time data
  # stay_time_tele <- cell_captures_tele |>
  #   dplyr::group_by(pass_i) |>
  #   dplyr::summarise(t_stay = max(t) - min(t),
  #                    lscape_index = lscape_index[1]) |>
  #   dplyr::select(lscape_index, t_stay) |>
  #   dplyr::group_by(lscape_index) |>
  #   dplyr::mutate(encounter = 1:n()) |>
  #   dplyr::ungroup() |>
  #   pivot_wider(names_from = encounter, values_from = t_stay)
  # 
  # stay_time_tele <- as.matrix(dplyr::left_join(
  #   data.frame(lscape_index = cam_locs$lscape_index),
  #   stay_time_tele,
  #   by = "lscape_index") |>
  #                               dplyr::select(-lscape_index)
  # )
  
  return(stay_time_raw_tele)
}   
    
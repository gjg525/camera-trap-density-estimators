########################################
# Collect data telemetry
########################################
Collect_tele_data = function(animalxy.all, study_design) {
  # Collect animal location in each time step
  cell_captures_tele <- animalxy.all %>% 
    dplyr::filter(t %in% seq(1, study_design$t_steps, study_design$dt)) %>% 
    dplyr::rename(Speed = lscape_type) %>% 
    dplyr::group_by(Speed) %>% 
    dplyr::count() %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(
      stay_prop = n / sum(n),
      # Standard deviation for log-transformed counts
      stay_sd = sqrt(stay_prop * (1 - stay_prop) / sum(n)) / stay_prop
    ) %>% 
    dplyr::arrange(desc(Speed))
  
  # # First, collect data within each grid cell to reduce subsequent calculations
  # # Creates IDs (pass_i) for each pass through a cell and landscape indices for the cell that the individual passes into (next_i)
  # cell_captures_tele <- animalxy.all |>
  #   ungroup() |> 
  #   dplyr::mutate(pass_i = cumsum(
  #     c(1, 
  #       (abs(ii[-length(ii)] - ii[-1]) > 1) +
  #         (abs(Animal_ID[-length(Animal_ID)] - Animal_ID[-1]) > 0) +                              (abs(lscape_index[-length(Animal_ID)] - lscape_index[-1]) > 0)))) |> 
  #   dplyr::group_by(pass_i) |>
  #   dplyr::mutate(next_i = ifelse(max(t) == study_design$t_steps * study_design$dt,
  #                                 max(ii),
  #                                 max(ii) + 1)) |>
  #   dplyr::ungroup()
  # 
  # # # Use animal's next step to capture all trajectories within cell
  # cell_captures_tele <- cell_captures_tele |>
  #   dplyr::add_row(animalxy.all[unique(cell_captures_tele$next_i), ] |>
  #                    dplyr::ungroup() |>
  #                    dplyr::mutate(pass_i = unique(cell_captures_tele$pass_i))) |>
  #   dplyr::arrange(ii) |>
  #   dplyr::select(-next_i)
  # 
  # # Remove repeat rows
  # cell_captures_tele <- dplyr::distinct(cell_captures_tele)
  # 
  # cell_check <- cell_captures_tele |>
  #   dplyr::group_by(pass_i) |>
  #   dplyr::mutate(num_animals = length(unique(Animal_ID)),
  #                 num_lscapes = length(unique(lscape_index)))
  # # Check that only 1 animal is included in each pass
  # if (any(cell_check$num_animals > 1)) {
  #   stop("More than 1 animal in a pass")
  # }
  # # Check that only 2 cells are included in each pass
  # if (any(cell_check$num_lscapes > 2)) {
  #   stop("More than 2 cells in a pass")
  # }
  # 
  # stay_time_raw_tele <- cell_captures_tele |>
  #   dplyr::group_by(pass_i) |>
  #   dplyr::summarise(t_stay = max(t) - min(t),
  #                    lscape_index = lscape_index[1],
  #                    speed = lscape_type[1],
  #                    .groups = 'drop')
  
  
  return(cell_captures_tele)
}   
    
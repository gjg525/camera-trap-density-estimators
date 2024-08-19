########################################
# Collect data for triangle viewshed cameras
########################################
# First, collect data within each grid cell to reduce subsequent calculations
# Creates IDs (pass_i) for each pass through a cell and landscape indices for the cell that the individual passes into (next_i)
cell_captures_tele <- animalxy.all |>
  ungroup() |> 
  dplyr::mutate(pass_i = cumsum(c(1, (abs(ii[-length(ii)] - ii[-1]) > 1) +
                                    (abs(ID[-length(ID)] - ID[-1]) > 0) +                                     (abs(XY_inds[-length(ID)] - XY_inds[-1]) > 0)))) |> 
  dplyr::group_by(pass_i) |>
  dplyr::mutate(next_i = ifelse(max(t) == t.steps,
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
  dplyr::mutate(num_animals = length(unique(ID)),
                num_lscapes = length(unique(XY_inds)))
# Check that only 1 animal is included in each pass
if (any(cell_check$num_animals > 1)) {
  stop("More than 1 animal in a pass")
}
# Check that only 2 cells are included in each pass
if (any(cell_check$num_lscapes > 2)) {
  stop("More than 2 cells in a pass")
}

# # Collect all camera viewshed captures
# # in_cam determines whether an individual starts within the camera viewshed
# cam_captures_tele <- cell_captures_tele |>
#   dplyr::group_by(pass_i) |>
#   dplyr::summarise(cam_samps = list(tri_cam_samps[tri_cam_samps$XY_inds %in% XY_inds[1], 3:4]),
#                    cam_intersects = list(calc_intersects(matrix(unlist(cam_samps), 
#                                                                 nrow = length(unlist(cam_samps))/2, ncol = 2),
#                                                          t(rbind(x, y)),
#                                                          speeds,
#                                                          t)),
#                    # intersect = c(lapply(cam_intersects, '[[', 1), NA),
#                    # t_stay = sum(unique(unlist(lapply(cam_intersects, '[[', 2)))),
#                    t_stay = c(unlist(lapply(cam_intersects, '[[', 2)), NA),
#                    in_cam = c(unlist(lapply(cam_intersects, '[[', 3)), NA),
#                    encounter = unlist(lapply(cam_intersects, '[[', 4)),
#                    t_in = c(unlist(lapply(cam_intersects, '[[', 5)), NA),
#                    t_out = c(unlist(lapply(cam_intersects, '[[', 6)), NA),
#                    XY_inds = XY_inds[1],
#                    ID = ID,
#                    t = t,
#                    ii = ii,
#                    x = x,
#                    y = y,
#                    lscape_type = lscape_type[1],
#                    road = road,
#                    .groups = 'drop'
#   )  |>
#   dplyr::filter(t_stay > 0) |>
#   dplyr::group_by(pass_i) |>
#   dplyr::mutate(pass_i = pass_i + 0.001 * cumsum(c(1, abs(encounter[-length(encounter)] - encounter[-1]) > 0))) |>
#   dplyr::select(-cam_intersects) |>
#   ungroup()


  ########################################
  # Collect stay time data for each encounter
  ########################################
  stay_time_raw_tele <- cell_captures_tele |>
    dplyr::group_by(pass_i) |>
    dplyr::summarise(t_stay = max(t) - min(t),
                     XY_inds = XY_inds[1],
                     speed = lscape_type[1],
                     .groups = 'drop')
  
  
  # Format staying time data
  stay_time_tele <- cell_captures_tele |>
    dplyr::group_by(pass_i) |>
    dplyr::summarise(t_stay = max(t) - min(t),
                     XY_inds = XY_inds[1]) |>
    # dplyr::add_row(XY_inds = cam.samps[cam.samps %notin% cam_captures_tele$XY_inds]) |>
    # dplyr::left_join(tri_cam_samps |> dplyr::filter(vertex == 1), # need to define tri_cam_samps with lists...
    #                  by = "XY_inds") |>
    dplyr::select(XY_inds, t_stay) |>
    dplyr::group_by(XY_inds) |>
    dplyr::mutate(encounter = 1:n()) |>
    dplyr::ungroup() |>
    pivot_wider(names_from = encounter, values_from = t_stay)
  
  stay_time_tele <- as.matrix(dplyr::left_join(data.frame(XY_inds = cam.samps),
                                               stay_time_tele,
                                               by = "XY_inds") |>
                                dplyr::select(-XY_inds)
  )
  
  
########################################
# Collect data for triangle viewshed cameras
########################################
# First, collect data within each grid cell to reduce subsequent calculations
# Creates IDs (pass_i) for each pass through a cell and landscape indices for the cell that the individual passes into (next_i)
cell_captures <- animalxy.all |>
  ungroup() |> 
  dplyr::filter(lscape_index %in% cam.samps) |>
  dplyr::mutate(pass_i = cumsum(c(1, (abs(ii[-length(ii)] - ii[-1]) > 1) +
                                    (abs(Animal_ID[-length(Animal_ID)] - Animal_ID[-1]) > 0) + 
                                    (abs(lscape_index[-length(Animal_ID)] - lscape_index[-1]) > 0)))) |> 
  dplyr::group_by(pass_i) |>
  dplyr::mutate(next_i = ifelse(max(t) == t.steps,
                                max(ii),
                                max(ii) + 1)) |>
  dplyr::ungroup()

# # Use animal's next step to capture all trajectories within cell
cell_captures <- cell_captures |>
  dplyr::add_row(animalxy.all[unique(cell_captures$next_i), ] |>
                   dplyr::ungroup() |>
                   dplyr::mutate(pass_i = unique(cell_captures$pass_i))) |>
  dplyr::arrange(ii) |>
  dplyr::select(-next_i)
# Remove repeat rows
cell_captures <- dplyr::distinct(cell_captures)

cell_check <- cell_captures |>
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

# # Separate passes when individuals cross adjacent cells with cameras in both cells
# cell_captures_temp <- cell_captures |>
#   dplyr::group_by(pass_i) |>
#   dplyr::mutate(num_lscapes = length(unique(lscape_index)),
#                 offset = c(0,cumsum(diff(lscape_index) != 0))) |>
#   dplyr::filter(num_lscapes > 2) |>
#   dplyr::group_by(pass_i, offset) |>
#   dplyr::mutate(next_i = ifelse(t == t.steps,
#                                 max(ii),
#                                 max(ii) + 1)) |>
#   do(dplyr::add_row(., animalxy.all[animalxy.all$ii == .$next_i[1],] |>
#                       dplyr::mutate(pass_i = .$pass_i[1],
#                                     offset = .$offset[1]))) |>
#   dplyr::ungroup() |>
#   dplyr::rowwise() |>
#   dplyr::mutate(pass_i_OG = pass_i,
#                 pass_i = pass_i + offset/100) |>
#   dplyr::group_by(pass_i_OG) |>
#   dplyr::filter(pass_i < max(pass_i)) |>
#   dplyr::ungroup() |>
#   dplyr::select(-c(num_lscapes, offset, next_i))

# # Include multi-cell passes to original cell capture data
# cell_captures <- cell_captures |>
#   dplyr::filter(pass_i %notin% cell_captures_temp$pass_i_OG)

# cell_captures <- cell_captures |>
#   dplyr::ungroup() |>
#   dplyr::add_row(cell_captures_temp |>
#                    dplyr::select(-pass_i_OG)) |>
#   dplyr::arrange(ii) |>
#   dplyr::group_by(pass_i) |>
#   dplyr::mutate(lscape_index = lscape_index[1]) |>
#   dplyr::ungroup()

# Collect all camera viewshed captures
# in_cam determines whether an individual starts within the camera viewshed
cam_captures <- cell_captures |>
  dplyr::group_by(pass_i) |>
  dplyr::summarise(cam_samps = list(tri_cam_samps[tri_cam_samps$lscape_index %in% lscape_index[1], 3:4]),
                   cam_intersects = list(calc_intersects(matrix(unlist(cam_samps), 
                                                                nrow = length(unlist(cam_samps))/2, ncol = 2),
                                                         t(rbind(x, y)),
                                                         speeds,
                                                         t)),
                   # intersect = c(lapply(cam_intersects, '[[', 1), NA),
                   # t_stay = sum(unique(unlist(lapply(cam_intersects, '[[', 2)))),
                   t_stay = c(unlist(lapply(cam_intersects, '[[', 2)), NA),
                   in_cam = c(unlist(lapply(cam_intersects, '[[', 3)), NA),
                   encounter = unlist(lapply(cam_intersects, '[[', 4)),
                   t_in = c(unlist(lapply(cam_intersects, '[[', 5)), NA),
                   t_out = c(unlist(lapply(cam_intersects, '[[', 6)), NA),
                   lscape_index = lscape_index[1],
                   Animal_ID = Animal_ID,
                   t = t,
                   ii = ii,
                   x = x,
                   y = y,
                   lscape_type = lscape_type[1],
                   road = road,
                   .groups = 'drop'
  )  |>
  dplyr::filter(t_stay > 0) |>
  dplyr::group_by(pass_i) |>
  dplyr::mutate(pass_i = pass_i + 0.001 * cumsum(c(1, abs(encounter[-length(encounter)] - encounter[-1]) > 0))) |>
  dplyr::select(-cam_intersects) |>
  ungroup()

########################################
# Snapshot count data
########################################
count_data <- cam_captures |>
  dplyr::filter(t %in% 0:(t.steps - 1) & in_cam == T) |> # Take `snapshots' of data
  # dplyr::rowwise() |>
  # dplyr::mutate(in_tri = calc_in_triangle(tri_cam_samps[tri_cam_samps$lscape_index %in% lscape_index, 3:4],
  #                                         t(rbind(x, y)))) |>
  # dplyr::filter(in_tri == TRUE) |>
  dplyr::group_by(lscape_index, t) |>
  dplyr::summarise(clump = n(),
                   x = x,
                   y = y,
                   .groups= 'drop') |>
  dplyr::group_by(lscape_index) |>
  dplyr::summarise(count = sum(clump),
                   mean_clump = mean(clump),
                   x = x[1],
                   y = y[1],
                   .groups = 'drop')

# Add zero count data for all cameras without counts
count_data <- dplyr::add_row(count_data,
                             lscape_index = cam.samps[cam.samps %notin% count_data$lscape_index],
                             count = 0)

# Change order of count data
count_data <- dplyr::left_join(data.frame(lscape_index = cam.samps),
                               count_data,
                               by = c("lscape_index")) |>
  dplyr::mutate(speed = lscape_speeds$Speed[lscape_index],
                road = lscape_speeds$Road[lscape_index])



# Only collect other data when counts are available
if(max(count_data$count) > 0) {
  
  ########################################
  # Collect observations at each camera
  ########################################
  # # Track the number of times an individual encounters each camera and time spent in front of the camera
  # cam_data <- cam_captures |>
  #   dplyr::filter(lscape_index %in% cam.samps) |>
  #   dplyr::group_by(pass_i) |>
  #   dplyr::summarise(t_stay = unique(t_stay),
  #                    lscape_index = unique(lscape_index),
  #                    .groups = 'drop') |>
  #   dplyr::ungroup() |>
  #   dplyr::group_by(lscape_index) |>
  #   dplyr::mutate(encounter = length(unique(pass_i))) |>
  #   dplyr::ungroup()
  # 
  # cam_data <- dplyr::add_row(cam_data,
  #                            lscape_index = cam.samps[cam.samps %notin% cam_data$lscape_index],
  #                            encounter = 0) |>
  #   dplyr::mutate(speed = lscape_speeds$Speed[lscape_index],
  #                 road = lscape_speeds$Road[lscape_index])
  
  ########################################
  # Collect encounter data for each camera
  ########################################
  encounter_data <- cam_captures |>
    dplyr::group_by(lscape_index) |>
    dplyr::summarise(encounter = length(unique(pass_i)),
                     .groups = 'drop')
  
  
  # Change order of encounter data
  encounter_data <- dplyr::left_join(data.frame(lscape_index = cam.samps),
                                     encounter_data,
                                     by = "lscape_index") |>
    dplyr::mutate(speed = lscape_speeds$Speed[lscape_index],
                  road = lscape_speeds$Road[lscape_index])
  
  encounter_data$encounter[is.na(encounter_data$encounter)] <- 0
  
  # encounter_data <- cam_data |>
  #   dplyr::group_by(lscape_index) |>
  #   dplyr::summarise(encounter = unique(encounter),
  #                    speed = unique(speed),
  #                    road = unique(road),
  #                    .groups = 'drop')
  # 
  # # Change order of encounter data
  # encounter_data <- dplyr::left_join(data.frame(lscape_index = cam.samps),
  #                                    encounter_data,
  #                                    by = "lscape_index"
  # )
  
  ########################################
  # Collect stay time data for each encounter
  ########################################
  stay_time_raw <- cam_captures |>
    dplyr::group_by(pass_i) |>
    dplyr::summarise(t_stay = sum(t_stay),
                     lscape_index = lscape_index[1],
                     speed = lscape_type[1],
                     .groups = 'drop')
  

  # Format staying time data
  stay_time_data <- cam_captures |>
    dplyr::group_by(pass_i) |>
    dplyr::summarise(t_stay = sum(t_stay),
                     lscape_index = lscape_index[1]) |>
    dplyr::add_row(lscape_index = cam.samps[cam.samps %notin% cam_captures$lscape_index]) |>
    dplyr::left_join(tri_cam_samps |> dplyr::filter(vertex == 1), # need to define tri_cam_samps with lists...
                     by = "lscape_index") |>
    dplyr::select(lscape_index, t_stay) |>
    dplyr::group_by(lscape_index) |>
    dplyr::mutate(encounter = 1:n()) |>
    dplyr::ungroup() |>
    pivot_wider(names_from = encounter, values_from = t_stay)
  
  stay_time_data <- as.matrix(dplyr::left_join(data.frame(lscape_index = cam.samps),
                                               stay_time_data,
                                               by = "lscape_index") |>
                                dplyr::select(-lscape_index)
  )
  
  
  ########################################
  # Collect TTE at each camera
  ########################################
  TTE_data <- data.frame(lscape_index = double(),
                         TTE = double(),
                         occasion = double())
  
  # Loop through all occasions
  for(jj in 1:num.occ){
    occ.times <- t.steps*c((jj-1)/num.occ, jj/num.occ)
    TTE_temp <- cam_captures |>
      dplyr::filter(lscape_index %in% cam.samps,
                    t >= occ.times[1],
                    t < occ.times[2])
    if(nrow(TTE_temp) > 0) {
      TTE_data <- rbind(TTE_data, TTE_temp |>
                          dplyr::group_by(lscape_index) |>
                          dplyr::summarise(TTE = min(t) - (jj-1)*JJ,
                                           occasion = jj))
    }
    TTE_data <- dplyr::add_row(TTE_data,
                               lscape_index = cam.samps[cam.samps %notin% TTE_temp$lscape_index],
                               TTE = NA,
                               occasion = jj)
  }
  TTE_data <- TTE_data |>
    dplyr::mutate(speed = lscape_speeds$Speed[lscape_index],
                  road = lscape_speeds$Road[lscape_index])
  
  TTE_data_raw <- TTE_data
  
  # Format TTE data
  TTE_data <- TTE_data |>
    dplyr::add_row(lscape_index = cam.samps[cam.samps %notin% TTE_data$lscape_index]) |>
    dplyr::left_join(tri_cam_samps,
                     by = "lscape_index") |>
    dplyr::select(lscape_index, TTE, occasion) |>
    dplyr::distinct() |>
    pivot_wider(names_from = occasion, values_from = TTE)
  
  TTE_data <- as.matrix(dplyr::left_join(data.frame(lscape_index = cam.samps),
                                         TTE_data,
                                         by = "lscape_index") |>
                          dplyr::select(-lscape_index)
  )
  TTE_data[TTE_data == 0] <- NA
  
  
  ########################################
  # Collect STE data
  ########################################
  # STE_data_long <- cell_captures |>
  #   dplyr::filter(t %in% 0:(t.steps - 1) & lscape_index %in% cam.samps) |>
  #   dplyr::rowwise() |>
  #   dplyr::mutate(in_tri = calc_in_triangle(tri_cam_samps[tri_cam_samps$lscape_index %in% lscape_index, 3:4],
  #                                           t(rbind(x, y)))) |>
  #   dplyr::filter(in_tri == TRUE) |>
  #   dplyr::select(Animal_ID, t, lscape_index) |>
  #   dplyr::group_by(t) |>
  #   dplyr::summarise(cam_inds = sample(c(unique(lscape_index), rep(0, ncam - length(unique(lscape_index))))),
  #                    .groups = 'drop')
  # 
  # STE_data <- STE_data_long |>
  #   dplyr::group_by(t) |>
  #   dplyr::summarise(STE = min(which(cam_inds > 0))*cam.A)
  # STE_data <- dplyr::add_row(STE_data, STE = rep(NA, t.steps - nrow(STE_data)))
  
  STE_data <- cam_captures |>
    dplyr::filter(t %in% 0:(t.steps - 1) & in_cam == T) |> # Take `snapshots' of data
    dplyr::select(Animal_ID, t, lscape_index) |>
    dplyr::group_by(t) |>
    dplyr::summarise(cam_inds = list(sample(c(unique(lscape_index), rep(0, ncam - length(unique(lscape_index)))))),
                     STE = min(which(unlist(cam_inds) > 0))*cam.A,
                     .groups = 'drop') |>
    dplyr::select(t, STE) 
  STE_data <- dplyr::add_row(STE_data, 
                             STE = rep(NA, t.steps - nrow(STE_data))) |> 
    dplyr::arrange(t)
    
  
  # Mean clump size for adjusted STE
  mean_clump <- animalxy.all |>
    dplyr::filter(t %in% 0:(t.steps - 1) & lscape_index %in% cam.samps) |> # Take `snapshots' of data
    dplyr::group_by(lscape_index, t) |>
    dplyr::summarise(clump = length(t),
                     .groups= 'drop') |>
    dplyr::summarise(mean_clump = mean(clump)) |>
    pull(mean_clump)
  
  ########################################
  # Collect space-use for each animal
  ########################################
  animal_data <- animalxy.all |>
    dplyr::group_by(Animal_ID) |>
    dplyr::mutate(t_diff = c(t[2:length(t)] - t[1:(length(t)-1)], 0)) |>
    dplyr::group_by(lscape_type) %>% 
    dplyr::summarise(t_spent = sum(t_diff)/nind,
                     .groups = 'drop') 
  
}

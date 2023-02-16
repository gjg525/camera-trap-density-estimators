########################################
# Collect data for triangle viewshed cameras
########################################
# First, collect data within each grid cell to reduce calculations later on
cell_captures <- animalxy.all |>
  dplyr::filter(lscape_index %in% cam.samps) |>
  dplyr::ungroup() |>
  dplyr::mutate(pass_i = cumsum(c(1, abs(ii[-length(ii)] - ii[-1]) > 1))) |>
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

# Split passes with animals spanning multiple cells
cell_captures_temp <- cell_captures |>
  dplyr::group_by(pass_i) |>
  dplyr::mutate(num_lscapes = length(unique(lscape_index)),
         offset = c(0,cumsum(diff(lscape_index) != 0))) |>
  dplyr::filter(num_lscapes > 2) |>
  dplyr::group_by(pass_i, offset) |>
  dplyr::mutate(next_i = ifelse(t == t.steps,
                         max(ii),
                         max(ii) + 1)) |>
  do(dplyr::add_row(., animalxy.all[animalxy.all$ii == .$next_i[1],] |>
               dplyr::mutate(pass_i = .$pass_i[1],
                      offset = .$offset[1]))) |>
  dplyr::ungroup() |>
  dplyr::rowwise() |>
  dplyr::mutate(pass_i_OG = pass_i,
         pass_i = pass_i + offset/100) |>
  dplyr::group_by(pass_i_OG) |>
  dplyr::filter(pass_i < max(pass_i)) |>
  dplyr::ungroup() |>
  dplyr::select(-c(num_lscapes, offset, next_i))

# Include multi-cell passes to original cell capture data
cell_captures <- cell_captures |>
  dplyr::filter(pass_i %notin% cell_captures_temp$pass_i_OG)

cell_captures <- cell_captures |>
  dplyr::ungroup() |>
  dplyr::add_row(cell_captures_temp |>
            dplyr::select(-pass_i_OG)) |>
  dplyr::arrange(ii) |>
  dplyr::group_by(pass_i) |>
  dplyr::mutate(lscape_index = lscape_index[1]) |>
  dplyr::ungroup()

# Collect all camera viewshed captures
cam_captures <- cell_captures |>
  dplyr::group_by(pass_i) |>
  dplyr::summarise(cam_samps = list(tri_cam_samps[tri_cam_samps$lscape_index %in% lscape_index[1], 3:4]),
           cam_intersects = list(calc_intersects(matrix(unlist(cam_samps), nrow = length(unlist(cam_samps))/2, ncol = 2),
                                            t(rbind(x, y)),
                                            speeds)),
           intersect = lapply(cam_intersects, '[[', 1),
           t_stay = sum(unique(unlist(lapply(cam_intersects, '[[', 2)))),
           lscape_index = lscape_index,
           Animal_ID = Animal_ID,
           t = t,
           ii = ii,
           x = x,
           y = y,
           lscape_type = lscape_type,
           road = road,
           .groups = 'drop'
            )  |>
  dplyr::filter(t_stay > 0) |>
  dplyr::select(-cam_intersects)

########################################
# Snapshot count data
########################################
count_data <- cam_captures |>
  dplyr::filter(t %in% 1:t.steps) |> # Take `snapshots' of data
  dplyr::rowwise() |>
  dplyr::mutate(in_tri = calc_in_triangle(tri_cam_samps[tri_cam_samps$lscape_index %in% lscape_index, 3:4],
                                   t(rbind(x, y)))) |>
  dplyr::filter(in_tri == TRUE) |>
  dplyr::group_by(lscape_index, t) |>
  dplyr::summarise(clump = length(t),
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
                      count = 0) |>
  dplyr::mutate(speed = lscape_speeds$Speed[lscape_index],
         road = lscape_speeds$Road[lscape_index])

# Change order of encounter data
count_data <- dplyr::left_join(data.frame(lscape_index = cam.samps),
                        count_data,
                        by = "lscape_index"
)

# Only collect other data when counts are available
if(max(count_data$count) > 0) {

########################################
# Collect observations at each camera
########################################
# Track the number of times an individual encounters each camera and time spent in front of the camera
cam_data <- cam_captures |>
  dplyr::filter(lscape_index %in% cam.samps) |>
  dplyr::group_by(Animal_ID, pass_i) |>
  dplyr::summarise(t_stay = unique(t_stay),
            lscape_index = unique(lscape_index),
            .groups = 'drop') |>
  dplyr::ungroup() |>
  dplyr::group_by(lscape_index) |>
  dplyr::mutate(encounter = length(unique(pass_i))) |>
  dplyr::ungroup()

cam_data <- dplyr::add_row(cam_data,
                    lscape_index = cam.samps[cam.samps %notin% cam_data$lscape_index],
                  encounter = 0) |>
  dplyr::mutate(speed = lscape_speeds$Speed[lscape_index],
       road = lscape_speeds$Road[lscape_index])

########################################
# Collect encounter data for each camera
########################################
encounter_data <- cam_data |>
  dplyr::group_by(lscape_index) |>
  dplyr::summarise(encounter = unique(encounter),
            speed = unique(speed),
            road = unique(road),
            .groups = 'drop')

# Change order of encounter data
encounter_data <- dplyr::left_join(data.frame(lscape_index = cam.samps),
                            encounter_data,
                            by = "lscape_index"
                            )

########################################
# Collect stay time data for each encounter
########################################
stay_time_raw <- cam_data[!is.na(cam_data$t_stay), ]

# Format staying time data
stay_time_data <- stay_time_raw |>
  dplyr::add_row(lscape_index = cam.samps[cam.samps %notin% stay_time_raw$lscape_index]) |>
  dplyr::left_join(tri_cam_samps) |>
  dplyr::select(lscape_index, t_stay) |>
  dplyr::distinct() |>
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
                  t <= occ.times[2])
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
  dplyr::left_join(tri_cam_samps) |>
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
STE_data_long <- cell_captures |>
  dplyr::filter(t %in% 1:t.steps & lscape_index %in% cam.samps) |>
  dplyr::rowwise() |>
  dplyr::mutate(in_tri = calc_in_triangle(tri_cam_samps[tri_cam_samps$lscape_index %in% lscape_index, 3:4],
                                   t(rbind(x, y)))) |>
  dplyr::filter(in_tri == TRUE) |>
  dplyr::select(Animal_ID, t, lscape_index) |>
  dplyr::group_by(t) |>
  dplyr::summarise(cam_inds = sample(c(unique(lscape_index), rep(0, ncam - length(unique(lscape_index))))),
            .groups = 'drop')

STE_data <- STE_data_long |>
  dplyr::group_by(t) |>
  dplyr::summarise(STE = min(which(cam_inds > 0))*cam.A)
STE_data <- dplyr::add_row(STE_data, STE = rep(NA, t.steps - nrow(STE_data)))

# Mean clump size for adjusted STE
mean_clump <- animalxy.all |>
  dplyr::filter(t %in% 1:t.steps & lscape_index %in% cam.samps) |> # Take `snapshots' of data
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
  dplyr::mutate(t_diff = c(0,t[2:length(t)] - t[1:(length(t)-1)])) |>
  dplyr::group_by(Animal_ID, road) |>
  dplyr::summarise(t_spent = sum(t_diff),
            .groups = 'drop') |>
  dplyr::mutate(area_cover = ifelse(road == "Off Trail",
                           (sum(lscape_speeds$Road == "Off Trail", na.rm = T)),
                           sum(lscape_speeds$Road == "On Trail", na.rm = T)))

}

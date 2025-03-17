library(tidyverse)
library(RColorBrewer)
library(lattice)
library(gridExtra)
library(doParallel)

source("./R/EDE_update/utils.R")
source("./R/EDE_update/plot_funs.R")
source("./R/EDE_update/ABM_sim.R")
source("./R/EDE_update/MCMC_functions.R")
source("./R/EDE_update/Create_landscape.R")
source("./R/EDE_update/Collect_data.R")
source("./R/EDE_update/Collect_tele_data.R")

# # Load animal GPS data
# load(file = "Sim_results/save_animal_data.RData")
# load(file = "Sim_results/save_lscape_defs.RData")

# Initializations
fig_colors <- c("#2ca25f", "#fc8d59", "#67a9cf", "#f768a1", "#bae4b3", "#fed98e")
options(ggplot2.discrete.colour = fig_colors)
options(ggplot2.discrete.fill = fig_colors)

# Study design
study_design <- tibble::tibble(
  q = 30^2, # Number grid cells
  dx = 1,  # Grid cell lengths
  dy = 1,
  t_steps = 500, # Number of time steps (3000 * 15 min ~ 31.25 days)
  dt = 1,     # Time step size
  bounds = list(c(0, dx * q ^ 0.5)), # Sampling area boundaries
  tot_A = (bounds[[1]][2] - bounds[[1]][1])^2,
  num_groups = 50,
  # num_groups = 25,
  group_sizes = list(rep(1, num_groups)),
  group_spread = 0, # Tightness of grouping behavior (relative to grid size)
  tot_animals = sum(unlist(group_sizes)),
  h_range_strength = NA,  # Home range size (not to scale)
  activity_sync = "sync",
  activity_prob = list(rep(1, t_steps)), # can be defined for all time steps
  # MCMC parms
  num_runs = 1000,
  n_iter = 40000,
  burn_in = 30000,
  # staying time censors
  t_censor = 10, 
  run_models = list(1:5),
  covariate_labels = list(c("On Trail", "Off Trail"))
)

# Landscape design
lscape_design <- tibble::tibble(
  lscape_tag = "grid", # "Random", #  "Random road", # "cross", # 
  default_kappa = 0, # Turning angle
  num_roads = 2,
  Speed_ID = list(c("Medium")),
  # Speed_mins = list(c(.05, .3, .9)),
  # Speed_maxes = list(c(.15, .5, 1.1)),
  # Speed_mins = list(c(.05)),
  # Speed_maxes = list(c(.15)),
  Speed_mins = list(c(.3)),
  Speed_maxes = list(c(.5)),
  Trail_ID = list(c("On Trail", "Off Trail")),
  Trail_speed = list(c("Medium", "Medium"))
)

# Cam designs
cam_design <- tibble::tibble(
  ncam = 250,
  Design = "Random",
  Props = list(c(1, 1)), # proportion of cameras placed on and off roads
  # Design = "Road Bias",
  # Props = list(c(1, 0)), # proportion of cameras placed on and off roads
  # Props = list(c(0, 1)), # proportion of cameras placed on and off roads
  cam_length = study_design$dx * 0.9, # length of all viewshed sides
  # cam_length = study_design$dx * 0.1, # length of all viewshed sides
  cam_A = cam_length ^ 2 / 2,
  tot_snaps = ncam * study_design$t_steps
)

num_models <- length(unlist(study_design$run_models))

# Initialize summary matrices
D.all <- tibble::tibble(
  iteration = rep(1:study_design$num_runs, each = num_models),
  Model = NA,
  Covariate = NA,
  Est = NA,
  SD = NA
  # all_results = NA
)

all_data <- tibble::tibble(
  iteration = 1:study_design$num_runs,
  cam_captures = NA,
  count_data = NA,
  encounter_data = NA,
  stay_time_all = NA,
  stay_time_data = NA
)

# # camera data summaries
cam_design <- cam_design %>%
  dplyr::mutate(percent_cam_coverage = ncam * cam_A / study_design$tot_A)

IS_mean <- c()
cell_mean <- c()
cell_mean_PR <- c()
PR_hab <- c()

# Multi-run simulations
for (run in 1:study_design$num_runs) {
  print(paste("Run", run, "of", study_design$num_runs))
  
  # # Create custom landscape (tag = "grid", "circ", "squares", "metapop")
  lscape_defs <- lscape_creator(study_design, lscape_design)
  
  # border_check <- border_check(
  #   lscape_defs, 
  #   max(lscape_defs$X), 
  #   max(lscape_defs$Y)
  # )
  # 
  # if(!is.null(border_check)) {
  #   warning("Check directions on borders")
  # }
  
  # # Run agent-based model
  animalxy.all <- ABM_sim(study_design,
                          lscape_defs)
  
  # Create covariate matrix with 0, 1 values
  study_design <- study_design %>% 
    dplyr::mutate(
      num_covariates = length(unlist(covariate_labels)),
      Z = list(create_covariate_mat(
        lscape_defs,
        study_design,
        unlist(covariate_labels),
        grouping = "Road")
      )
    )
  
  tele_summary <- Collect_tele_data(
    animalxy.all, 
    study_design, 
    grouping = "Road"
  )
  
  # Use smallest stay time as reference category
  # ref_cat_idx <- which(tele_summary$stay_prop == min(tele_summary$stay_prop))
  ref_cat_idx <- which(tele_summary$stay_prop == min(tele_summary$stay_prop))
  
  # Set reference category for intercept
  study_design$Z[[1]][, ref_cat_idx] <- 1
  
  # Subtract reference category from stay time proportion
  prop_adjust <- tele_summary$stay_prop /
    tele_summary$stay_prop[ref_cat_idx]
  prop_adjust[ref_cat_idx] <- tele_summary$stay_prop[ref_cat_idx]
  kappa.prior.mu.adj <- log(prop_adjust)
  kappa.prior.mu <- log(tele_summary$stay_prop)
  kappa.prior.var <- tele_summary$stay_sd^2 # stay_time_summary$cell_sd ^ 2
  
  # Place cameras on study area
  cam_locs <- create_cam_samp_design(study_design,
                                     lscape_defs,
                                     cam_design)
  
  ################################
  # Collect data
  ################################
  "%notin%" <- Negate("%in%")
  all_data$cam_captures[run] <- list(get_cam_captures(animalxy.all %>%
                                                        dplyr::filter(t != 0)))
  
  # habitat_summary <- lscape_defs %>%
  #   dplyr::group_by(Road) %>%
  #   dplyr::summarise(
  #     n_lscape = dplyr::n()
  #   ) %>%
  #   dplyr::ungroup() %>%
  #   dplyr::left_join(
  #     cam_locs %>%
  #       dplyr::group_by(Road) %>%
  #       dplyr::summarise(
  #         ncams = dplyr::n(),
  #         .groups = 'drop'
  #       ),
  #     by = dplyr::join_by(Road)
  #   ) %>%
  #   dplyr::mutate(
  #     prop_cams = ncams / sum(ncams, na.rm = T)
  #   ) %>%
  #   replace(is.na(.), 0) %>% 
  #   dplyr::select(Road, n_lscape, prop_cams)
  # 
  # habitat_summary <- habitat_summary[order(unlist(study_design$covariate_labels)),]
  
  # Run models only if any data points were collected
  all_data[run, ] <- all_data[run, ] %>% 
    dplyr::mutate(
      count_data = list(get_count_data(cam_locs, cam_captures[[1]], animalxy.all %>% 
                                         dplyr::filter(t != 0))),
      encounter_data = list(get_encounter_data(cam_locs, cam_captures[[1]])),
      stay_time_all = list(get_stay_time_data(cam_locs, cam_captures[[1]])),
      # stay_time_raw = list(stay_time_all[[1]][[1]]),
      stay_time_data = list(stay_time_all[[1]][[2]])
    )
  
  count_data <- all_data$count_data[[run]]

  habitat_summary <- lscape_defs %>%
    dplyr::group_by(Road) %>%
    dplyr::summarise(
      n_lscape = dplyr::n()
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      prop_lscape = n_lscape / sum(n_lscape)
    ) %>%
    dplyr::left_join(
      cam_locs %>%
        dplyr::group_by(Road) %>%
        dplyr::summarise(
          ncams = dplyr::n(),
          .groups = 'drop'
        ),
      by = dplyr::join_by(Road)
    ) %>%
    dplyr::left_join(
      tele_summary, #%>%
      # dplyr::mutate(
      #   tot_n = sum(stay_prop),
      #   stay_prop = stay_prop / tot_n
      # ),
      by = dplyr::join_by(Road)
    )  %>%
    dplyr::mutate(
      prop_cams = ncams / sum(ncams, na.rm = T),
      d_coeff = n_lscape * prop_cams / stay_prop / (cam_design$cam_A * study_design$t_steps),
      cell_coeff = n_lscape * prop_lscape / stay_prop / (study_design$t_steps)
    ) %>%
    dplyr::arrange(desc(Road)) %>%
    replace(is.na(.), 0)
  
  # cc <- get_count_data(cam_locs, all_data$cam_captures[[run]], animalxy.all %>%
  #                        dplyr::filter(t != 0))
  
  n_adj <- habitat_summary %>%
    dplyr::left_join(
      count_data %>%
        dplyr::group_by(Road) %>%
        dplyr::summarise(
          mean_count = mean(count, na.rm = T),
          .groups = 'drop'
        ),
      by = dplyr::join_by(Road)
    ) %>%
    dplyr::mutate(
      mean_count = tidyr::replace_na(mean_count, 0),
      ncams = tidyr::replace_na(ncams, 0)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      # stay_prop_adj = 1 / stay_prop * ncams / sum(ncams, na.rm = T),
      n_lscape = mean_count * n_lscape / (cam_design$cam_A * study_design$t_steps),
      n_habitat = mean_count * d_coeff
      # n_full = n_lscape / stay_prop
      # big_D = mean_count * sum(stay_prop) / stay_prop
    ) %>%
    dplyr::select(Road, prop_lscape, mean_count, ncams, stay_prop, prop_cams, d_coeff, n_lscape, n_habitat, cell_coeff)
  
  # This method is more reliable when one or more lscape type is missing
  PR_hab[run] <- sum(n_adj$n_habitat)
  
  IS_mean[run] <- mean(count_data$count) * study_design$tot_A / 
    (cam_design$cam_A * study_design$t_steps)
  
  cell_mean[run] <- animalxy.all %>% 
    dplyr::filter(t %in% 1:500) %>% 
    dplyr::group_by(lscape_index, Road) %>% 
    dplyr::count() %>% 
    dplyr::ungroup() %>%
    # dplyr::group_by(Road) %>% 
    dplyr::summarise(mm = mean(n)) %>%  
    dplyr::mutate(tot_N = mm * study_design$tot_A / (study_design$t_steps)) %>% 
    dplyr::pull(tot_N)
  
  cell_mean_PR[run] <- animalxy.all %>% 
    dplyr::filter(t %in% 1:500) %>% 
    dplyr::group_by(lscape_index, Road) %>% 
    dplyr::count() %>% 
    # dplyr::ungroup() %>%
    dplyr::group_by(Road) %>%
    dplyr::summarise(mm = mean(n)) %>% 
    dplyr::left_join(habitat_summary %>% 
                       dplyr::select(Road, cell_coeff),
                     by = dplyr::join_by(Road)) %>% 
    dplyr::mutate(tot_N = mm * cell_coeff) %>% 
    dplyr::ungroup() %>% 
    dplyr::summarise(tot_N = sum(tot_N)) %>% 
    dplyr::pull(tot_N)
  
  
}


# # # Remove outlier estimates
D.all$Est[D.all$Est > 5 * study_design$tot_animals] <- NA
D.all$SD[D.all$SD > 5 * study_design$tot_animals] <- NA

# # Omit all_results column (too much data)
# D.all <- D.all %>% 
#   dplyr::select(-all_results)  
#   
# # Omit TTE data
# all_data <- all_data %>% 
#   dplyr::select(-c(stay_time_all, TTE_data_all, TTE_data_raw, TTE_data))
# 
# D.all$Model <- factor(D.all$Model, levels = c("TDST", "REST", "TTE", "PR", "IS", "STE", "SECR"))
# 
# NA_summary <- D.all %>%
#   group_by(Model) %>%
#   summarise(num_NAs = sum(is.na(Est)))

###########################
# Plots
###########################
plot_multirun_means(study_design, D.all %>% 
                      dplyr::filter(is.finite(Est)))
plot_multirun_sds(D.all %>% 
                    dplyr::filter(is.finite(Est)))
# plot_multirun_mape(D.all %>%
#                       dplyr::filter(is.finite(Est)),
#                    study_design$tot_animals)
# plot_multirun_CV(D.all %>% 
#                    dplyr::filter(is.finite(Est)))
# plot_multirun_hist(D.all)

plot_ABM(study_design,
         cam_design,
         cam_locs,
         animalxy.all)
plot_space_use(study_design,
               animalxy.all)
#
plot_count_data(count_data_in = all_data$count_data[[1]],
                fill = "Road")
# plot_encounter_data(encounter_data_in = all_data$encounter_data[[1]],
#                 fill = "Speed")
# plot_staytime_data(stay_time_raw_in = all_data$stay_time_raw[[1]],
#                    fill = "Speed")

lscape_defs %>% 
  dplyr::group_by(Road) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(props = n / sum(n))

animalxy.all %>% 
  dplyr::filter(t %in% 1:500) %>% 
  dplyr::group_by(Road) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(props = n / sum(n))

animalxy.all %>% 
  dplyr::filter(t %in% 1:500) %>% 
  dplyr::group_by(Road, lscape_index) %>%
  dplyr::count() %>% 
  dplyr::group_by(Road) %>% 
  dplyr::summarise(mm = mean(n)) %>% 
  dplyr:::ungroup() %>% 
  dplyr::mutate(props = mm / sum(mm))


count_data %>% 
  dplyr::group_by(Road) %>% 
  dplyr::summarise(tot_count = sum(count)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(props = tot_count / sum(tot_count))

animalxy.all %>% 
  dplyr::filter(t %in% 1:500) %>% 
  dplyr::group_by(Road, lscape_index) %>%
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  ggplot(aes(x = n, fill = Road)) + 
  geom_histogram(position = "identity", alpha = 0.4)


# # Plot camera time series
# cam_caps <- all_data$cam_captures[[1]]
# cam_caps <- cam_caps %>%
#   dplyr::filter(in_cam == T) %>%
#     dplyr::mutate(t_hour = ceiling(t)) %>%
#   dplyr::group_by(t_hour) %>%
#     dplyr::arrange(t_hour) %>% 
#   dplyr::summarise(
#     group_size = n(),
#     .groups = "drop"
#   )
# 
# cam_caps <- cam_caps %>%
#   dplyr::bind_rows(tibble::tibble(
#     t_hour = which(seq(study_design$dt, 
#                        study_design$t_steps * study_design$dt,
#                        by = study_design$dt) %notin% cam_caps$t_hour),
#       group_size = 0
#     ) ) %>% 
#   arrange(t_hour)
# 
# xy_cams <- unnest(cam_locs, cols = c(x, y, vertex)) %>%
#   group_by(lscape_index) %>%
#   dplyr::mutate(
#     x = x[1],
#     y = y[1]
#   ) %>%
#   dplyr::select(cam_ID, lscape_index, x, y) %>%
#   distinct()

# # Plot spatial cam captures
# cd <- all_data$count_data[[1]]
# 
# cd %>%
#   ggplot() +
#   geom_point(aes(xy_cams$x, xy_cams$y),
#              color = "black",
#              # shape = 2,
#              size = cd$count*5,
#              stroke = 2) +
#   labs(x = "x", y = "y") +
#   theme_minimal(base_size = 25)
# 
# animal_stats <- animalxy.all %>% 
#   dplyr::group_by(Animal_ID) %>% 
#   dplyr::summarise(
#     tot_dist = sum(trav_dist, na.rm = T)
#   )
# 
Data_summary <- all_data %>%
  group_by(iteration) %>%
  dplyr::summarise(
    num_counts = sum(count_data[[1]]$count),
    num_det = cam_captures[[1]] %>%
      dplyr::filter(in_cam == T) %>%
      dplyr::group_by(lscape_index, t) %>%
      dplyr::summarise(
        group = n(),
        .groups = "drop"
      ) %>%
      dplyr::summarise(num_det = sum(group > 0)) %>%
      pull(num_det),
    det_prob = num_det / (study_design$t_steps * cam_design$ncam),
    num_encounters = sum(encounter_data[[1]]$encounter),
    mean_stay_time = mean(as.matrix(stay_time_data[[1]]), na.rm = T),
    max_stay_time = max(stay_time_data[[1]], na.rm = T),
    num_stay_time_outliers_1 = sum(stay_time_data[[1]] > study_design$dt, na.rm = T),
    num_stay_time_outliers_2 = sum(stay_time_data[[1]] > 2 * study_design$dt, na.rm = T),
    num_stay_time_outliers_3 = sum(stay_time_data[[1]] > 3 * study_design$dt, na.rm = T)
  )


# results_road_slow_cam <- list(
#   # save_animal_data,
#   study_design,
#   cam_design,
#   lscape_design,
#   all_data,
#   D.all
# )
# 
# save(results_road_slow_cam, file = "Sim_results/results_road_slow_cam.RData")


# save(save_animal_data, file = "Sim_results/save_animal_data.RData")
# save(save_lscape_defs, file = "Sim_results/save_lscape_defs.RData")

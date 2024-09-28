# NOTES
# All cams on slow or all cams on fast with abm in for loop: pr habitat works perfectly
# 60% fast cams overestimates (mean = 103)
# 60% slow cams slight underestimate (mean = 99.5)
# random cams slight overestimate (mean = 102.7), and slight underestimate? (mean = 98.9)
# random cams with 500 cams slight underestimate (mean = 99.1)

# # PUT LSCAPE_DESIGN IN FOR LOOP
# No bias in random sample design
# median = 100.8 and mean = 101.2 for 80% fast
# median = 99.7 and mean = 99.5 for 80% slow

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
  num_groups = 100,
  group_sizes = list(rep(1, num_groups)),
  group_spread = 0, # Tightness of grouping behavior (relative to grid size)
  tot_animals = sum(unlist(group_sizes)),
  h_range_strength = NA,  # Home range size (not to scale)
  activity_sync = "sync",
  activity_prob = list(rep(1, t_steps)), # can be defined for all time steps
  # MCMC parms
  num_runs = 1000,
  n_iter = 5000,
  burn_in = 3000,
  # staying time censors
  t_censor = 10, 
  # num_occ = t_steps / TTE_censor,
  run_models = list(1:5),
  covariate_labels = list(c("Slow", "Medium", "Fast"))
)

# Landscape design
lscape_design <- tibble::tibble(
  lscape_tag = "Random",
  default_kappa = 0, # Turning angle
  num_roads = 0,
  Speed_ID = list(c("Slow", "Medium", "Fast")),
  # Speed_maxes = list(c(.15, .5, 1.1)),
  # Speed_mins = list(c(.04, .19, .4)),
  # Speed_maxes = list(c(.05, .21, .6)),
  Speed_mins = list(c(.19, .4, .7)),
  Speed_maxes = list(c(.21, .5, .9)),
  # Speed_mins = list(c(.19, .19, .19)),
  # Speed_maxes = list(c(.21, .21, .21)),
  Trail_ID = list(c("On Trail", "Off Trail")),
  Trail_speed = list(c("Medium", "Medium"))
)

# Cam designs
cam_design <- tibble::tibble(
  ncam = 250,
  # Design = "Random",
  # Props = list(c(1, 1, 1)), # proportion of cameras placed on and off roads
  Design = "Bias",
  Props = list(c(0.6, 0.2, 0.2)), # proportion of cameras placed on and off roads
  # Design = "Bias",
  # Props = list(c(0, 0, 1)), # proportion of cameras placed on and off roads
  cam_length = study_design$dx * 0.3, # length of all viewshed sides
  # cam_length = study_design$dx * 0.005, # length of all viewshed sides (12.5 m^2)
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
  SD = NA,
  all_results = NA
)

save_animal_data <- tibble::tibble(
  iteration = 1:study_design$num_runs,
  data = NA
)

all_data <- tibble::tibble(
  iteration = 1:study_design$num_runs,
  cam_captures = NA,
  count_data = NA,
  encounter_data = NA,
  stay_time_all = NA,
  stay_time_raw = NA,
  stay_time_data = NA,
  TTE_data_all = NA,
  TTE_data_raw = NA,
  TTE_data = NA,
  STE_data = NA
)

# # Calculate total area using available space
# study_design <- study_design %>%
#   dplyr::mutate(
#     tot_A = unique(
#       study_design %>%
#         dplyr::reframe(
#           bounds = unlist(bounds),
#           tot_A = (bounds[2] - bounds[1])^2 * sum(!is.na(lscape_defs$Road)) / q
#         ) %>%
#         pull(tot_A)
#     )
#   )

# # camera data summaries
cam_design <- cam_design %>%
  dplyr::mutate(percent_cam_coverage = ncam * cam_A / study_design$tot_A)

dd_og <- c()
dd_test <- c()
dd_multi_N <- c()
dd <- c()

# Multi-run simulations
for (run in 1:study_design$num_runs) {
  print(paste("Run", run, "of", study_design$num_runs))
  
  # Create custom landscape (tag = "grid", "circ", "squares", "metapop")
  lscape_defs <- lscape_creator(study_design, lscape_design)
  
  # Create covariate matrix with 0, 1 values
  study_design <- study_design %>% 
    dplyr::mutate(
      num_covariates = length(unlist(covariate_labels)),
      Z = list(create_covariate_mat(
        lscape_defs,
        study_design,
        unlist(covariate_labels)))
    )
  
  # Run agent-based model
  animalxy.all <- ABM_sim(study_design,
                          lscape_defs)
  
  stay_time_tele <- Collect_tele_data(animalxy.all, study_design)
  
  # If running new ABM simulation on each run, save
  # save_animal_data$data[run] <- list(animalxy.all)  
  
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


  stay_time_summary <- stay_time_tele %>%
    dplyr::group_by(speed) %>%
    dplyr::summarise(
      # cell_mu = mean(t_stay),
      # cell_sd = sd(t_stay),
      cam_mu = mean(t_stay * cam_design$cam_A / (study_design$dx * study_design$dx)),
      cam_sd = sd(t_stay * cam_design$cam_A / (study_design$dx * study_design$dx)),
      .groups = 'drop'
    ) %>%
    dplyr::mutate(stay_prop = cam_mu / sum(cam_mu)) %>%
    dplyr::rename(Speed = speed) %>%
    dplyr::arrange(desc(Speed))

  kappa.prior.mu.tdst <- log(stay_time_summary$cam_mu)
  kappa.prior.var.tdst <- log(stay_time_summary$cam_sd)

  habitat_summary <- lscape_defs %>%
    dplyr::group_by(Speed) %>%
    dplyr::summarise(
      n_lscape = dplyr::n()
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      prop_lscape = n_lscape / sum(n_lscape)
    ) %>%
    dplyr::left_join(
      cam_locs %>%
        dplyr::group_by(Speed) %>%
        dplyr::summarise(
          # mean_count = mean(count, na.rm = T),
          ncams = dplyr::n(),
          .groups = 'drop'
        ),
      by = dplyr::join_by(Speed)
    ) %>%
    dplyr::left_join(
      stay_time_summary,
      by = dplyr::join_by(Speed)
    )  %>%
    dplyr::mutate(
      prop_cams = ncams / sum(ncams, na.rm = T),
      # stay_prop_full = stay_prop / sum(stay_prop),
      d_coeff = n_lscape * prop_cams / stay_prop / (cam_design$cam_A * study_design$t_steps),
      d_test = prop_cams / prop_lscape * stay_prop
    ) %>%
    dplyr::arrange(desc(Speed)) %>%
    replace(is.na(.), 0)
  #
  cc <- get_count_data(cam_locs, all_data$cam_captures[[run]], animalxy.all %>%
                         dplyr::filter(t != 0))
  # 
  n_adj <- habitat_summary %>%
    dplyr::left_join(
      cc %>%
        dplyr::group_by(Speed) %>%
        dplyr::summarise(
          mean_count = mean(count, na.rm = T),
          .groups = 'drop'
        ),
      by = dplyr::join_by(Speed)
    ) %>%
    dplyr::mutate(
      mean_count = tidyr::replace_na(mean_count, 0),
      ncams = tidyr::replace_na(ncams, 0)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      # stay_prop_adj = 1 / stay_prop * ncams / sum(ncams, na.rm = T),
      n_hab = mean_count * n_lscape / (cam_design$cam_A * study_design$t_steps),
      n_habitat = mean_count * d_coeff,
      n_full = n_hab / stay_prop,
      # big_D = mean_count * sum(stay_prop) / stay_prop
    ) #%>%
    # dplyr::select(Speed, prop_lscape, mean_count, ncams, cam_mu, stay_prop, prop_cams, d_coeff, n_lscape, n_habitat)

  dd_og[run] <- mean(cc$count)  * study_design$tot_A / (cam_design$cam_A * study_design$t_steps)
  
  dd[run] <- sum(n_adj$n_habitat)

  dd_multi_N[run] <- mean(n_adj$n_full) 
  
  dd_test[run] <- dd[run] / sum(habitat_summary$d_test)
  
}

  mean(dd)
  
  median(dd)
  
  boxplot(dd)
  lines(c(0,3), c(100, 100))

  # plot_ABM(study_design,
  #          cam_design,
  #          cam_locs,
  #          animalxy.all)

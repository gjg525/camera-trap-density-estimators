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
  n_iter = 20000,
  burn_in = 10000,
  # staying time censors
  t_censor = 10, 
  run_models = list(1:5),
  covariate_labels = list(c("Slow", "Medium", "Fast"))
)

# Landscape design
lscape_design <- tibble::tibble(
  lscape_tag = "Random",
  default_kappa = 0, # Turning angle
  num_roads = 0,
  Speed_ID = list(c("Slow", "Medium", "Fast")),
  # Speed_mins = list(c(.05, .3, .9)),
  # Speed_maxes = list(c(.15, .5, 1.1)),
  Speed_mins = list(c(.19, .4, .9)),
  Speed_maxes = list(c(.21, .5, 1.1)),
  Trail_ID = list(c("On Trail", "Off Trail")),
  Trail_speed = list(c("Medium", "Medium"))
)

# Cam designs
cam_design <- tibble::tibble(
  ncam = 200,
  Design = "Random",
  Props = list(c(1, 1, 1)), # proportion of cameras placed on and off roads

  # Design = "Bias",
  # Props = list(c(0, 0, 1)), # proportion of cameras placed on and off roads
  # cam_length = study_design$dx * 0.3, # length of all viewshed sides
  cam_length = study_design$dx * 0.1, # length of all viewshed sides
  cam_A = cam_length ^ 2 / 2,
  tot_snaps = ncam * study_design$t_steps
)
# Design = "Bias",
# Props = list(c(0.8, 0.1, 0.1)), # proportion of cameras placed on and off roads  # Design = "Bias",
  # Props = list(c(1, 0, 0)), # proportion of cameras placed on and off roads
  # Design = "Bias",
  # Props = list(c(0.1, 0.8, 0.1)), # proportion of cameras placed on and off roads
  # Design = "Bias",
  # Props = list(c(0.1, 0.1, 0.8)), # proportion of cameras placed on and off roads
  # Design = "Bias",
  # Props = list(c(0, 1, 0)), # proportion of cameras placed on and off roads

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
  # stay_time_raw = NA,
  stay_time_data = NA
  # TTE_data_all = NA,
  # TTE_data_raw = NA,
  # TTE_data = NA,
  # STE_data = NA
)

# # Calculate total area using available space
# study_design <- study_design %>%
#   dplyr::mutate(
#     tot_A = unique(
#       study_design %>%
#         dplyr::reframe(
#           bounds = unlist(bounds),
#           tot_A = (bounds[2] - bounds[1])^2 * sum(!is.na(lscape_defs$Road)) / q
#           ) %>%
#         pull(tot_A)
#     )
#   )

# # camera data summaries
cam_design <- cam_design %>%
  dplyr::mutate(percent_cam_coverage = ncam * cam_A / study_design$tot_A)

# Multi-run simulations
for (run in 1:study_design$num_runs) {
  print(paste("Run", run, "of", study_design$num_runs))

  # # Create custom landscape (tag = "grid", "circ", "squares", "metapop")
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
  
  # Set first column as reference category for intercept
  study_design$Z[[1]][, 1] <- 1
  
  # Run agent-based model
  animalxy.all <- ABM_sim(study_design,
                          lscape_defs)
  
  stay_time_tele <- Collect_tele_data(animalxy.all, study_design)
  
  stay_time_summary <- stay_time_tele %>% 
    dplyr::group_by(speed) %>% 
    dplyr::summarise(
      cell_mu = mean(t_stay),
      cell_sd = sd(t_stay),
      # cam_mu = mean(t_stay * cam_design$cam_A / (study_design$dx * study_design$dx)),
      # cam_sd = sd(t_stay * cam_design$cam_A / (study_design$dx * study_design$dx)),
      .groups = 'drop'
    ) %>% 
    dplyr::mutate(
      sum_stay = sum(cell_mu),
      stay_prop = cell_mu / sum_stay,
      cell_sd = cell_sd / sum_stay
    ) %>% 
    dplyr::rename(Speed = speed) %>% 
    dplyr::arrange(desc(Speed))
  
  kappa.prior.mu.tdst <- log(stay_time_summary$stay_prop)
  kappa.prior.var.tdst <- stay_time_summary$cell_sd

  # If running new ABM simulation on each run, save
  save_animal_data$data[run] <- list(animalxy.all)
  
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
  # cam_captures <- get_cam_captures(animalxy)
  
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
      d_coeff = n_lscape * prop_cams / stay_prop / (cam_design$cam_A * study_design$t_steps)
    ) %>%
    dplyr::arrange(desc(Speed)) %>%
    replace(is.na(.), 0)
  # 
  # cc <- get_count_data(cam_locs, all_data$cam_captures[[run]], animalxy.all %>%
  #                        dplyr::filter(t != 0))
  # 
  # n_adj <- habitat_summary %>%
  #   dplyr::left_join(
  #     cc %>%
  #       dplyr::group_by(Speed) %>%
  #       dplyr::summarise(
  #         mean_count = mean(count, na.rm = T),
  #         .groups = 'drop'
  #       ),
  #     by = dplyr::join_by(Speed)
  #   ) %>%
  #   dplyr::mutate(
  #     mean_count = tidyr::replace_na(mean_count, 0),
  #     ncams = tidyr::replace_na(ncams, 0)
  #   ) %>%
  #   dplyr::ungroup() %>%
  #   dplyr::mutate(
  #     # stay_prop_adj = 1 / stay_prop * ncams / sum(ncams, na.rm = T),
  #     n_lscape = mean_count * n_lscape / (cam_design$cam_A * study_design$t_steps),
  #     n_habitat = mean_count * d_coeff
  #     # n_full = n_lscape / stay_prop
  #     # big_D = mean_count * sum(stay_prop) / stay_prop
  #   ) %>%
  #   dplyr::select(Speed, prop_lscape, mean_count, ncams, cell_mu, stay_prop, prop_cams, d_coeff, n_lscape, n_habitat)

  # # The sum over the adjustments should equal total abundance when cameras are placed in each lscape type
  # # Does not work if cameras aren't placed in every landscape type
  # sum(n_adj$n_lscape, na.rm = T)
  #
  # # This method is more reliable when one or more lscape type is missing
  # sum(n_adj$n_habitat)

  # Run models only if any data points were collected
  if (nrow(all_data$cam_captures[[run]]) > 0) {
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
    encounter_data <- all_data$encounter_data[[run]]
    stay_time_data <- all_data$stay_time_data[[run]]
    
    # n_cores <- parallel::detectCores() $ Check number of cores available
    my_cluster <- parallel::makeCluster(3, type = "PSOCK")
    # register cluster to be used by %dopar%
    doParallel::registerDoParallel(cl = my_cluster)
    
    D.chain <- foreach::foreach(iter = unlist(study_design$run_models)) %dopar% {
      ###################################
      # TDST no priors, camera stay time data
      ###################################
      if (iter == 1) {
        if (sum(count_data$count) == 0) {
          D.TDST.MCMC <- NA
          SD.TDST.MCMC <- NA
        } else {
          chain.TDST <- fit.model.mcmc.TDST.cov(
            study_design = study_design,
            cam_design = cam_design,
            cam_locs = cam_locs,
            gamma_start = log(mean(count_data$count)),
            kappa_start = rep(log(mean(as.matrix(stay_time_data), na.rm = T)), study_design$num_covariates),
            gamma_prior_var = 10^4,
            kappa_prior_var = 10^4,
            gamma_tune = -1,
            kappa_tune = rep(-1, study_design$num_covariates),
            count_data_in = count_data$count,
            stay_time_data_in = as.matrix(stay_time_data)
          )
          
          # ## Posterior summaries
          # pop.ind.TDST <- which(names(chain.TDST) == "u")
          # MCMC.parms.TDST.cov <- coda::as.mcmc(do.call(cbind, chain.TDST[-pop.ind.TDST])[-c(1:study_design$burn_in), ])
          # summary(MCMC.parms.TDST.cov)
          
          # plot(chain.TDST$tot_u[study_design$burn_in:study_design$n_iter])
          D.TDST.MCMC <- mean(chain.TDST$tot_u[study_design$burn_in:study_design$n_iter])
          SD.TDST.MCMC <- sd(chain.TDST$tot_u[study_design$burn_in:study_design$n_iter])
          
          if (any(colMeans(chain.TDST$accept[study_design$burn_in:study_design$n_iter, ]) < 0.2) ||
              any(colMeans(chain.TDST$accept[study_design$burn_in:study_design$n_iter, ]) > 0.7)) {
            warning(("TDST accept rate OOB"))
            D.TDST.MCMC <- NA
            SD.TDST.MCMC <- NA
          }
        }
        D.chain <- tibble::tibble(
          iteration = run,
          Model = "TDST",
          Covariate = "Covariate",
          Est = D.TDST.MCMC,
          SD = SD.TDST.MCMC
          # all_results = list(chain.TDST)
        )
      }
      
      ###################################
      # TDST w/priors, tele stay time
      ###################################
      if (iter == 2) {
        if (sum(count_data$count) == 0) {
          D.TDST.MCMC <- NA
          SD.TDST.MCMC <- NA
        } else {
          chain.TDST.priors <- fit.model.mcmc.TDST.cov(
            study_design = study_design,
            cam_design = cam_design,
            cam_locs = cam_locs,
            gamma_start = log(mean(count_data$count)),
            kappa_start = kappa.prior.mu.tdst,
            gamma_prior_var = 10^4,
            kappa_prior_mu = kappa.prior.mu.tdst,
            kappa_prior_var = kappa.prior.var.tdst,
            gamma_tune = -1,
            kappa_tune = rep(-1, study_design$num_covariates),
            count_data_in = count_data$count,
            stay_time_data_in = NULL
          )
          
          # ## Posterior summaries
          # pop.ind.TDST <- which(names(chain.TDST.priors) == "u")
          # MCMC.parms.TDST.cov <- coda::as.mcmc(do.call(cbind, chain.TDST.priors[-pop.ind.TDST])[-c(1:study_design$burn_in), ])
          # summary(MCMC.parms.TDST.cov)
          
          # plot(chain.TDST.priors$tot_u[study_design$burn_in:study_design$n_iter])
          D.TDST.MCMC <- mean(chain.TDST.priors$tot_u[study_design$burn_in:study_design$n_iter])
          SD.TDST.MCMC <- sd(chain.TDST.priors$tot_u[study_design$burn_in:study_design$n_iter])
          
          if (any(colMeans(chain.TDST.priors$accept[study_design$burn_in:study_design$n_iter, ]) < 0.2) ||
              any(colMeans(chain.TDST.priors$accept[study_design$burn_in:study_design$n_iter, ]) > 0.7)) {
            warning(("TDST accept rate OOB"))
            D.TDST.MCMC <- NA
            SD.TDST.MCMC <- NA
          }
        }
        D.chain <- tibble::tibble(
          iteration = run,
          Model = "TDST w/ Priors",
          Covariate = "Covariate",
          Est = D.TDST.MCMC,
          SD = SD.TDST.MCMC
          # all_results = list(chain.TDST.priors)
        )
      }
      
      ########################################
      ## PR no covariates
      ########################################
      if (iter == 3) {
        if (sum(count_data$count) == 0) {
          D.PR.MCMC <- NA
          SD.PR.MCMC <- NA
        } else {
          chain.PR <- fit.model.mcmc.PR(
            study_design = study_design,
            cam_design = cam_design,
            gamma_start = log(mean(count_data$count)),
            gamma_prior_var = 10^4,
            gamma_tune = -1,
            count_data_in = count_data$count
          )
          
          # ## Posterior summaries
          # MCMC.parms.PR <- coda::as.mcmc(do.call(cbind, chain.PR)[-c(1:study_design$burn_in), ])
          # summary(MCMC.parms.PR)
          
          # plot(chain.PR$tot_u[study_design$burn_in:study_design$n_iter])
          D.PR.MCMC <- mean(chain.PR$tot_u[study_design$burn_in:study_design$n_iter])
          SD.PR.MCMC <- sd(chain.PR$tot_u[study_design$burn_in:study_design$n_iter])
          
          if (mean(chain.PR$accept[study_design$burn_in:study_design$n_iter, ]) < 0.2 || mean(chain.PR$accept[study_design$burn_in:study_design$n_iter, ]) > 0.7) {
            warning(("Mean Count accept rate OOB"))
            D.PR.MCMC <- NA
            SD.PR.MCMC <- NA
          }
        }
        D.chain <- tibble::tibble(
          iteration = run,
          Model = "PR",
          Covariate = "Non-Covariate",
          Est = D.PR.MCMC,
          SD = SD.PR.MCMC
          # all_results = list(chain.PR)
        )
      }
      
      ########################################
      ## PR w/ covariates
      ########################################
      if (iter == 4) {
        if (sum(count_data$count) == 0) {
          D.PR.MCMC.cov <- NA
          SD.PR.MCMC.cov <- NA
        } else {
          chain.PR.cov <- fit.model.mcmc.PR.cov(
            study_design = study_design,
            cam_design = cam_design,
            cam_locs = cam_locs,
            gamma_start = rep(log(mean(count_data$count)), study_design$num_covariates),
            gamma_prior_var = 10^4,
            gamma_tune = rep(-1, study_design$num_covariates),
            count_data_in = count_data$count
          )
          
          # ## Posterior summaries
          # pop.ind.PR <- which(names(chain.PR.cov) == "u")
          # MCMC.parms.PR.cov <- coda::as.mcmc(do.call(cbind, chain.PR.cov[-pop.ind.PR])[-c(1:study_design$burn_in), ])
          # summary(MCMC.parms.PR.cov)
          
          # plot(chain.PR.cov$tot_u[study_design$burn_in:study_design$n_iter])
          D.PR.MCMC.cov <- mean(chain.PR.cov$tot_u[study_design$burn_in:study_design$n_iter])
          SD.PR.MCMC.cov <- sd(chain.PR.cov$tot_u[study_design$burn_in:study_design$n_iter])
          
          if (mean(chain.PR.cov$accept[study_design$burn_in:study_design$n_iter, ]) < 0.2 || mean(chain.PR.cov$accept[study_design$burn_in:study_design$n_iter, ]) > 0.7) {
            warning(("Mean Count accept rate OOB"))
            D.PR.MCMC.cov <- NA
            SD.PR.MCMC.cov <- NA
          }
        }
        D.chain <- tibble::tibble(
          iteration = run,
          Model = "PR",
          Covariate = "Covariate",
          Est = D.PR.MCMC.cov,
          SD = SD.PR.MCMC.cov
          # all_results = list(chain.PR.cov)
        )
      }
      
      ########################################
      ## PR habitat model
      ########################################
      if (iter == 5) {
        if (sum(count_data$count) == 0) {
          D.PR.MCMC.cov <- NA
          SD.PR.MCMC.cov <- NA
        } else {
          chain.PR.habitat <- fit.model.mcmc.PR.habitat(
            study_design = study_design,
            cam_design = cam_design,
            cam_locs = cam_locs,
            gamma_start = rep(log(mean(count_data$count)), study_design$num_covariates),
            gamma_prior_var = 10^4,
            gamma_tune = rep(-1, study_design$num_covariates),
            kappa_start = kappa.prior.mu.tdst,
            kappa_prior_mu = kappa.prior.mu.tdst,
            kappa_prior_var = kappa.prior.var.tdst,
            kappa_tune = rep(-1, study_design$num_covariates),
            count_data_in = count_data,
            habitat_summary
          )
          
          # ## Posterior summaries
          # pop.ind.PR.habitat <- which(names(chain.PR.habitat) == "u")
          # MCMC.parms.PR.habitat <- coda::as.mcmc(do.call(cbind, chain.PR.habitat[-pop.ind.PR])[-c(1:study_design$burn_in), ])
          # summary(MCMC.parms.PR.habitat)
          
          # plot(chain.PR.habitat$tot_u[study_design$burn_in:study_design$n_iter])
          D.PR.MCMC.habitat <- mean(chain.PR.habitat$tot_u[study_design$burn_in:study_design$n_iter])
          SD.PR.MCMC.habitat <- sd(chain.PR.habitat$tot_u[study_design$burn_in:study_design$n_iter])
          
          if (any(colMeans(chain.PR.habitat$accept[study_design$burn_in:study_design$n_iter, ]) < 0.2) || any(colMeans(chain.PR.habitat$accept[study_design$burn_in:study_design$n_iter, ]) > 0.7)) {
            warning(("Mean Count accept rate OOB"))
            D.PR.MCMC.habitat <- NA
            SD.PR.MCMC.habitat <- NA
          }
        }
        D.chain <- tibble::tibble(
          iteration = run,
          Model = "PR Habitat",
          Covariate = "Non-Covariate",
          Est = D.PR.MCMC.habitat,
          SD = SD.PR.MCMC.habitat
          # all_results = list(chain.PR.habitat)
        )
      }
      
      D.chain <- D.chain
      
      # D.all[D.all$iteration == run, ] <- D.mcmc
    }
    stopCluster(my_cluster)
  }
  
  D.all[(num_models * (run - 1) + 1):(num_models * run), ] <- dplyr::bind_rows(D.chain)
}

# # # Remove outlier estimates
# D.all$Est[D.all$Est > 5 * study_design$tot_animals] <- NA
# D.all$SD[D.all$SD > 5 * study_design$tot_animals] <- NA

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
plot_multirun_mape(D.all %>%
                      dplyr::filter(is.finite(Est)),
                   study_design$tot_animals)
plot_multirun_CV(D.all %>% 
                   dplyr::filter(is.finite(Est)))
# plot_multirun_hist(D.all)

# plot_ABM(study_design,
#          cam_design,
#          cam_locs,
#          animalxy.all)
# plot_space_use(study_design,
#                animalxy.all)
#
# plot_count_data(count_data_in = all_data$count_data[[1]],
#                 fill = "Speed")
# plot_encounter_data(encounter_data_in = all_data$encounter_data[[1]],
#                 fill = "Speed")
# plot_staytime_data(stay_time_raw_in = all_data$stay_time_raw[[1]],
#                    fill = "Speed")

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


# results_fast_all_cam <- list(
#   save_animal_data,
#   study_design,
#   cam_design,
#   lscape_design,
#   all_data,
#   D.all
# )
# 
# save(results_fast_all_cam, file = "Sim_results/results_fast_all_cam.RData")
# 

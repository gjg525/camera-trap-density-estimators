################################################################################
Run_mcmc_models <- function(study_design,
                            all_data,
                            cam_design,
                            cam_locs,
                            run) {

  count_data <- all_data$count_data[[1]]
  encounter_data <- all_data$encounter_data[[1]]
  stay_time_data <- all_data$stay_time_data[[1]]

  ################################
  # MCMC methods
  ################################
  # Run models in parallel
  run_iter <- get_mcmc_runs(study_design)
  # get_mcmc_runs <- function(study_design) {
  #   iter_legend <- tibble::tibble(
  #     Model = c("TDST", "REST", "REST", "TTE", "TTE", "PR", "PR", "STE"),
  #     Covariate = c("Covariate", "Non-Covariate", "Covariate", "Non-Covariate",
  #                   "Covariate", "Non-Covariate", "Covariate","Non-Covariate"),
  #     iter = 1:8
  #   )
  #
  #   run_iter <- study_design |>
  #     tidyr::expand(Model = unlist(run_models),
  #            Covariate = unlist(run_covariates)) |>
  #     dplyr::left_join(iter_legend, by = c("Model", "Covariate")) |>
  #     dplyr::summarise(iter = iter[!is.na(iter)]) |>
  #     dplyr::arrange(iter) |>
  #     dplyr::pull(iter)
  # }

  # n_cores <- parallel::detectCores() $ Check number of cores available
  my_cluster <- parallel::makeCluster(3, type = "PSOCK")
  # register cluster to be used by %dopar%
  doParallel::registerDoParallel(cl = my_cluster)

  D.chain <- foreach::foreach(iter = run_iter) %dopar% {
    ###################################
    # TDST
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
          gamma_prior_var = 10^6,
          kappa_prior_var = 10^6,
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
        SD = SD.TDST.MCMC,
        all_results = list(chain.TDST)
      )
    }

    ###################################
    # REST no covariates
    ###################################
    if (iter == 2) {
      if (sum(encounter_data$encounter) == 0) {
        D.REST.MCMC <- NA
        SD.REST.MCMC <- NA
      } else {
        chain.REST <- fit.model.mcmc.REST(
          study_design = study_design,
          cam_design = cam_design,
          gamma_start = log(mean(encounter_data$encounter)),
          kappa_start = log(mean(as.matrix(stay_time_data), na.rm = T)),
          gamma_prior_var = 10^6,
          kappa_prior_var = 10^6,
          gamma_tune = -1,
          kappa_tune = -1,
          encounter_data_in = encounter_data$encounter,
          stay_time_data_in = as.matrix(stay_time_data)
        )

        # ## Posterior summaries
        # MCMC.parms.REST <- coda::as.mcmc(do.call(cbind, chain.REST)[-c(1:study_design$burn_in), ])
        # summary(MCMC.parms.REST)

        # plot(chain.REST$tot_u[study_design$burn_in:study_design$n_iter])
        D.REST.MCMC <- mean(chain.REST$tot_u[study_design$burn_in:study_design$n_iter])
        SD.REST.MCMC <- sd(chain.REST$tot_u[study_design$burn_in:study_design$n_iter])

        if (any(colMeans(chain.REST$accept[study_design$burn_in:study_design$n_iter, ]) < 0.2) || any(colMeans(chain.REST$accept[study_design$burn_in:study_design$n_iter, ]) > 0.7)) {
          warning(("REST accept rate OOB"))
          D.REST.MCMC <- NA
          SD.REST.MCMC <- NA
        }
      }

      D.chain <- tibble::tibble(
        iteration = run,
        Model = "REST",
        Covariate = "Non-Covariate",
        Est = D.REST.MCMC,
        SD = SD.REST.MCMC,
        all_results = list(chain.REST)
      )
    }

    ###################################
    # REST w/ covariates
    ###################################
    if (iter == 3) {
      if (sum(encounter_data$encounter) == 0) {
        D.REST.MCMC.cov <- NA
        SD.REST.MCMC.cov <- NA
      } else {
        chain.REST.cov <- fit.model.mcmc.REST.cov(
          study_design = study_design,
          cam_design = cam_design,
          cam_locs = cam_locs,
          gamma_start = rep(log(mean(encounter_data$encounter)), study_design$num_covariates),
          kappa_start = rep(log(mean(as.matrix(stay_time_data), na.rm = T)), study_design$num_covariates),
          gamma_prior_var = 10^6,
          kappa_prior_var = 10^6,
          gamma_tune = rep(-1, study_design$num_covariates),
          kappa_tune = rep(-1, study_design$num_covariates),
          encounter_data_in = encounter_data$encounter,
          stay_time_data_in = as.matrix(stay_time_data)
        )

        # ## Posterior summaries
        # pop.ind.REST <- which(names(chain.REST.cov) == "u")
        # MCMC.parms.REST.cov <- coda::as.mcmc(do.call(cbind, chain.REST.cov[-pop.ind.REST])[-c(1:study_design$burn_in), ])
        # summary(MCMC.parms.REST.cov)

        # plot(chain.REST$tot_u[study_design$burn_in:study_design$n_iter])
        D.REST.MCMC.cov <- mean(chain.REST.cov$tot_u[study_design$burn_in:study_design$n_iter])
        SD.REST.MCMC.cov <- sd(chain.REST.cov$tot_u[study_design$burn_in:study_design$n_iter])

        if (any(colMeans(chain.REST.cov$accept[study_design$burn_in:study_design$n_iter, ]) < 0.2) ||
          any(colMeans(chain.REST.cov$accept[study_design$burn_in:study_design$n_iter, ]) > 0.7)) {
          warning(("REST accept rate OOB"))
          D.REST.MCMC.cov <- NA
          SD.REST.MCMC.cov <- NA
        }
      }
      D.chain <- tibble::tibble(
        iteration = run,
        Model = "REST",
        Covariate = "Covariate",
        Est = D.REST.MCMC.cov,
        SD = SD.REST.MCMC.cov,
        all_results = list(chain.REST.cov)
      )
    }

    ########################################
    ## TTE no covariates
    ########################################
    if (iter == 4) {
      if (all(is.na(TTE_data))) {
        D.TTE.MCMC <- NA
        SD.TTE.MCMC <- NA
      } else {
        chain.TTE <- fit.model.mcmc.TTE(
          study_design = study_design,
          cam_design = cam_design,
          gamma_start = log(mean(TTE_data, na.rm = T)),
          kappa_start = log(mean(as.matrix(stay_time_data), na.rm = T)),
          gamma_prior_var = 10^6,
          kappa_prior_var = 10^6,
          gamma_tune = -1,
          kappa_tune = -1,
          TTE_data_in = TTE_data,
          stay_time_data_in = as.matrix(stay_time_data)
        )

        # ## Posterior summaries
        # MCMC.parms.TTE <- coda::as.mcmc(do.call(cbind, chain.TTE)[-c(1:study_design$burn_in), ])
        # summary(MCMC.parms.TTE)

        # plot(chain.TTE$tot_u[study_design$burn_in:study_design$n_iter])
        D.TTE.MCMC <- mean(chain.TTE$tot_u[study_design$burn_in:study_design$n_iter])
        SD.TTE.MCMC <- sd(chain.TTE$tot_u[study_design$burn_in:study_design$n_iter])

        if (any(colMeans(chain.TTE$accept[study_design$burn_in:study_design$n_iter, ]) < 0.2) || any(colMeans(chain.TTE$accept[study_design$burn_in:study_design$n_iter, ]) > 0.7)) {
          warning(("TTE accept rate OOB"))
          D.TTE.MCMC <- NA
          SD.TTE.MCMC <- NA
        }
      }
      D.chain <- tibble::tibble(
        iteration = run,
        Model = "TTE",
        Covariate = "Non-Covariate",
        Est = D.TTE.MCMC,
        SD = SD.TTE.MCMC,
        all_results = list(chain.TTE)
      )
    }

    ########################################
    ## TTE w/ covariates
    ########################################
    if (iter == 5) {
      if (all(is.na(TTE_data))) {
        D.TTE.MCMC.cov <- NA
        SD.TTE.MCMC.cov <- NA
      } else {
        chain.TTE.cov <- fit.model.mcmc.TTE.cov(
          study_design = study_design,
          cam_design = cam_design,
          cam_locs = cam_locs,
          gamma_start = rep(log(mean(TTE_data, na.rm = T)), study_design$num_covariates),
          kappa_start = rep(log(mean(as.matrix(stay_time_data), na.rm = T)), study_design$num_covariates),
          gamma_prior_var = 10^6,
          kappa_prior_var = 10^6,
          gamma_tune = rep(-1, study_design$num_covariates),
          kappa_tune = rep(-1, study_design$num_covariates),
          TTE_data_in = TTE_data,
          stay_time_data_in = as.matrix(stay_time_data)
        )

        # ## Posterior summaries
        # pop.ind.TTE <- which(names(chain.TTE.cov) == "u")
        # MCMC.parms.TTE.cov <- coda::as.mcmc(do.call(cbind, chain.TTE.cov[-pop.ind.TTE])[-c(1:study_design$burn_in), ])
        # summary(MCMC.parms.TTE.cov)

        # plot(chain.TTE.cov$tot_u[study_design$burn_in:study_design$n_iter])
        D.TTE.MCMC.cov <- mean(chain.TTE.cov$tot_u[study_design$burn_in:study_design$n_iter])
        SD.TTE.MCMC.cov <- sd(chain.TTE.cov$tot_u[study_design$burn_in:study_design$n_iter])

        if (any(colMeans(chain.TTE.cov$accept[study_design$burn_in:study_design$n_iter, ]) < 0.2) ||
          any(colMeans(chain.TTE.cov$accept[study_design$burn_in:study_design$n_iter, ]) > 0.7)) {
          warning(("TTE accept rate OOB"))
          D.TTE.MCMC.cov <- NA
          SD.TTE.MCMC.cov <- NA
        }
      }
      D.chain <- tibble::tibble(
        iteration = run,
        Model = "TTE",
        Covariate = "Covariate",
        Est = D.TTE.MCMC.cov,
        SD = SD.TTE.MCMC.cov,
        all_results = list(chain.TTE.cov)
      )
    }

    ########################################
    ## PR no covariates
    ########################################
    if (iter == 6) {
      if (sum(count_data$count) == 0) {
        D.PR.MCMC <- NA
        SD.PR.MCMC <- NA
      } else {
        chain.PR <- fit.model.mcmc.PR(
          study_design = study_design,
          cam_design = cam_design,
          gamma_start = log(mean(count_data$count)),
          gamma_prior_var = 10^6,
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
        SD = SD.PR.MCMC,
        all_results = list(chain.PR)
      )
    }

    ########################################
    ## PR w/ covariates
    ########################################
    if (iter == 7) {
      if (sum(count_data$count) == 0) {
        D.PR.MCMC.cov <- NA
        SD.PR.MCMC.cov <- NA
      } else {
        chain.PR.cov <- fit.model.mcmc.PR.cov(
          study_design = study_design,
          cam_design = cam_design,
          cam_locs = cam_locs,
          gamma_start = rep(log(mean(count_data$count)), study_design$num_covariates),
          gamma_prior_var = 10^6,
          gamma_tune = rep(-1, study_design$num_covariates),
          count_data_in = count_data$count
        )

        # ## Posterior summaries
        # pop.ind.PR <- which(names(chain.PR.cov) == "u")
        # MCMC.parms.PR.cov <- coda::as.mcmc(do.call(cbind, chain.PR.cov[-pop.ind.PR])[-c(1:study_design$burn_in), ])
        # summary(MCMC.parms.PR.cov)

        # plot(chain.PR$tot_u[study_design$burn_in:study_design$n_iter])
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
        SD = SD.PR.MCMC.cov,
        all_results = list(chain.PR.cov)
      )
    }

    ########################################
    ## STE no covariates
    ########################################
    if (iter == 8) {
      if (all(is.na(STE_data$STE))) {
        D.STE.MCMC <- NA
        SD.STE.MCMC <- NA
      } else {
        chain.STE <- fit.model.mcmc.STE(
          study_design = study_design,
          cam_design = cam_design,
          gamma_start = mean(STE_data$STE, na.rm = T),
          gamma_prior_var = 10^6,
          gamma_tune = -1,
          STE_data_in = STE_data$STE
        )

        # ## Posterior summaries
        # # plot(chain.STE$tot_u[study_design$burn_in:study_design$n_iter])
        # MCMC.parms.STE <- coda::as.mcmc(do.call(cbind, chain.STE)[-c(1:study_design$burn_in), ])
        # summary(MCMC.parms.STE)

        # plot(chain.STE$tot_u[study_design$burn_in:study_design$n_iter])
        D.STE.MCMC <- mean(chain.STE$tot_u[study_design$burn_in:study_design$n_iter])
        SD.STE.MCMC <- sd(chain.STE$tot_u[study_design$burn_in:study_design$n_iter])

        if (mean(chain.STE$accept[study_design$burn_in:study_design$n_iter]) < 0.2 || mean(chain.STE$accept[study_design$burn_in:study_design$n_iter, ]) > 0.7) {
          warning(("STE accept rate OOB"))
          D.STE.MCMC <- NA
          SD.STE.MCMC <- NA
        }
      }

      D.chain <- tibble::tibble(
        iteration = run,
        Model = "STE",
        Covariate = "Non-Covariate",
        Est = D.STE.MCMC,
        SD = SD.STE.MCMC,
        all_results = list(chain.STE)
      )
    }
    D.chain <- D.chain
  }

  # Stop cluster
  stopCluster(my_cluster)

  D.all <- dplyr::bind_rows(D.chain)
  # D.all <- D.all |>
  #   dplyr::add_row(dplyr::bind_rows(D.chain))
}

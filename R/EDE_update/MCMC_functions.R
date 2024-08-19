# Metropolis-within-Gibbs sampler
########################################
# MCMC for TDST w/ covariates
########################################
################################################################################
#' @export
#'
fit.model.mcmc.TDST.cov <- function(study_design,
                                    cam_design,
                                    cam_locs,
                                    gamma_start,
                                    kappa_start,
                                    gamma_prior_var = 10^6,
                                    kappa_prior_var = 10^6,
                                    gamma_tune,
                                    kappa_tune,
                                    count_data_in,
                                    stay_time_data_in
                                    ) {
  # Initialize variables
  n_iter <- study_design$n_iter
  t_steps <- study_design$t_steps
  cam_A <- cam_design$cam_A
  cell_A <- study_design$dx * study_design$dy
  censor <- study_design$t_censor
  ncam <- cam_design$ncam
  covariate_labels <- unlist(study_design$covariate_labels)
  num_covariates <- study_design$num_covariates
  Z <- matrix(unlist(study_design$Z), study_design$q, study_design$num_covariates)


  # Variables that will be saved
  gamma <- matrix(NA, n_iter + 1, 1)
  kappa <- matrix(NA, n_iter + 1, num_covariates)
  tot_u <- matrix(NA, n_iter + 1, 1)
  accept <- matrix(NA, n_iter + 1, 1 + num_covariates)
  gamma[1] <- gamma_start
  colnames(gamma) <- "gamma"
  kappa[1, ] <- kappa_start
  colnames(kappa) <- paste0("kappa.", covariate_labels)
  colnames(tot_u) <- "Total estimate"
  colnames(accept) <- c(
    "accept.rate.gamma",
    paste0("accept.rate.kappa.", covariate_labels)
  )

  # Account for censored times
  stay_time_data_all <- stay_time_data_in
  stay_time_data_all[stay_time_data_all >= censor] <- NA
  stay_time_data_censor <- stay_time_data_in
  stay_time_data_censor[stay_time_data_censor < censor] <- NA

  # Initialize with landscape-scale covariates
  beta <- exp(gamma[1])
  phi <- exp_na_covs(Z, kappa[1, ])
  sum_phi <- sum(phi, na.rm = T)
  if (is.infinite(sum_phi)) {
    u <- phi * 0
  } else {
    u <- beta * phi / sum_phi
  }

  tune_check <- 100
  batch_n <- 0

  # Begin MCMC loop
  # prog_bar <- txtProgressBar(min = 0,max = n_iter,style = 3,width = 50,char = "=")
  for (i in 1:n_iter) {
    # setTxtProgressBar(prog_bar, i)
    # Sample gamma
    gamma_star <- rnorm(1, gamma[i], exp(2 * gamma_tune))
    beta_star <- exp(gamma_star)
    sum_phi <- sum(phi, na.rm = T)
    if (is.infinite(sum_phi)) {
      u_star <- phi * 0
    } else {
      u_star <- beta_star * phi / sum_phi
    }

    # repeat estimated parms for fitting
    u_star_cams <- u_star[cam_locs$lscape_index] * cam_A / cell_A
    u_cams <- u[cam_locs$lscape_index] * cam_A / cell_A

    if (all(u_star_cams > 0) & all(!is.infinite(u_star_cams))) {
      mh1 <- sum(dpois(count_data_in, u_star_cams, log = TRUE), na.rm = TRUE) +
        sum(dnorm(gamma_star, 0, gamma_prior_var^0.5, log = TRUE))
      mh2 <- sum(dpois(count_data_in, u_cams, log = TRUE), na.rm = TRUE) +
        sum(dnorm(gamma[i, ], 0, gamma_prior_var^0.5, log = TRUE))
      mh <- exp(mh1 - mh2)

      if (mh > runif(1) & !is.na(mh)) {
        gamma[i + 1] <- gamma_star
        accept[i + 1, 1] <- 1
        u <- u_star
        beta <- beta_star
      } else {
        gamma[i + 1] <- gamma[i]
        accept[i + 1, 1] <- 0
      }
    } else {
      gamma[i + 1] <- gamma[i]
      accept[i + 1, 1] <- 0
    }
    beta <- exp(gamma[i + 1])

    # Sample kappa
    kappa_temp <- kappa[i, ]
    for (kk in 1:num_covariates) {
      kappa_star <- kappa[i, ]
      kappa_star[kk] <- rnorm(1, kappa[i, kk], exp(2 * kappa_tune[kk]))
      phi_star <- exp_na_covs(Z, kappa_star)
      sum_phi_star <- sum(phi_star, na.rm = T)
      if (is.infinite(sum_phi_star)) {
        u_star <- phi_star * 0
      } else {
        u_star <- beta * phi_star / sum_phi_star
      }

      # # repeat estimated parms for fitting
      phi_star_cams <- phi_star[cam_locs$lscape_index] * cam_A / cell_A
      phi_star_cam_rep <- matrix(rep(phi_star_cams, dim(stay_time_data_all)[2]),
        nrow = ncam, ncol = dim(stay_time_data_all)[2]
      )
      phi_all_cams <- phi[cam_locs$lscape_index] * cam_A / cell_A
      phi_all_cam_rep <- matrix(rep(phi_all_cams, dim(stay_time_data_all)[2]),
        nrow = ncam, ncol = dim(stay_time_data_all)[2]
      )


      if (all(1 / phi_star_cams > 0) & all(!is.infinite(1 / phi_star_cams))) {
        mh1 <- sum(dexp(stay_time_data_all, 1 / phi_star_cam_rep, log = TRUE), na.rm = TRUE) +
          sum(pexp(stay_time_data_censor, 1 / phi_star_cam_rep, lower.tail = F, log = TRUE), na.rm = TRUE) +
          sum(dnorm(kappa_star, 0, kappa_prior_var^0.5, log = TRUE))
        mh2 <- sum(dexp(stay_time_data_all, 1 / phi_all_cam_rep, log = TRUE), na.rm = TRUE) +
          sum(pexp(stay_time_data_censor, 1 / phi_all_cam_rep, lower.tail = F, log = TRUE), na.rm = TRUE) +
          sum(dnorm(kappa[i, ], 0, kappa_prior_var^0.5, log = TRUE))
        mh <- exp(mh1 - mh2)

        # Alternate: fit staying time with gamma distribution
        # if(all(1/phi_star_cams>0) & all(!is.infinite(1/phi_star_cams))){
        #   mh1 <- sum(dgamma(stay_time_data_all, phi_star_cam_rep,log=TRUE),na.rm=TRUE) +
        #     sum(pgamma(stay_time_data_censor, phi_star_cam_rep, lower.tail = F, log = TRUE),na.rm=TRUE) +
        #     sum(dnorm(kappa_star,0,kappa_prior_var^0.5,log=TRUE))
        #   mh2 <- sum(dgamma(stay_time_data_all, phi_all_cam_rep,log=TRUE),na.rm=TRUE) +
        #     sum(pgamma(stay_time_data_censor, phi_all_cam_rep, lower.tail = F, log = TRUE),na.rm=TRUE) +
        #     sum(dnorm(kappa[i,],0,kappa_prior_var^0.5,log=TRUE))
        #   mh <- exp(mh1-mh2)

        if (mh > runif(1) & !is.na(mh)) {
          kappa[i, ] <- kappa_star
          accept[i + 1, 1 + kk] <- 1
          u <- u_star
          phi <- phi_star
        } else {
          kappa[i, ] <- kappa[i, ]
          accept[i + 1, 1 + kk] <- 0
        }
      } else {
        kappa[i, ] <- kappa[i, ]
        accept[i + 1, 1 + kk] <- 0
      }
    }
    kappa[i + 1, ] <- kappa[i, ]
    phi <- exp_na_covs(Z, kappa[i + 1, ])

    tot_u[i + 1] <- sum(u, na.rm = T) / t_steps

    # Update tuning parms
    if (i %% tune_check == 0) {
      batch_n <- batch_n + 1
      delta_n <- batch_n^-1
      # delta_n <- min(0.01,batch_n^-1)
      accept_gamma_check <- mean(accept[(i - tune_check + 1):i, 1], na.rm = T)
      gamma_tune[which(accept_gamma_check > 0.44)] <- gamma_tune[which(accept_gamma_check > 0.44)] + delta_n
      gamma_tune[which(accept_gamma_check <= 0.44)] <- gamma_tune[which(accept_gamma_check <= 0.44)] - delta_n
      accept_kappa_check <- colMeans(
        array(
          accept[
            (i - tune_check + 1):i,
            2:(1 + num_covariates)
          ],
          dim = c(tune_check, num_covariates)
        ),
        na.rm = T
      )
      kappa_tune[which(accept_kappa_check > 0.44)] <- kappa_tune[which(accept_kappa_check > 0.44)] + delta_n
      kappa_tune[which(accept_kappa_check <= 0.44)] <- kappa_tune[which(accept_kappa_check <= 0.44)] - delta_n
    }
  }
  # print("MCMC complete")

  list(accept = accept, gamma = gamma, kappa = kappa, tot_u = tot_u, u = u)
}


########################################
# MCMC for REST no covariates
########################################
################################################################################
#' @export
#'
fit.model.mcmc.REST <- function(study_design,
                                cam_design,
                                gamma_start,
                                kappa_start,
                                gamma_prior_var = 10^6,
                                kappa_prior_var = 10^6,
                                gamma_tune = -1,
                                kappa_tune = -1,
                                encounter_data_in,
                                stay_time_data_in
                                ) {
  n_iter <- study_design$n_iter
  t_steps <- study_design$t_steps
  cam_A <- cam_design$cam_A
  censor <- study_design$t_censor
  tot_A <- study_design$tot_A
  dt <- study_design$dt

  # Variables that will be saved
  gamma <- matrix(NA, n_iter + 1, 1)
  kappa <- matrix(NA, n_iter + 1, 1)
  tot_u <- matrix(NA, n_iter + 1, 1)
  accept <- matrix(NA, n_iter + 1, 2)

  gamma[1] <- gamma_start
  colnames(gamma) <- "gamma"
  kappa[1] <- kappa_start
  colnames(kappa) <- "kappa"
  colnames(tot_u) <- "Total estimate"
  colnames(accept) <- c("accept.rate.gamma", "accept.rate.kappa")

  # Account for censored times
  stay_time_data_all <- stay_time_data_in
  stay_time_data_all[stay_time_data_all > censor] <- NA
  stay_time_data_censor <- stay_time_data_in
  stay_time_data_censor[stay_time_data_censor <= censor] <- NA

  u <- exp(gamma[1])
  phi <- exp(kappa[1])
  eta <- u * dt * t_steps * cam_A / phi
  tot_u <- tot_A * u

  tune_check <- 100
  batch_n <- 0

  # Begin MCMC loop
  # prog_bar <- txtProgressBar(min = 0,max = n_iter,style = 3,width = 50,char = "=")
  for (i in 1:n_iter) {
    # setTxtProgressBar(prog_bar, i)
    # Sample gamma
    gamma_star <- rnorm(1, gamma[i], exp(2 * gamma_tune))
    u_star <- exp(gamma_star)
    eta_star <- u_star * dt * t_steps * cam_A / phi

    if (all(eta_star > 0)) {
      mh1 <- sum(dpois(encounter_data_in, eta_star, log = TRUE), na.rm = TRUE) +
        sum(dnorm(gamma_star, 0, gamma_prior_var^0.5, log = TRUE))
      mh2 <- sum(dpois(encounter_data_in, eta, log = TRUE), na.rm = TRUE) +
        sum(dnorm(gamma[i], 0, gamma_prior_var^0.5, log = TRUE))
      mh <- exp(mh1 - mh2)

      if (mh > runif(1) & !is.na(mh)) {
        gamma[i + 1] <- gamma_star
        accept[i + 1, 1] <- 1
        eta <- eta_star
      } else {
        gamma[i + 1] <- gamma[i]
        accept[i + 1, 1] <- 0
      }
    } else {
      gamma[i + 1] <- gamma[i]
      accept[i + 1, 1] <- 0
    }

    # Sample kappa
    kappa_star <- rnorm(1, kappa[i], exp(2 * kappa_tune))
    phi_star <- exp(kappa_star)

    if (1 / phi_star > 0 & !is.infinite(1 / phi_star)) {
      mh1 <- sum(dexp(stay_time_data_all, 1 / phi_star, log = TRUE), na.rm = TRUE) +
        sum(pexp(stay_time_data_censor, 1 / phi_star, lower.tail = F, log = TRUE), na.rm = TRUE) +
        sum(dnorm(kappa_star, 0, kappa_prior_var^0.5, log = TRUE))
      mh2 <- sum(dexp(stay_time_data_all, 1 / phi, log = TRUE), na.rm = TRUE) +
        sum(pexp(stay_time_data_censor, 1 / phi, lower.tail = F, log = TRUE), na.rm = TRUE) +
        sum(dnorm(kappa[i], 0, kappa_prior_var^0.5, log = TRUE))
      mh <- exp(mh1 - mh2)

      if (mh > runif(1) & !is.na(mh)) {
        kappa[i + 1] <- kappa_star
        accept[i + 1, 2] <- 1
        phi <- phi_star
      } else {
        kappa[i + 1] <- kappa[i]
        accept[i + 1, 2] <- 0
      }
    } else {
      kappa[i + 1] <- kappa[i]
      accept[i + 1, 2] <- 0
    }

    u <- exp(gamma[i + 1])
    phi <- exp(kappa[i + 1])
    eta <- u * dt * t_steps * cam_A / phi
    tot_u[i + 1] <- tot_A * u

    # Update tuning parms
    if (i %% tune_check == 0) {
      batch_n <- batch_n + 1
      delta_n <- batch_n^-1
      # delta_n <- min(0.01,batch_n^-1)
      accept_gamma_check <- mean(accept[(i - tune_check + 1):i, 1], na.rm = T)
      gamma_tune[accept_gamma_check > 0.44] <- gamma_tune[accept_gamma_check > 0.44] + delta_n
      gamma_tune[accept_gamma_check <= 0.44] <- gamma_tune[accept_gamma_check <= 0.44] - delta_n
      accept_kappa_check <- mean(accept[(i - tune_check + 1):i, 2], na.rm = T)
      kappa_tune[accept_kappa_check > 0.44] <- kappa_tune[accept_kappa_check > 0.44] + delta_n
      kappa_tune[accept_kappa_check <= 0.44] <- kappa_tune[accept_kappa_check <= 0.44] - delta_n
    }
  }
  # print("MCMC complete")

  list(accept = accept, gamma = gamma, kappa = kappa, tot_u = tot_u)
}


########################################
# MCMC for REST w/ covariates
########################################
################################################################################
#' @export
#'
fit.model.mcmc.REST.cov <- function(study_design,
                                    cam_design,
                                    cam_locs,
                                    gamma_start,
                                    kappa_start,
                                    gamma_prior_var = 10^6,
                                    kappa_prior_var = 10^6,
                                    gamma_tune,
                                    kappa_tune,
                                    encounter_data_in,
                                    stay_time_data_in
                                    ) {

  n_iter <- study_design$n_iter
  t_steps <- study_design$t_steps
  cam_A <- cam_design$cam_A
  cell_A <- study_design$dx * study_design$dy
  censor <- study_design$t_censor
  dt <- study_design$dt
  ncam <- cam_design$ncam
  covariate_labels <- unlist(study_design$covariate_labels)
  num_covariates <- study_design$num_covariates
  Z <- matrix(unlist(study_design$Z), study_design$q, study_design$num_covariates)

  # Variables that will be saved
  gamma <- matrix(NA, n_iter + 1, num_covariates)
  kappa <- matrix(NA, n_iter + 1, num_covariates)
  tot_u <- matrix(NA, n_iter + 1, 1)
  accept <- matrix(NA, n_iter + 1, 2 * num_covariates)

  gamma[1, ] <- gamma_start
  colnames(gamma) <- paste0("gamma.", covariate_labels)
  kappa[1, ] <- kappa_start
  colnames(kappa) <- paste0("kappa.", covariate_labels)
  colnames(tot_u) <- "Total estimate"
  colnames(accept) <- c(
    paste0("accept.rate.gamma.", covariate_labels),
    paste0("accept.rate.kappa.", covariate_labels)
  )

  # Account for censored times
  stay_time_data_all <- stay_time_data_in
  stay_time_data_all[stay_time_data_all > censor] <- NA
  stay_time_data_censor <- stay_time_data_in
  stay_time_data_censor[stay_time_data_censor <= censor] <- NA

  # Initialize with landscape-scale covariates
  u <- exp_na_covs(Z, gamma[1, ])
  phi <- exp_na_covs(Z, kappa[1, ])
  eta <- u * dt * t_steps / phi * cell_A
  tot_u[1] <- sum(u, na.rm = T) * cell_A

  tune_check <- 100
  batch_n <- 0

  # Begin MCMC loop
  # prog_bar <- txtProgressBar(min = 0,max = n_iter,style = 3,width = 50,char = "=")
  for (i in 1:n_iter) {
    # setTxtProgressBar(prog_bar, i)
    # Sample gamma
    for (gg in 1:num_covariates) {
      gamma_star <- gamma[i, ]
      gamma_star[gg] <- rnorm(1, gamma[i, gg], exp(2 * gamma_tune[gg]))
      u_star <- exp_na_covs(Z, gamma_star)
      eta_star <- u_star * dt * t_steps / phi * cell_A

      eta_star_cams <- eta_star[cam_locs$lscape_index]
      eta_all_cams <- eta[cam_locs$lscape_index]

      if (all(eta_star_cams > 0)) {
        mh1 <- sum(dpois(encounter_data_in, eta_star_cams, log = TRUE), na.rm = TRUE) +
          sum(dnorm(gamma_star, 0, gamma_prior_var^0.5, log = TRUE))
        mh2 <- sum(dpois(encounter_data_in, eta_all_cams, log = TRUE), na.rm = TRUE) +
          sum(dnorm(gamma[i, ], 0, gamma_prior_var^0.5, log = TRUE))
        mh <- exp(mh1 - mh2)

        if (mh > runif(1) & !is.na(mh)) {
          gamma[i, ] <- gamma_star
          accept[i + 1, gg] <- 1
          eta <- eta_star
        } else {
          gamma[i, ] <- gamma[i, ]
          accept[i + 1, gg] <- 0
        }
      } else {
        gamma[i, ] <- gamma[i, ]
        accept[i + 1, gg] <- 0
      }
    }
    gamma[i + 1, ] <- gamma[i, ]
    u <- exp_na_covs(Z, gamma[i + 1, ])
    eta <- u * dt * t_steps / phi * cell_A

    # Sample kappa
    for (kk in 1:num_covariates) {
      kappa_star <- kappa[i, ]
      kappa_star[kk] <- rnorm(1, kappa[i, kk], exp(2 * kappa_tune[kk]))
      phi_star <- exp_na_covs(Z, kappa_star)

      # repeat estimated parms for fitting
      phi_star_cams <- phi_star[cam_locs$lscape_index] * cam_A / cell_A
      phi_star_cam_rep <- matrix(rep(phi_star_cams, dim(stay_time_data_all)[2]),
        nrow = ncam, ncol = dim(stay_time_data_all)[2]
      )
      phi_all_cams <- phi[cam_locs$lscape_index] * cam_A / cell_A
      phi_all_cam_rep <- matrix(rep(phi_all_cams, dim(stay_time_data_all)[2]),
        nrow = ncam, ncol = dim(stay_time_data_all)[2]
      )

      # Alternate: fit staying time with gamma distribution
      # if(all(1/phi_star_cams>0) & all(!is.infinite(1/phi_star_cams))){
      #   mh1 <- sum(dgamma(stay_time_data_all,phi_star_cam_rep,log=TRUE),na.rm=TRUE) +
      #     sum(pgamma(stay_time_data_censor, phi_star_cam_rep, lower.tail = F, log = TRUE),na.rm=TRUE) +
      #     sum(dnorm(kappa_star,0,kappa_prior_var^0.5,log=TRUE))
      #   mh2 <- sum(dgamma(stay_time_data_all,phi_all_cam_rep,log=TRUE),na.rm=TRUE) +
      #     sum(pgamma(stay_time_data_censor, phi_all_cam_rep, lower.tail = F, log = TRUE),na.rm=TRUE) +
      #     sum(dnorm(kappa[i,],0,kappa_prior_var^0.5,log=TRUE))
      if (all(1 / phi_star_cams > 0) & all(!is.infinite(1 / phi_star_cams))) {
        mh1 <- sum(dexp(stay_time_data_all, 1 / phi_star_cam_rep, log = TRUE), na.rm = TRUE) +
          sum(pexp(stay_time_data_censor, 1 / phi_star_cam_rep, lower.tail = F, log = TRUE), na.rm = TRUE) +
          sum(dnorm(kappa_star, 0, kappa_prior_var^0.5, log = TRUE))
        mh2 <- sum(dexp(stay_time_data_all, 1 / phi_all_cam_rep, log = TRUE), na.rm = TRUE) +
          sum(pexp(stay_time_data_censor, 1 / phi_all_cam_rep, lower.tail = F, log = TRUE), na.rm = TRUE) +
          sum(dnorm(kappa[i, ], 0, kappa_prior_var^0.5, log = TRUE))
        mh <- exp(mh1 - mh2)

        if (mh > runif(1) & !is.na(mh)) {
          kappa[i, ] <- kappa_star
          accept[i + 1, num_covariates + kk] <- 1
          phi <- phi_star
        } else {
          kappa[i, ] <- kappa[i, ]
          accept[i + 1, num_covariates + kk] <- 0
        }
      } else {
        kappa[i, ] <- kappa[i, ]
        accept[i + 1, num_covariates + kk] <- 0
      }
    }
    kappa[i + 1, ] <- kappa[i, ]
    phi <- exp_na_covs(Z, kappa[i + 1, ])

    tot_u[i + 1] <- sum(u, na.rm = T) * cell_A

    # Update tuning parms
    if (i %% tune_check == 0) {
      batch_n <- batch_n + 1
      delta_n <- batch_n^-1
      # delta_n <- min(0.01,batch_n^-1)
      accept_gamma_check <- colMeans(
        array(
          accept[
            (i - tune_check + 1):i,
            1:num_covariates
          ],
          dim = c(tune_check, num_covariates)
        ),
        na.rm = T
      )
      gamma_tune[which(accept_gamma_check > 0.44)] <- gamma_tune[which(accept_gamma_check > 0.44)] + delta_n
      gamma_tune[which(accept_gamma_check <= 0.44)] <- gamma_tune[which(accept_gamma_check <= 0.44)] - delta_n
      accept_kappa_check <- colMeans(
        array(
          accept[
            (i - tune_check + 1):i,
            (num_covariates + 1):(2 * num_covariates)
          ],
          dim = c(tune_check, num_covariates)
        ),
        na.rm = T
      )
      kappa_tune[which(accept_kappa_check > 0.44)] <- kappa_tune[which(accept_kappa_check > 0.44)] + delta_n
      kappa_tune[which(accept_kappa_check <= 0.44)] <- kappa_tune[which(accept_kappa_check <= 0.44)] - delta_n
    }
  }
  # print("MCMC complete")

  list(accept = accept, gamma = gamma, kappa = kappa, tot_u = tot_u, u = u)
}


########################################
# MCMC for TTE no covariates
########################################
################################################################################
#' @export
#'
fit.model.mcmc.TTE <- function(study_design,
                               cam_design,
                               gamma_start,
                               kappa_start,
                               gamma_prior_var = 10^6,
                               kappa_prior_var = 10^6,
                               gamma_tune = -1,
                               kappa_tune = -1,
                               TTE_data_in,
                               stay_time_data_in
                               ) {
  n_iter <- study_design$n_iter
  censor <- study_design$t_steps
  TTE_censor <- study_design$TTE_censor
  cam_A <- cam_design$cam_A
  tot_A <- study_design$tot_A
  num_occ <- study_design$num_occ

  # Variables that will be saved
  gamma <- matrix(NA, n_iter + 1, 1)
  kappa <- matrix(NA, n_iter + 1, 1)
  tot_u <- matrix(NA, n_iter + 1, 1)
  accept <- matrix(NA, n_iter + 1, 2)

  gamma[1] <- gamma_start
  colnames(gamma) <- "gamma"
  kappa[1] <- kappa_start
  colnames(kappa) <- "kappa"
  colnames(tot_u) <- "Total estimate"
  colnames(accept) <- c("accept.rate.gamma", "accept.rate.kappa")

  # Define censored time to events
  TTE_data_in[TTE_data_in == TTE_censor] <- NA
  data.censor <- TTE_data_in
  data.censor[] <- 0
  data.censor[which(is.na(TTE_data_in))] <- TTE_censor

  # Account for censored times
  stay_time_data_all <- stay_time_data_in
  stay_time_data_all[stay_time_data_all > censor] <- NA
  stay_time_data_censor <- stay_time_data_in
  stay_time_data_censor[stay_time_data_censor <= censor] <- NA

  u <- exp(gamma[1])
  phi <- exp(kappa[1])
  tot_u <- u * tot_A
  eta <- u * cam_A / phi

  tune_check <- 100
  batch_n <- 0

  # Begin MCMC loop
  # prog_bar <- txtProgressBar(min = 0,max = n_iter,style = 3,width = 50,char = "=")
  for (i in 1:n_iter) {
    # setTxtProgressBar(prog_bar, i)
    # Sample gamma
    gamma_star <- rnorm(1, gamma[i], exp(2 * gamma_tune))
    u_star <- exp(gamma_star)
    eta_star <- u_star * cam_A / phi

    if (all(eta_star > 0)) {
      mh1 <- sum(dexp(TTE_data_in, eta_star, log = TRUE), na.rm = TRUE) +
        sum(pexp(data.censor, eta_star, lower.tail = F, log = TRUE)) +
        sum(dnorm(gamma_star, 0, gamma_prior_var^0.5, log = TRUE))
      mh2 <- sum(dexp(TTE_data_in, eta, log = TRUE), na.rm = TRUE) +
        sum(pexp(data.censor, eta, lower.tail = F, log = TRUE)) +
        sum(dnorm(gamma[i], 0, gamma_prior_var^0.5, log = TRUE))
      mh <- exp(mh1 - mh2)

      if (mh > runif(1) & !is.na(mh)) {
        gamma[i + 1] <- gamma_star
        accept[i + 1, 1] <- 1
        eta <- eta_star
      } else {
        gamma[i + 1] <- gamma[i]
        accept[i + 1, 1] <- 0
      }
    } else {
      gamma[i + 1] <- gamma[i]
      accept[i + 1, 1] <- 0
    }

    # Sample kappa
    kappa_star <- rnorm(1, kappa[i], exp(2 * kappa_tune))
    phi_star <- exp(kappa_star)

    if (1 / phi_star > 0 & !is.infinite(1 / phi_star)) {
      mh1 <- sum(dexp(stay_time_data_all, 1 / phi_star, log = TRUE), na.rm = TRUE) +
        sum(pexp(stay_time_data_censor, 1 / phi_star, lower.tail = F, log = TRUE), na.rm = TRUE) +
        sum(dnorm(kappa_star, 0, kappa_prior_var^0.5, log = TRUE))
      mh2 <- sum(dexp(stay_time_data_all, 1 / phi, log = TRUE), na.rm = TRUE) +
        sum(pexp(stay_time_data_censor, 1 / phi, lower.tail = F, log = TRUE), na.rm = TRUE) +
        sum(dnorm(kappa[i], 0, kappa_prior_var^0.5, log = TRUE))
      mh <- exp(mh1 - mh2)

      if (mh > runif(1) & !is.na(mh)) {
        kappa[i + 1] <- kappa_star
        accept[i + 1, 2] <- 1
        phi <- phi_star
      } else {
        kappa[i + 1] <- kappa[i]
        accept[i + 1, 2] <- 0
      }
    } else {
      kappa[i + 1] <- kappa[i]
      accept[i + 1, 2] <- 0
    }

    u <- exp(gamma[i + 1])
    phi <- exp(kappa[i + 1])
    tot_u[i + 1] <- u * tot_A
    eta <- u * cam_A / phi

    # Update tuning parms
    if (i %% tune_check == 0) {
      batch_n <- batch_n + 1
      delta_n <- batch_n^-1
      # delta_n <- min(0.01,batch_n^-1)
      accept_gamma_check <- mean(accept[(i - tune_check + 1):i, 1], na.rm = T)
      gamma_tune[accept_gamma_check > 0.44] <- gamma_tune[accept_gamma_check > 0.44] + delta_n
      gamma_tune[accept_gamma_check <= 0.44] <- gamma_tune[accept_gamma_check <= 0.44] - delta_n
      accept_kappa_check <- mean(accept[(i - tune_check + 1):i, 2], na.rm = T)
      kappa_tune[accept_kappa_check > 0.44] <- kappa_tune[accept_kappa_check > 0.44] + delta_n
      kappa_tune[accept_kappa_check <= 0.44] <- kappa_tune[accept_kappa_check <= 0.44] - delta_n
    }
  }
  # print("MCMC complete")

  list(accept = accept, gamma = gamma, kappa = kappa, tot_u = tot_u)
}


########################################
# MCMC for TTE w/ covariates
########################################
################################################################################
#' @export
#'
fit.model.mcmc.TTE.cov <- function(study_design,
                                   cam_design,
                                   cam_locs,
                                   gamma_start,
                                   kappa_start,
                                   gamma_prior_var = 10^6,
                                   kappa_prior_var = 10^6,
                                   gamma_tune,
                                   kappa_tune,
                                   TTE_data_in,
                                   stay_time_data_in,
                                   covariate_labels,
                                   num_covariates
                                   ) {
  n_iter <- study_design$n_iter
  censor <- study_design$t_steps
  TTE_censor <- study_design$TTE_censor
  cam_A <- cam_design$cam_A
  cell_A <- study_design$dx * study_design$dy
  num_occ <- study_design$num_occ
  ncam <- cam_design$ncam
  covariate_labels <- unlist(study_design$covariate_labels)
  num_covariates <- study_design$num_covariates
  Z <- matrix(unlist(study_design$Z), study_design$q, study_design$num_covariates)

  # Variables that will be saved
  gamma <- matrix(NA, n_iter + 1, num_covariates)
  kappa <- matrix(NA, n_iter + 1, num_covariates)
  tot_u <- matrix(NA, n_iter + 1, 1)
  accept <- matrix(NA, n_iter + 1, 2 * num_covariates)
  gamma[1, ] <- gamma_start
  colnames(gamma) <- paste0("gamma.", covariate_labels)
  kappa[1, ] <- kappa_start
  colnames(kappa) <- paste0("kappa.", covariate_labels)
  colnames(tot_u) <- "Total estimate"
  colnames(accept) <- c(
    paste0("accept.rate.gamma.", covariate_labels),
    paste0("accept.rate.kappa.", covariate_labels)
  )

  # Define censored time to events
  TTE_data_in[TTE_data_in == TTE_censor] <- NA
  data.censor <- TTE_data_in
  data.censor[] <- 0
  data.censor[which(is.na(TTE_data_in))] <- TTE_censor

  # Account for censored staying times
  stay_time_data_all <- stay_time_data_in
  stay_time_data_all[stay_time_data_all > censor] <- NA
  stay_time_data_censor <- stay_time_data_in
  stay_time_data_censor[stay_time_data_censor <= censor] <- NA

  # Initialize with landscape-scale covariates
  u <- exp_na_covs(Z, gamma[1, ])
  phi <- exp_na_covs(Z, kappa[1, ])
  tot_u[1] <- sum(u, na.rm = T) * cell_A
  eta <- u / phi * cell_A

  tune_check <- 100
  batch_n <- 0

  # Begin MCMC loop
  # prog_bar <- txtProgressBar(min = 0,max = n_iter,style = 3,width = 50,char = "=")
  for (i in 1:n_iter) {
    # setTxtProgressBar(prog_bar, i)
    # Sample gamma
    for (gg in 1:num_covariates) {
      gamma_star <- gamma[i, ]
      gamma_star[gg] <- rnorm(1, gamma[i, gg], exp(2 * gamma_tune[gg]))
      u_star <- exp_na_covs(Z, gamma_star)
      eta_star <- u_star / phi * cell_A

      # repeat estimated parms for fitting
      eta_star_cams <- eta_star[cam_locs$lscape_index]
      eta_star_cam_rep <- matrix(rep(eta_star_cams, num_occ), nrow = ncam, ncol = num_occ)
      eta_all_cams <- eta[cam_locs$lscape_index]
      eta_all_cam_rep <- matrix(rep(eta_all_cams, num_occ), nrow = ncam, ncol = num_occ)

      if (all(c(eta_star_cam_rep, eta_all_cam_rep) > 1e-10)) {
        mh1 <- sum(dexp(TTE_data_in, eta_star_cam_rep, log = TRUE), na.rm = TRUE) +
          sum(pexp(data.censor, eta_star_cam_rep, lower.tail = F, log = TRUE)) +
          sum(dnorm(gamma_star, 0, gamma_prior_var^0.5, log = TRUE))
        mh2 <- sum(dexp(TTE_data_in, eta_all_cam_rep, log = TRUE), na.rm = TRUE) +
          sum(pexp(data.censor, eta_all_cam_rep, lower.tail = F, log = TRUE)) +
          sum(dnorm(gamma[i, ], 0, gamma_prior_var^0.5, log = TRUE))
        mh <- exp(mh1 - mh2)

        if (mh > runif(1) & !is.na(mh)) {
          gamma[i, ] <- gamma_star
          accept[i + 1, gg] <- 1
          u <- u_star
          eta <- eta_star
        } else {
          gamma[i, ] <- gamma[i, ]
          accept[i + 1, gg] <- 0
        }
      } else {
        gamma[i, ] <- gamma[i, ]
        accept[i + 1, gg] <- 0
      }
    }
    gamma[i + 1, ] <- gamma[i, ]
    u <- exp_na_covs(Z, gamma[i + 1, ])
    eta <- u / phi * cell_A

    # Sample kappa
    for (kk in 1:num_covariates) {
      kappa_star <- kappa[i, ]
      kappa_star[kk] <- rnorm(1, kappa[i, kk], exp(2 * kappa_tune[kk]))
      phi_star <- exp_na_covs(Z, kappa_star)

      # repeat estimated parms for fitting
      phi_star_cams <- phi_star[cam_locs$lscape_index] * cam_A / cell_A
      phi_star_cam_rep <- matrix(rep(phi_star_cams, dim(stay_time_data_all)[2]), nrow = ncam, ncol = dim(stay_time_data_all)[2])
      phi_all_cams <- phi[cam_locs$lscape_index] * cam_A / cell_A
      phi_all_cam_rep <- matrix(rep(phi_all_cams, dim(stay_time_data_all)[2]), nrow = ncam, ncol = dim(stay_time_data_all)[2])

      # Alternate: fit staying time with gamma distribution
      # if(all(1/phi_star_cams>0) & all(!is.infinite(1/phi_star_cams))){
      # mh1 <- sum(dgamma(stay_time_data_all,phi_star_cam_rep,log=TRUE),na.rm=TRUE) +
      #   sum(pgamma(stay_time_data_censor, phi_star_cam_rep, lower.tail = F, log = TRUE),na.rm=TRUE) +
      #   sum(dnorm(kappa_star,0,kappa_prior_var^0.5,log=TRUE))
      # mh2 <- sum(dgamma(stay_time_data_all,phi_all_cam_rep,log=TRUE),na.rm=TRUE) +
      #   sum(pgamma(stay_time_data_censor, phi_all_cam_rep, lower.tail = F, log = TRUE),na.rm=TRUE) +
      #   sum(dnorm(kappa[i,],0,kappa_prior_var^0.5,log=TRUE))
      # mh <- exp(mh1-mh2)
      if (all(1 / phi_star_cams > 0) & all(!is.infinite(1 / phi_star_cams))) {
        mh1 <- sum(dexp(stay_time_data_all, 1 / phi_star_cam_rep, log = TRUE), na.rm = TRUE) +
          sum(pexp(stay_time_data_censor, 1 / phi_star_cam_rep, lower.tail = F, log = TRUE), na.rm = TRUE) +
          sum(dnorm(kappa_star, 0, kappa_prior_var^0.5, log = TRUE))
        mh2 <- sum(dexp(stay_time_data_all, 1 / phi_all_cam_rep, log = TRUE), na.rm = TRUE) +
          sum(pexp(stay_time_data_censor, 1 / phi_all_cam_rep, lower.tail = F, log = TRUE), na.rm = TRUE) +
          sum(dnorm(kappa[i, ], 0, kappa_prior_var^0.5, log = TRUE))
        mh <- exp(mh1 - mh2)

        if (is.infinite(mh)) {
          hh <- 0
          # message("Warning: Inf vals")
        }
        if (is.na(mh)) {
          message("ERROR: na vals")
        }

        if (mh > runif(1) & !is.na(mh)) {
          kappa[i, ] <- kappa_star
          accept[i + 1, num_covariates + kk] <- 1
          phi <- phi_star
        } else {
          kappa[i, ] <- kappa[i, ]
          accept[i + 1, num_covariates + kk] <- 0
        }
      } else {
        kappa[i, ] <- kappa[i, ]
        accept[i + 1, num_covariates + kk] <- 0
      }
    }
    kappa[i + 1, ] <- kappa[i, ]
    phi <- exp_na_covs(Z, kappa[i + 1, ])

    tot_u[i + 1] <- sum(u, na.rm = T) * cell_A

    # Update tuning parms
    if (i %% tune_check == 0) {
      batch_n <- batch_n + 1
      delta_n <- batch_n^-1
      # delta_n <- min(0.01,batch_n^-1)
      accept_gamma_check <- colMeans(
        array(
          accept[
            (i - tune_check + 1):i,
            1:num_covariates
          ],
          dim = c(tune_check, num_covariates)
        ),
        na.rm = T
      )
      gamma_tune[which(accept_gamma_check > 0.44)] <- gamma_tune[which(accept_gamma_check > 0.44)] + delta_n
      gamma_tune[which(accept_gamma_check <= 0.44)] <- gamma_tune[which(accept_gamma_check <= 0.44)] - delta_n
      accept_kappa_check <- colMeans(
        array(
          accept[
            (i - tune_check + 1):i,
            (num_covariates + 1):(2 * num_covariates)
          ],
          dim = c(tune_check, num_covariates)
        ),
        na.rm = T
      )
      kappa_tune[which(accept_kappa_check > 0.44)] <- kappa_tune[which(accept_kappa_check > 0.44)] + delta_n
      kappa_tune[which(accept_kappa_check <= 0.44)] <- kappa_tune[which(accept_kappa_check <= 0.44)] - delta_n
    }
  }
  # print("MCMC complete")

  list(accept = accept, gamma = gamma, kappa = kappa, tot_u = tot_u, u = u)
}



########################################
# MCMC for PR no covariates
########################################
################################################################################
#' @export
#'
fit.model.mcmc.PR <- function(study_design,
                              cam_design,
                              gamma_start,
                              gamma_prior_var = 10^6,
                              gamma_tune = -1,
                              count_data_in
                              ) {

  n_iter <- study_design$n_iter
  cam_A <- cam_design$cam_A
  tot_A <- study_design$tot_A
  t_steps <- study_design$t_steps

  # Variables that will be saved
  gamma <- matrix(NA, n_iter + 1, 1)
  tot_u <- matrix(NA, n_iter + 1, 1)
  accept <- matrix(NA, n_iter + 1, 1)

  gamma[1] <- gamma_start
  colnames(gamma) <- "gamma"
  colnames(tot_u) <- "Total estimate"
  colnames(accept) <- "accept.rate.gamma"

  d <- exp(gamma[1])
  tot_u <- d * tot_A / (cam_A * t_steps)

  tune_check <- 100
  batch_n <- 0

  # Begin MCMC loop
  # prog_bar <- txtProgressBar(min = 0,max = n_iter,style = 3,width = 50,char = "=")
  for (i in 1:n_iter) {
    # setTxtProgressBar(prog_bar, i)
    # Sample gamma
    gamma_star <- rnorm(1, gamma[i], exp(2 * gamma_tune))
    d_star <- exp(gamma_star)

    if (all(d_star > 0)) {
      mh1 <- sum(dpois(count_data_in, d_star, log = TRUE), na.rm = TRUE) +
        sum(dnorm(gamma_star, 0, gamma_prior_var^0.5, log = TRUE))
      mh2 <- sum(dpois(count_data_in, d, log = TRUE), na.rm = TRUE) +
        sum(dnorm(gamma[i], 0, gamma_prior_var^0.5, log = TRUE))
      mh <- exp(mh1 - mh2)

      if (mh > runif(1) & !is.na(mh)) {
        gamma[i + 1] <- gamma_star
        accept[i + 1] <- 1
        d <- d_star
      } else {
        gamma[i + 1] <- gamma[i]
        accept[i + 1] <- 0
      }
    } else {
      gamma[i + 1] <- gamma[i]
      accept[i + 1] <- 0
    }

    d <- exp(gamma[i + 1])
    tot_u[i + 1] <- d * tot_A / (cam_A * t_steps)

    # Update tuning parms
    if (i %% tune_check == 0) {
      batch_n <- batch_n + 1
      delta_n <- batch_n^-1
      # delta_n <- min(0.01,batch_n^-1)
      accept_gamma_check <- mean(accept[(i - tune_check + 1):i], na.rm = T)
      gamma_tune[accept_gamma_check > 0.44] <- gamma_tune[accept_gamma_check > 0.44] + delta_n
      gamma_tune[accept_gamma_check <= 0.44] <- gamma_tune[accept_gamma_check <= 0.44] - delta_n
    }
  }
  # print("MCMC complete")

  list(accept = accept, gamma = gamma, tot_u = tot_u)
}


########################################
# MCMC for PR w/ covariates
########################################
################################################################################
#' @export
#'
fit.model.mcmc.PR.cov <- function(study_design,
                                  cam_design,
                                  cam_locs,
                                  gamma_start,
                                  gamma_prior_var = 10^6,
                                  gamma_tune,
                                  count_data_in,
                                  stay_time_data_in
                                  ) {

  n_iter <- study_design$n_iter
  cam_A <- cam_design$cam_A
  cell_A <- study_design$dx * study_design$dy
  t_steps <- study_design$t_steps
  covariate_labels <- unlist(study_design$covariate_labels)
  num_covariates <- study_design$num_covariates
  Z <- matrix(unlist(study_design$Z), study_design$q, study_design$num_covariates)

  # Variables that will be saved
  gamma <- matrix(NA, n_iter + 1, num_covariates)
  tot_u <- matrix(NA, n_iter + 1, 1)
  accept <- matrix(NA, n_iter + 1, num_covariates)
  gamma[1, ] <- gamma_start
  colnames(gamma) <- paste0("gamma.", covariate_labels)
  colnames(tot_u) <- "Total estimate"
  colnames(accept) <- c(paste0(
    "accept.rate.gamma.",
    covariate_labels
  ))
  # Initialize with landscape-scale covariates
  d <- exp_na_covs(Z, gamma[1, ]) # expected densities
  u <- d / t_steps
  tot_u[1] <- sum(u, na.rm = T)

  tune_check <- 100
  batch_n <- 0

  # Begin MCMC loop
  # prog_bar <- txtProgressBar(min = 0,max = n_iter,style = 3,width = 50,char = "=")
  for (i in 1:n_iter) {
    # setTxtProgressBar(prog_bar, i)
    # Sample gamma
    for (gg in 1:num_covariates) {
      gamma_star <- gamma[i, ]
      gamma_star[gg] <- rnorm(1, gamma[i, gg], exp(2 * gamma_tune[gg]))
      d_star <- exp_na_covs(Z, gamma_star)
      u_star <- d_star / t_steps

      # repeat estimated parms for fitting
      d_star_cams <- d_star[cam_locs$lscape_index] * cam_A / cell_A
      d_all_cams <- d[cam_locs$lscape_index] * cam_A / cell_A

      if (all(d_star_cams > 0)) {
        mh1 <- sum(dpois(count_data_in, d_star_cams, log = TRUE), na.rm = TRUE) +
          sum(dnorm(gamma_star, 0, gamma_prior_var^0.5, log = TRUE))
        mh2 <- sum(dpois(count_data_in, d_all_cams, log = TRUE), na.rm = TRUE) +
          sum(dnorm(gamma[i, ], 0, gamma_prior_var^0.5, log = TRUE))
        mh <- exp(mh1 - mh2)

        if (mh > runif(1) & !is.na(mh)) {
          gamma[i, ] <- gamma_star
          accept[i + 1, gg] <- 1
          u <- u_star
          d <- d_star
        } else {
          gamma[i, ] <- gamma[i, ]
          accept[i + 1, gg] <- 0
        }
      } else {
        gamma[i, ] <- gamma[i, ]
        accept[i + 1, gg] <- 0
      }
    }
    gamma[i + 1, ] <- gamma[i, ]
    d <- exp_na_covs(Z, gamma[i + 1, ])
    u <- d / t_steps
    tot_u[i + 1] <- sum(u, na.rm = T)

    # Update tuning parms
    if (i %% tune_check == 0) {
      batch_n <- batch_n + 1
      delta_n <- batch_n^-1
      # delta_n <- min(0.01,batch_n^-1)
      accept_gamma_check <- colMeans(
        array(
          accept[
            (i - tune_check + 1):i,
            1:num_covariates
          ],
          dim = c(tune_check, num_covariates)
        ),
        na.rm = T
      )
      gamma_tune[which(accept_gamma_check > 0.44)] <- gamma_tune[which(accept_gamma_check > 0.44)] + delta_n
      gamma_tune[which(accept_gamma_check <= 0.44)] <- gamma_tune[which(accept_gamma_check <= 0.44)] - delta_n
    }
  }
  # print("MCMC complete")

  list(accept = accept, gamma = gamma, tot_u = tot_u, u = u)
}

########################################
# MCMC for STE no covariates
########################################
################################################################################
#' @export
#'
fit.model.mcmc.STE <- function(study_design,
                               cam_design,
                               gamma_start,
                               gamma_prior_var = 10^6,
                               gamma_tune = -1,
                               STE_data_in
                               ) {

  n_iter <- study_design$n_iter
  censor <- cam_design$ncam * cam_design$cam_A
  tot_A <- study_design$tot_A
  cam_A <- cam_design$cam_A

  # Variables that will be saved
  gamma <- matrix(NA, n_iter + 1, 1)
  tot_u <- matrix(NA, n_iter + 1, 1)
  accept <- matrix(NA, n_iter + 1, 1)

  gamma[1] <- gamma_start
  colnames(gamma) <- "gamma"
  colnames(tot_u) <- "Total estimate"
  colnames(accept) <- "accept.rate.gamma"

  # Define censored space to event
  STE_data_in[STE_data_in == censor] <- NA
  STE.censor <- STE_data_in
  STE.censor[] <- 0
  STE.censor[which(is.na(STE_data_in))] <- censor

  u <- exp(gamma[1])
  tot_u <- u * tot_A

  tune_check <- 100
  batch_n <- 0

  # Begin MCMC loop
  # prog_bar <- txtProgressBar(min = 0,max = n_iter,style = 3,width = 50,char = "=")
  for (i in 1:n_iter) {
    # setTxtProgressBar(prog_bar, i)
    # Sample gamma
    gamma_star <- rnorm(1, gamma[i], exp(2 * gamma_tune))
    u_star <- exp(gamma_star)

    if (all(u_star > 0)) {
      mh1 <- sum(dexp(STE_data_in, u_star, log = TRUE), na.rm = TRUE) +
        sum(pexp(STE.censor, u_star, lower.tail = F, log = TRUE)) +
        sum(dnorm(gamma_star, 0, gamma_prior_var^0.5, log = TRUE))
      mh2 <- sum(dexp(STE_data_in, u, log = TRUE), na.rm = TRUE) +
        sum(pexp(STE.censor, u, lower.tail = F, log = TRUE)) +
        sum(dnorm(gamma[i], 0, gamma_prior_var^0.5, log = TRUE))
      mh <- exp(mh1 - mh2)

      if (mh > runif(1) & !is.na(mh)) {
        gamma[i + 1] <- gamma_star
        accept[i + 1] <- 1
        u <- u_star
      } else {
        gamma[i + 1] <- gamma[i]
        accept[i + 1] <- 0
      }
    } else {
      gamma[i + 1] <- gamma[i]
      accept[i + 1] <- 0
    }

    u <- exp(gamma[i + 1])
    tot_u[i + 1] <- u * tot_A

    # Update tuning parms
    if (i %% tune_check == 0) {
      batch_n <- batch_n + 1
      delta_n <- batch_n^-1
      # delta_n <- min(0.01,batch_n^-1)
      accept_gamma_check <- mean(accept[(i - tune_check + 1):i], na.rm = T)
      gamma_tune[accept_gamma_check > 0.44] <- gamma_tune[accept_gamma_check > 0.44] + delta_n
      gamma_tune[accept_gamma_check <= 0.44] <- gamma_tune[accept_gamma_check <= 0.44] - delta_n
    }
  }
  # print("MCMC complete")

  list(accept = accept, gamma = gamma, tot_u = tot_u)
}


########################################
# MCMC for SECR no covariates
########################################
################################################################################
#' @export
#'
fit.model.mcmc.SECR <- function(SECR.dat,
                                psi.start,
                                sigma.start,
                                lambda_0.start,
                                S.tune,
                                sigma.tune,
                                lambda_0.tune,
                                X,
                                M,
                                x,
                                y,
                                n_iter) {
  # Variables that will be saved
  psi <- matrix(NA, n_iter + 1, 1)
  sigma <- matrix(NA, n_iter + 1, 1)
  lambda_0 <- matrix(NA, n_iter + 1, 1)
  tot_u <- matrix(NA, n_iter + 1, 1)
  accept <- matrix(NA, n_iter + 1, 2)

  psi[1] <- psi.start
  sigma[1] <- sigma.start
  lambda_0[1] <- lambda_0.start
  accept[1, ] <- 0
  colnames(psi) <- "psi"
  colnames(sigma) <- "sigma"
  colnames(lambda_0) <- "lambda_0"
  colnames(tot_u) <- "Total estimate"
  colnames(accept) <- c("accept.rate.sigma", "accept.rate.lambda_0")

  tune_check <- 50
  batch_n <- 10

  # create initial values
  # home ranges
  S <- cbind(runif(M, x[1], x[2]), runif(M, y[1], y[2]))
  z <- rbinom(M, 1, psi[1])
  seen <- apply(SECR.dat > 0, 1, any)
  z[seen] <- 1 # set seen individuals' z=1
  tot_u[1] <- sum(z)

  # initiate distance matrix and lambda_ij matrix
  d <- e2dist(S, X)
  lambda <- lambda_0[1] * exp(-(d^2) / (2 * sigma[1]^2))

  # start iterations of the chain
  for (iter in 1:n_iter) {
    # update sigma
    sigma.star <- rnorm(1, sigma[iter], exp(2 * sigma.tune))
    if (sigma.star > 0 & sigma.star < 10) { # automatically reject sigma.star that are < 10
      lambda.star <- lambda_0[iter] * exp(-(d^2) / (2 * sigma.star^2))
      mh1 <- sum(dpois(SECR.dat, lambda.star * z, log = TRUE))
      mh2 <- sum(dpois(SECR.dat, lambda * z, log = TRUE))
      mh <- exp(mh1 - mh2)
      if (runif(1) < mh & !is.na(mh)) {
        lambda <- lambda.star
        sigma[iter + 1] <- sigma.star
        accept[iter + 1, 1] <- 1
      } else {
        sigma[iter + 1] <- sigma[iter]
        accept[iter + 1, 1] <- 0
      }
    } else {
      sigma[iter + 1] <- sigma[iter]
      accept[iter + 1, 1] <- 0
    }

    # update lambda_0
    lambda_0.star <- rnorm(1, lambda_0[iter], exp(2 * lambda_0.tune))
    if (lambda_0.star > 0) { # automatically reject lambda_0.star that are < 10
      lambda.star <- lambda_0.star * exp(-(d^2) / (2 * sigma[iter + 1]^2))
      mh1 <- sum(dpois(SECR.dat, lambda.star * z, log = TRUE))
      mh2 <- sum(dpois(SECR.dat, lambda * z, log = TRUE))
      mh <- exp(mh1 - mh2)
      if (runif(1) < mh & !is.na(mh)) {
        lambda <- lambda.star
        lambda_0[iter + 1] <- lambda_0.star
        accept[iter + 1, 2] <- 1
      } else {
        lambda_0[iter + 1] <- lambda_0[iter]
        accept[iter + 1, 2] <- 0
      }
    } else {
      lambda_0[iter + 1] <- lambda_0[iter]
      accept[iter + 1, 2] <- 0
    }

    # update z
    Zups <- 0
    for (i in 1:M) {
      if (!seen[i]) { # no need to update seen individuals, since their z =1
        z.star <- ifelse(z[i] == 0, 1, 0)
        mh1 <- sum(dpois(SECR.dat[i, ], lambda[i, ] * z.star, log = TRUE)) +
          dbinom(z.star, 1, psi[iter], log = TRUE)
        mh2 <- sum(dpois(SECR.dat[i, ], lambda[i, ] * z[i], log = TRUE)) +
          dbinom(z[i], 1, psi[iter], log = TRUE)
        mh <- exp(mh1 - mh2)

        if (runif(1) < mh & !is.na(mh)) {
          z[i] <- z.star
          Zups <- Zups + 1
        }
      }
    } # end M loop
    # print(mean(Zups))

    # update psi
    psi[iter + 1] <- rbeta(1, 1 + sum(z), 1 + M - sum(z))

    # update S
    Ssups <- 0
    for (i in 1:M) {
      S.star <- c(rnorm(1, S[i, 1], S.tune), rnorm(1, S[i, 2], S.tune))

      # Check to make sure S.star is within bounds
      inbox <-
        S.star[1] < x[2] &
          S.star[1] > x[1] &
          S.star[2] < y[2] &
          S.star[2] > y[1]

      if (inbox) {
        dtmp <- sqrt((S.star[1] - X[, 1])^2 + (S.star[2] - X[, 2])^2)
        lambda.star <- lambda_0[iter + 1] * exp(-(dtmp^2) / (2 * sigma[iter + 1]^2))

        mh1 <- sum(dpois(SECR.dat[i, ], lambda.star * z[i], log = TRUE))
        mh2 <- sum(dpois(SECR.dat[i, ], lambda[i, ] * z[i], log = TRUE))
        mh <- exp(mh1 - mh2)
        if (runif(1) < mh & !is.na(mh)) {
          S[i, ] <- S.star
          lambda[i, ] <- lambda.star
          d[i, ] <- dtmp
          Ssups <- Ssups + 1
        }
      }
    } # end M loop
    # print(Ssups)

    tot_u[iter + 1] <- sum(z)

    # Update tuning parms
    if (iter %% tune_check == 0) {
      batch_n <- batch_n + 1
      delta_n <- batch_n^-1
      # delta_n <- min(0.01,batch_n^-1)
      accept.sigma.check <- mean(accept[(iter - tune_check + 1):iter, 1], na.rm = T)
      sigma.tune[accept.sigma.check > 0.44] <- sigma.tune[accept.sigma.check > 0.44] + delta_n
      sigma.tune[accept.sigma.check <= 0.44] <- sigma.tune[accept.sigma.check <= 0.44] - delta_n
      accept.lambda_0.check <- mean(accept[(iter - tune_check + 1):iter, 2], na.rm = T)
      lambda_0.tune[accept.lambda_0.check > 0.44] <- lambda_0.tune[accept.lambda_0.check > 0.44] + delta_n
      lambda_0.tune[accept.lambda_0.check <= 0.44] <- lambda_0.tune[accept.lambda_0.check <= 0.44] - delta_n
    }
  }
  # print("MCMC complete")

  list(accept = accept, sigma = sigma, psi = psi, lambda_0 = lambda_0, tot_u = tot_u)
} # end of function call

################################################################################
#' @export
#'
SCR0pois <-
  function(y, X, M, xl, xu, yl, yu, delta, niter) {
    # create initial values
    S <- cbind(runif(M, xl, xu), runif(M, yl, yu))
    sigma <- runif(1, 0.5, 5)
    lam0 <- runif(1, 0.1, 1)
    psi <- runif(1, 0.2, 0.8)
    z <- rbinom(M, 1, psi)
    seen <- apply(y > 0, 1, any)
    z[seen] <- 1 # set seen individuals' z=1

    # initiate distance matrix and lamij matrix
    d <- e2dist(S, X)
    lam <- lam0 * exp(-(d * d) / (2 * sigma * sigma))

    # set up matrix to hold results
    out <- matrix(nrow = niter, ncol = 4)
    colnames(out) <- c("sigma", "lam0", "psi", "N")

    # have R print starting values
    cat("\nstarting values =", c(sigma, lam0, psi, sum(z), "\n\n"))

    # start iterations of the chain

    for (iter in 1:niter) {
      # have R output the time and parameter estimates at every 100th iteration
      if (iter %% 100 == 0) {
        cat("iter", iter, format(Sys.time(), "%H:%M:%S"), "\n")
        cat("   current =", out[iter - 1, ], "\n")
      }

      # update sigma
      sig.cand <- rnorm(1, sigma, delta[1])
      if (sig.cand > 0) { # automatically reject sig.cand that are <0
        lam.cand <- lam0 * exp(-(d * d) / (2 * sig.cand * sig.cand))
        ll <- sum(dpois(y, lam * z, log = TRUE))
        llcand <- sum(dpois(y, lam.cand * z, log = TRUE))
        if (runif(1) < exp(llcand - ll)) {
          ll <- llcand
          lam <- lam.cand
          sigma <- sig.cand
        }
      }

      # update lam0
      lam0.cand <- rnorm(1, lam0, delta[2])
      if (lam0.cand > 0) { # automatically reject lam0.cand that are <0
        lam.cand <- lam0.cand * exp(-(d * d) / (2 * sigma * sigma))
        ll <- sum(dpois(y, lam * z, log = TRUE))
        llcand <- sum(dpois(y, lam.cand * z, log = TRUE))
        if (runif(1) < exp(llcand - ll)) {
          ll <- llcand
          lam <- lam.cand
          lam0 <- lam0.cand
        }
      }

      # update z
      zUps <- 0 # set counter to monitor acceptance rate
      for (i in 1:M) {
        if (seen[i]) { # no need to update seen individuals, since their z =1
          next
        }
        zcand <- ifelse(z[i] == 0, 1, 0)
        llz <- sum(dpois(y[i, ], lam[i, ] * z[i], log = TRUE))
        llcand <- sum(dpois(y[i, ], lam[i, ] * zcand, log = TRUE))

        prior <- dbinom(z[i], 1, psi, log = TRUE)
        prior.cand <- dbinom(zcand, 1, psi, log = TRUE)
        if (runif(1) < exp((llcand + prior.cand) - (llz + prior))) {
          z[i] <- zcand
          zUps <- zUps + 1
        }
      } # end M loop

      # update psi
      psi <- rbeta(1, 1 + sum(z), 1 + M - sum(z))

      # update s
      Sups <- 0
      for (i in 1:M) {
        Scand <- c(rnorm(1, S[i, 1], delta[3]), rnorm(1, S[i, 2], delta[3]))
        inbox <- Scand[1] < xu & Scand[1] > xl & Scand[2] < yu & Scand[2] > yl
        if (inbox) {
          dtmp <- sqrt((Scand[1] - X[, 1])^2 + (Scand[2] - X[, 2])^2)
          lam.cand <- lam0 * exp(-(dtmp * dtmp) / (2 * sigma * sigma))

          llS <- sum(dpois(y[i, ], lam[i, ] * z[i], log = TRUE))
          llcand <- sum(dpois(y[i, ], lam.cand * z[i], log = TRUE))
          if (runif(1) < exp(llcand - llS)) {
            S[i, ] <- Scand
            lam[i, ] <- lam.cand
            d[i, ] <- dtmp
            Sups <- Sups + 1
          }
        }
      } # end M loop

      # prompt R to output acceptance rates of z and S
      if (iter %% 100 == 0) {
        cat("   Acceptance rates\n")
        cat("     z =", zUps / M, "\n")
        cat("     S =", Sups / M, "\n")
      }

      out[iter, ] <- c(sigma, lam0, psi, sum(z))
    } # end of iteration loop

    return(out)
  } # end of function call

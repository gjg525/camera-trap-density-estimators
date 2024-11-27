# Metropolis-within-Gibbs sampler
########################################
# MCMC for TDST w/ covariates
########################################
################################################################################
fit.model.mcmc.TDST.cov <- function(study_design,
                                    cam_design,
                                    cam_locs,
                                    gamma_start,
                                    kappa_start,
                                    gamma_prior_var = 10^4,
                                    kappa_prior_mu = 0,
                                    kappa_prior_var = 10^4,
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
  accept <- matrix(NA, n_iter + 1, 2)
  gamma[1] <- gamma_start
  colnames(gamma) <- "gamma"
  kappa[1, ] <- kappa_start
  colnames(kappa) <- paste0("kappa.", covariate_labels)
  colnames(tot_u) <- "Total estimate"
  colnames(accept) <- c(
    "accept.rate.gamma",
    # paste0("accept.rate.kappa.", covariate_labels)
    "accept.rate.kappa"
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
    kappa_star <- kappa[i, ]
    kappa_star <- rnorm(3, kappa[i, ], exp(2 * kappa_tune))

    phi_star <- exp_na_covs(Z, kappa_star)
    sum_phi_star <- sum(phi_star, na.rm = T)
    if (is.infinite(sum_phi_star)) {
      u_star <- phi_star * 0
    } else {
      u_star <- beta * phi_star / sum_phi_star
    }
    
    # repeat estimated parms for fitting
    u_star_cams <- u_star[cam_locs$lscape_index] * cam_A / cell_A
    u_cams <- u[cam_locs$lscape_index] * cam_A / cell_A
    
    if (is.null(stay_time_data_in)) {
      # No stay time data from cameras
      mh1 <- sum(dpois(count_data_in, u_star_cams, log = TRUE), na.rm = TRUE) +
        sum(dnorm(kappa_star,kappa_prior_mu,kappa_prior_var^0.5,log=TRUE))
      mh2 <- sum(dpois(count_data_in, u_cams, log = TRUE), na.rm = TRUE) +
        sum(dnorm(kappa[i,],kappa_prior_mu,kappa_prior_var^0.5,log=TRUE))
      mh <- exp(mh1-mh2)
    } else{
      # # repeat estimated parms for fitting
      phi_star_cams <- phi_star[cam_locs$lscape_index] * cam_A / cell_A
      phi_star_cam_rep <- matrix(rep(phi_star_cams, dim(stay_time_data_all)[2]),
                                 nrow = ncam, ncol = dim(stay_time_data_all)[2]
      )
      phi_all_cams <- phi[cam_locs$lscape_index] * cam_A / cell_A
      phi_all_cam_rep <- matrix(rep(phi_all_cams, dim(stay_time_data_all)[2]),
                                nrow = ncam, ncol = dim(stay_time_data_all)[2]
      )
      
      # repeat estimated parms for fitting
      u_star_cams <- u_star[cam_locs$lscape_index] * cam_A / cell_A
      u_cams <- u[cam_locs$lscape_index] * cam_A / cell_A
      
      if (all(1 / phi_star_cams > 0) & all(!is.infinite(1 / phi_star_cams))) {
        mh1 <- sum(dexp(stay_time_data_all, 1 / phi_star_cam_rep, log = TRUE), na.rm = TRUE) +
          sum(pexp(stay_time_data_censor, 1 / phi_star_cam_rep, lower.tail = F, log = TRUE), na.rm = TRUE) +
          sum(dpois(count_data_in, u_star_cams, log = TRUE), na.rm = TRUE) +
          sum(dnorm(kappa_star, kappa_prior_mu, kappa_prior_var^0.5, log = TRUE))
        mh2 <- sum(dexp(stay_time_data_all, 1 / phi_all_cam_rep, log = TRUE), na.rm = TRUE) +
          sum(pexp(stay_time_data_censor, 1 / phi_all_cam_rep, lower.tail = F, log = TRUE), na.rm = TRUE) +
          sum(dpois(count_data_in, u_cams, log = TRUE), na.rm = TRUE) +
          sum(dnorm(kappa[i, ], kappa_prior_mu, kappa_prior_var^0.5, log = TRUE))
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
        
      } else {
        kappa[i, ] <- kappa[i, ]
        accept[i + 1, 2] <- 0
      }
    }
    if (mh > runif(1) & !is.na(mh)) {
      kappa[i, ] <- kappa_star
      accept[i + 1, 2] <- 1
      u <- u_star
      phi <- phi_star
    } else {
      kappa[i, ] <- kappa[i, ]
      accept[i + 1, 2] <- 0
    }
    
    # for (kk in 1:num_covariates) {
    #   kappa_star <- kappa[i, ]
    #   kappa_star[kk] <- rnorm(1, kappa[i, kk], exp(2 * kappa_tune[kk]))
    #   phi_star <- exp_na_covs(Z, kappa_star)
    #   sum_phi_star <- sum(phi_star, na.rm = T)
    #   if (is.infinite(sum_phi_star)) {
    #     u_star <- phi_star * 0
    #   } else {
    #     u_star <- beta * phi_star / sum_phi_star
    #   }
    # 
    #   if (is.null(stay_time_data_in)) {
    #     # No stay time data from cameras
    #     mh1 <- sum(dpois(count_data_in, u_star_cams, log = TRUE), na.rm = TRUE) +
    #       sum(dnorm(kappa_star,kappa_prior_mu,kappa_prior_var^0.5,log=TRUE))
    #     mh2 <- sum(dpois(count_data_in, u_cams, log = TRUE), na.rm = TRUE) +
    #       sum(dnorm(kappa[i,],kappa_prior_mu,kappa_prior_var^0.5,log=TRUE))
    #     mh <- exp(mh1-mh2)
    #   } else{
    #     # # repeat estimated parms for fitting
    #     phi_star_cams <- phi_star[cam_locs$lscape_index] * cam_A / cell_A
    #     phi_star_cam_rep <- matrix(rep(phi_star_cams, dim(stay_time_data_all)[2]),
    #       nrow = ncam, ncol = dim(stay_time_data_all)[2]
    #     )
    #     phi_all_cams <- phi[cam_locs$lscape_index] * cam_A / cell_A
    #     phi_all_cam_rep <- matrix(rep(phi_all_cams, dim(stay_time_data_all)[2]),
    #       nrow = ncam, ncol = dim(stay_time_data_all)[2]
    #     )
    # 
    #     # repeat estimated parms for fitting
    #     u_star_cams <- u_star[cam_locs$lscape_index] * cam_A / cell_A
    #     u_cams <- u[cam_locs$lscape_index] * cam_A / cell_A
    #     
    #     if (all(1 / phi_star_cams > 0) & all(!is.infinite(1 / phi_star_cams))) {
    #       mh1 <- sum(dexp(stay_time_data_all, 1 / phi_star_cam_rep, log = TRUE), na.rm = TRUE) +
    #         sum(pexp(stay_time_data_censor, 1 / phi_star_cam_rep, lower.tail = F, log = TRUE), na.rm = TRUE) +
    #         sum(dpois(count_data_in, u_star_cams, log = TRUE), na.rm = TRUE) +
    #         sum(dnorm(kappa_star, kappa_prior_mu, kappa_prior_var^0.5, log = TRUE))
    #       mh2 <- sum(dexp(stay_time_data_all, 1 / phi_all_cam_rep, log = TRUE), na.rm = TRUE) +
    #         sum(pexp(stay_time_data_censor, 1 / phi_all_cam_rep, lower.tail = F, log = TRUE), na.rm = TRUE) +
    #         sum(dpois(count_data_in, u_cams, log = TRUE), na.rm = TRUE) +
    #         sum(dnorm(kappa[i, ], kappa_prior_mu, kappa_prior_var^0.5, log = TRUE))
    #       mh <- exp(mh1 - mh2)
    # 
    #       # Alternate: fit staying time with gamma distribution
    #       # if(all(1/phi_star_cams>0) & all(!is.infinite(1/phi_star_cams))){
    #       #   mh1 <- sum(dgamma(stay_time_data_all, phi_star_cam_rep,log=TRUE),na.rm=TRUE) +
    #       #     sum(pgamma(stay_time_data_censor, phi_star_cam_rep, lower.tail = F, log = TRUE),na.rm=TRUE) +
    #       #     sum(dnorm(kappa_star,0,kappa_prior_var^0.5,log=TRUE))
    #       #   mh2 <- sum(dgamma(stay_time_data_all, phi_all_cam_rep,log=TRUE),na.rm=TRUE) +
    #       #     sum(pgamma(stay_time_data_censor, phi_all_cam_rep, lower.tail = F, log = TRUE),na.rm=TRUE) +
    #       #     sum(dnorm(kappa[i,],0,kappa_prior_var^0.5,log=TRUE))
    #       #   mh <- exp(mh1-mh2)
    # 
    #     } else {
    #       kappa[i, ] <- kappa[i, ]
    #       accept[i + 1, 1 + kk] <- 0
    #     }
    #   }
    #   if (mh > runif(1) & !is.na(mh)) {
    #     kappa[i, ] <- kappa_star
    #     accept[i + 1, 1 + kk] <- 1
    #     u <- u_star
    #     phi <- phi_star
    #   } else {
    #     kappa[i, ] <- kappa[i, ]
    #     accept[i + 1, 1 + kk] <- 0
    #   }
    # }
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
      
      accept_kappa_check <- mean(accept[(i - tune_check + 1):i, 2], na.rm = T)
      kappa_tune[accept_kappa_check > 0.44] <- kappa_tune[accept_kappa_check > 0.44] + delta_n
      kappa_tune[accept_kappa_check <= 0.44] <- kappa_tune[accept_kappa_check <= 0.44] - delta_n
      
      # accept_kappa_check <- colMeans(
      #   array(
      #     accept[
      #       (i - tune_check + 1):i,
      #       2:(1 + num_covariates)
      #     ],
      #     dim = c(tune_check, num_covariates)
      #   ),
      #   na.rm = T
      # )
      # kappa_tune[which(accept_kappa_check > 0.44)] <- kappa_tune[which(accept_kappa_check > 0.44)] + delta_n
      # kappa_tune[which(accept_kappa_check <= 0.44)] <- kappa_tune[which(accept_kappa_check <= 0.44)] - delta_n
    }
  }
  # print("MCMC complete")

  list(accept = accept, gamma = gamma, kappa = kappa, tot_u = tot_u, u = u)
}

########################################
# MCMC for PR no covariates
########################################
fit.model.mcmc.PR <- function(study_design,
                              cam_design,
                              gamma_start,
                              gamma_prior_var = 10^4,
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
fit.model.mcmc.PR.cov <- function(study_design,
                                  cam_design,
                                  cam_locs,
                                  gamma_start,
                                  gamma_prior_var = 10^4,
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
# MCMC for PR habitat
########################################
fit.model.mcmc.PR.habitat <- function(study_design,
                                      cam_design,
                                      cam_locs,
                                      gamma_start,
                                      gamma_prior_var = 10^4,
                                      gamma_tune,
                                      kappa_start,
                                      kappa_prior_mu,
                                      kappa_prior_var,
                                      kappa_tune,
                                      count_data_in,
                                      habitat_summary) {
  
  n_iter <- study_design$n_iter
  cam_A <- cam_design$cam_A
  cell_A <- study_design$dx * study_design$dy
  t_steps <- study_design$t_steps
  covariate_labels <- unlist(study_design$covariate_labels)
  num_covariates <- study_design$num_covariates
  Z <- matrix(unlist(study_design$Z), study_design$q, study_design$num_covariates)
  
  # Variables that will be saved
  gamma <- matrix(, n_iter + 1, num_covariates)
  kappa <- matrix(, n_iter + 1, num_covariates)
  tot_u <- matrix(, n_iter + 1, 1)
  # accept <- matrix(, n_iter + 1, 2 * num_covariates)
  accept <- matrix(, n_iter + 1, num_covariates + 1)
  gamma[1, ] <- gamma_start
  colnames(gamma) <- paste0("gamma.", covariate_labels)
  kappa[1, ] <- kappa_start
  colnames(kappa) <- paste0("kappa.", covariate_labels)
  colnames(tot_u) <- "Total estimate"
  # colnames(accept) <- c(
  #   paste0("accept.rate.gamma.", covariate_labels),
  #   paste0("accept.rate.kappa.", covariate_labels)
  # )
  colnames(accept) <- c(
    paste0("accept.rate.gamma.", covariate_labels),
    "accept.rate.kappa."
  )
  
  d <- exp(gamma[1, ])
  u <- d * habitat_summary$n_lscape / exp(kappa[1, ]) * 
    habitat_summary$prop_cams / (cam_design$cam_A * study_design$t_steps)
  
  tot_u <- sum(u)
  
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
      d_star <- exp(gamma_star)
     
      count_data_in_h <- count_data_in$count[count_data_in$Speed == covariate_labels[gg]]
      
      if (length(count_data_in_h) == 0) {
        count_data_in_h <- 0
      }
      
      if (all(d_star > 0)) {
        mh1 <- sum(dpois(count_data_in_h, d_star[gg], log = TRUE), na.rm = TRUE) +
          sum(dnorm(gamma_star[gg], 0, gamma_prior_var^0.5, log = TRUE))
        mh2 <- sum(dpois(count_data_in_h, d[gg], log = TRUE), na.rm = TRUE) +
          sum(dnorm(gamma[i, gg], 0, gamma_prior_var^0.5, log = TRUE))
        mh <- exp(mh1 - mh2)
        
        if (mh > runif(1)) {
          gamma[i, ] <- gamma_star
          accept[i + 1, gg] <- 1
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
    d <- exp(gamma[i + 1, ])
    
    # Sample kappa
    kappa_star <- kappa[i, ]
    kappa_star <- rnorm(3, kappa[i, ], exp(2 * kappa_tune))
    kappa_star <- log(exp(kappa_star) / sum(exp(kappa_star)))
    
    # Proportional staying time defined by priors
    mh1 <- sum(dnorm(kappa_star,kappa_prior_mu,kappa_prior_var^0.5,log=TRUE))
    mh2 <- sum(dnorm(kappa[i,],kappa_prior_mu,kappa_prior_var^0.5,log=TRUE))
    mh <- exp(mh1-mh2)
    
    if (mh > runif(1) & !is.na(mh)) {
      kappa[i, ] <- kappa_star
      accept[i + 1, num_covariates + 1] <- 1
      # u <- u_star
    } else {
      kappa[i, ] <- kappa[i, ]
      accept[i + 1, num_covariates + 1] <- 0
    }
    
    # for (kk in 1:num_covariates) {
    #   kappa_star <- kappa[i, ]
    #   kappa_star[kk] <- rnorm(1, kappa[i, kk], exp(2 * kappa_tune[kk]))
    #   kappa_star <- log(exp(kappa_star) / sum(exp(kappa_star)))
    #   
    #   # u_star <- d * habitat_summary$n_lscape * habitat_summary$prop_cams * 
    #   #   exp(kappa_star) / (cam_design$cam_A * study_design$t_steps)
    #   
    #   # Proportional staying time defined by priors
    #   mh1 <- sum(dnorm(kappa_star,kappa_prior_mu,kappa_prior_var^0.5,log=TRUE))
    #   mh2 <- sum(dnorm(kappa[i,],kappa_prior_mu,kappa_prior_var^0.5,log=TRUE))
    #   mh <- exp(mh1-mh2)
    #   
    #   if (mh > runif(1) & !is.na(mh)) {
    #     kappa[i, ] <- kappa_star
    #     accept[i + 1, num_covariates + kk] <- 1
    #     # u <- u_star
    #   } else {
    #     kappa[i, ] <- kappa[i, ]
    #     accept[i + 1, num_covariates + kk] <- 0
    #   }
    #   
    # }
    
    kappa[i + 1, ] <- kappa[i, ]

    u <- d * habitat_summary$n_lscape * 
      habitat_summary$prop_cams / exp(kappa[i + 1, ]) / 
      (cam_design$cam_A * study_design$t_steps)
    tot_u[i + 1] <- sum(u)
    
    # Update tuning parms
    if (i %% tune_check == 0) {
      batch_n <- batch_n + 1
      delta_n <- batch_n^-1
      # delta_n <- min(0.01,batch_n^-1)
      accept_gamma_check <- colMeans(accept[(i - tune_check + 1):i, 1:num_covariates], na.rm = T)
      #mean(accept[(i - tune_check + 1):i], na.rm = T)
      gamma_tune[accept_gamma_check > 0.44] <- gamma_tune[accept_gamma_check > 0.44] + delta_n
      gamma_tune[accept_gamma_check <= 0.44] <- gamma_tune[accept_gamma_check <= 0.44] - delta_n
      
      accept_kappa_check <- mean(accept[(i - tune_check + 1):i, num_covariates + 1], na.rm = T)
      kappa_tune[accept_kappa_check > 0.44] <- kappa_tune[accept_kappa_check > 0.44] + delta_n
      kappa_tune[accept_kappa_check <= 0.44] <- kappa_tune[accept_kappa_check <= 0.44] - delta_n
      
      # accept_kappa_check <- colMeans(
      #   array(
      #     accept[
      #       (i - tune_check + 1):i,
      #       (num_covariates + 1):(2 * num_covariates)
      #     ],
      #     dim = c(tune_check, num_covariates)
      #   ),
      #   na.rm = T
      # )
      # kappa_tune[which(accept_kappa_check > 0.44)] <- kappa_tune[which(accept_kappa_check > 0.44)] + delta_n
      # kappa_tune[which(accept_kappa_check <= 0.44)] <- kappa_tune[which(accept_kappa_check <= 0.44)] - delta_n
    }
  }
  # print("MCMC complete")
  
  list(accept = accept, gamma = gamma, kappa = kappa, tot_u = tot_u)
}


########################################
# MCMC for REST no covariates
########################################
fit.model.mcmc.REST <- function(study_design,
                                cam_design,
                                gamma_start,
                                kappa_start,
                                gamma_prior_var = 10^4,
                                kappa_prior_var = 10^4,
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
fit.model.mcmc.REST.cov <- function(study_design,
                                    cam_design,
                                    cam_locs,
                                    gamma_start,
                                    kappa_start,
                                    gamma_prior_var = 10^4,
                                    kappa_prior_var = 10^4,
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
fit.model.mcmc.TTE <- function(study_design,
                               cam_design,
                               gamma_start,
                               kappa_start,
                               gamma_prior_var = 10^4,
                               kappa_prior_var = 10^4,
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
fit.model.mcmc.TTE.cov <- function(study_design,
                                   cam_design,
                                   cam_locs,
                                   gamma_start,
                                   kappa_start,
                                   gamma_prior_var = 10^4,
                                   kappa_prior_var = 10^4,
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
# MCMC for STE no covariates
########################################
fit.model.mcmc.STE <- function(study_design,
                               cam_design,
                               gamma_start,
                               gamma_prior_var = 10^4,
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

# Workflow for agent-based simulations for repeated simulations

library(truncnorm)
library(Rfast)
library(raster)
library(tidyverse)
library(RColorBrewer)
library(lattice)
library(gridExtra)
library(AHMbook)
library(doParallel)

source("./R/utils.R")
source("./R/plot_funs.R")
source("./R/ABM_sim.R")
source("./R/MLE_functions.R")
source("./R/MCMC_functions.R")
source("./R/Create_landscape.R")

################################################################################
# Initializations
################################################################################

fig_colors <- c("#2ca25f", "#fc8d59", "#67a9cf", "#f768a1", "#bae4b3", "#fed98e")

# Simulation variations
sim_num <- 1

sim_vars <- data.frame(
  sim_names = c("Random", 
                "Slow_cam_bias", 
                "Medium_cam_bias", 
                "Fast_cam_bias", 
                "Slow_cam_all", 
                "Medium_cam_all", 
                "Fast_cam_all"),
  lscape_tag = "Random",
  all_speed = 1
)

lscape_tag <- sim_vars$lscape_tag[sim_num]
all_speed <- sim_vars$all_speed[sim_num] # 1: Slow, 2: Medium, 3: Fast

# Cam sample designs (1: random, 2-4: biased slow, medium, fast)
cam.dist.prop.all <- tibble::tibble(
  set_ID = c("Random", 
             "Slow_cam_bias", 
             "Medium_cam_bias", 
             "Fast_cam_bias", 
             "Slow_cam_all", 
             "Medium_cam_all", 
             "Fast_cam_all"),
  cam_prop = c(
    list(NULL),
    list(c(0.8, 0.1, 0.1)),
    list(c(0.1, 0.8, 0.1)),
    list(c(0.1, 0.1, 0.8)),
    list(c(1, 0, 0)),
    list(c(0, 1, 0)),
    list(c(0, 0, 1))
  )
)

cam.dist.prop <- cam.dist.prop.all %>% 
  dplyr::filter(set_ID == sim_vars$sim_names[sim_num])

num_runs <- 50

# Define number of clumps
num.clumps <- 100

# Define clump sizes for every clump
clump_sizes <- rep(1,num.clumps)
nind <- sum(clump_sizes)

# Landscape parms
q <- 30^2             # Number grid cells
bounds <- c(0, q^0.5) # Sampling area boundaries
t.steps <- 500        # Number of time steps
dt <- 1               # Time step size

# Grid cell lengths
dx <- (bounds[2]-bounds[1])/q^0.5
dy <- (bounds[2]-bounds[1])/q^0.5

# TTE and staying time censors
JJ <- 20  # Occasion length
t.censor <- JJ
num.occ <- t.steps/JJ

# Random walk parms
# 0 is uncorrelated random walk, inf is ideal gas model (5 is good correlation)
default_kappa <- 0
clump.rad <- dx/2 # Tightness of clumping behavior

# Camera specs
ncam <- 250
cam_length <- dx*.3 # length of all viewshed sides
cam.A <- cam_length ^ 2 / 2 # Assumes equilateral triangle viewsheds

# MCMC parms
# n.iter <- 40000
# burn.in <- 30000
n.iter <- 5000
burn.in <- 3000

################################################################################
# Define movement speeds for each cell across landscape
################################################################################

animalxy.list <- list()
tri_cam_list <- list()
num_models <- 5
# D.all <- data.frame(Model = NA,
#                     Covariate = NA,
#                     Est = NA,
#                     SD = NA,
#                     Prop_speeds = NA
# )
D.all <- tibble::tibble(
  Model = rep(NA, num_runs * num_models),
  Covariate = rep(NA, num_runs * num_models),
  Est = rep(NA, num_runs * num_models),
  SD = rep(NA, num_runs * num_models),
  Prop_speeds = rep(NA, num_runs * num_models)
)

# Repeat Simulations
for (run in 1:num_runs) {
  print(paste("Run", run, "of", num_runs))
  
  # # Create custom landscape (tag = "Random", "Homogeneous")
  lscape_speeds <- lscape_creator(tag = lscape_tag, speed_index = all_speed)
  
  # Calculate total area using available space
  tot.A <- (bounds[2]-bounds[1])^2 * sum(!is.na(lscape_speeds$Road))/q
  
  # indices for all covariates including intercept in spatial.covariates
  # Assumes 3 covariate types (slow, medium, fast)
  covariates.index <- c(0,rep(1,3))
  covariate_labels <- c("Slow", "Medium", "Fast")
  
  # Define landscape covariates
  slow_inds <- which(lscape_speeds$Speed == "Slow")
  med_inds <- which(lscape_speeds$Speed == "Medium")
  fast_inds <- which(lscape_speeds$Speed == "Fast")
  # Covariates are either 0 or 1 for indices that represent movement speeds
  Z <- matrix(0, nrow = q, ncol = length(covariate_labels))
  Z[slow_inds, 1] <- 1
  Z[med_inds, 2] <- 1
  Z[fast_inds, 3] <- 1
  
  # # Run agent-based model
  animalxy.all <- ABM_sim(bounds = bounds,
                          t.steps = t.steps,
                          speeds = matrix(lscape_speeds$Value, q^0.5, q^0.5),
                          direction = matrix(lscape_speeds$Direction, q^0.5, q^0.5),
                          kappa = matrix(lscape_speeds$Kappa, q^0.5, q^0.5),
                          road = matrix(lscape_speeds$Road, q^0.5, q^0.5),
                          clump_sizes = clump_sizes,
                          clump.rad = clump.rad,
                          dx = dx,
                          dy = dy,
                          q = q,
                          dt = dt)
  
  animalxy.list[[run]] <- animalxy.all
  
  # Create camera sample designs
  cam.samps <- sample_speeds(cam.dist.prop)
  
  # Create triangle camera viewsheds
  tri_cam_samps <- data.frame(cam_ID = 1:ncam,
                              lscape_index = cam.samps,
                              x = ((cam.samps-1) %% q^0.5 + 1)*dx - dx,
                              y = ceiling(cam.samps/q^0.5)*dx - dx) |>
    dplyr::group_by(cam_ID) |>
    dplyr::summarise(lscape_index = lscape_index,
                     x = runif(1, x + cam_length*0.5 + 0.01*dx,
                               x + dx - cam_length*0.5 - 0.01*dx) + c(0, -0.5*cam_length, 0.5*cam_length),
                     y = runif(1, y + .01*dy, y + dx - cam_length - 0.01*dy) + c(0, cam_length, cam_length),
                     vertex = c(1, 2, 3),
                     cam_area = calc_tri_area(x, y),
                     .groups = 'drop'
    )
  
  tri_cam_list[[run]] <- tri_cam_samps
  
  ################################
  # Collect data
  ################################
  '%notin%' <- Negate('%in%')
  source("./R/Collect_tele_data.R")
  source("./R/Collect_data.R")
  
  # # Quick calculation of stay time priors
  # stay_time_summary_log <- stay_time_raw_tele %>%
  #   dplyr::group_by(speed) %>%
  #   dplyr::summarise(mu = mean(log(t_stay * cam.A / dx / dy)),
  #                    sd = sd(log(t_stay * cam.A / dx / dy))) %>%
  #   dplyr::mutate(speed = factor(speed, levels = c("Slow", "Medium", "Fast")),
  #                 mu = ifelse(is.na(mu), 1, mu),
  #                 sd = ifelse(is.na(sd), 0, sd)) %>%
  #   dplyr::arrange(speed)
  # 
  # kappa.prior.mu.tdst <- stay_time_summary_log %>%
  #   dplyr::pull(mu)
  # kappa.prior.var.tdst <- stay_time_summary_log %>%
  #   dplyr::pull(sd)
  
  # # Quick math check for full cam bias, no density bias estimate
  # activity_proportion <- exp(kappa.prior.mu.tdst) / sum(exp(kappa.prior.mu.tdst))
  
  ################################################################################
  # Do it without log transform
  stay_time_summary <- stay_time_raw_tele %>% 
    dplyr::group_by(speed) %>% 
    dplyr::summarise(
      cam_mu = mean(t_stay),
      cam_sd = sd(t_stay),
      mu = mean(t_stay / cam.A * dx * dy),
      sd = sd(t_stay / cam.A * dx * dy),
      .groups = 'drop'
    ) %>% 
    dplyr::mutate(stay_prop = mu / sum(mu)) %>% 
    dplyr::rename(Speed = speed) %>% 
    dplyr::arrange(desc(Speed))

  kappa.prior.mu.tdst <- log(stay_time_summary$cam_mu)
  kappa.prior.var.tdst <- log(stay_time_summary$cam_sd)
  
  habitat_summary <- lscape_speeds %>% 
    dplyr::group_by(Speed) %>% 
    dplyr::summarise(
      n_lscape = dplyr::n()
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(
      prop_lscape = n_lscape / sum(n_lscape)
    ) %>% 
    dplyr::left_join(
      count_data %>% 
        dplyr::group_by(speed) %>% 
        dplyr::summarise(
          mean_count = mean(count),
          ncams = dplyr::n(),
          .groups = 'drop'
        ) %>% 
        dplyr::rename(Speed = speed),
      by = dplyr::join_by(Speed)
    ) %>%
    dplyr::left_join(
      stay_time_summary,
      by = dplyr::join_by(Speed)
    )  %>% 
    dplyr::mutate(
      prop_cams = ncams / sum(ncams, na.rm = T),
      stay_prop_full = sum(stay_prop) / stay_prop,
      d_coeff = n_lscape * prop_cams * stay_prop_full / (cam.A * t.steps)
    ) %>% 
    dplyr::arrange(desc(Speed))
  
  # n_adj <- habitat_summary %>%
  #   dplyr::mutate(
  #     mean_count = tidyr::replace_na(mean_count, 0),
  #     ncams = tidyr::replace_na(ncams, 0)
  #   ) %>%
  #   dplyr::left_join(
  #     stay_time_summary,
  #     by = dplyr::join_by(Speed)
  #   ) %>%
  #   dplyr::ungroup() %>%
  #   dplyr::mutate(
  #     stay_prop_full = sum(stay_prop) / stay_prop,
  #     stay_prop_adj = stay_prop_full * ncams / sum(ncams, na.rm = T),
  #     n_ltype = mean_count * n_lscape / (cam.A * t.steps),
  #     n_full = n_ltype * stay_prop_full,
  #     n_adj = n_full * prop_cams
  #     # big_D = mean_count * sum(stay_prop) / stay_prop
  #   )
  # 
  # # The sum over the adjustments should equal total abundance when cameras are placed in each lscape type
  # # Does not work if cameras aren't placed in every landscape type
  # sum(n_adj$n_ltype, na.rm = T)
  # 
  # # This method is more reliable when one or more lscape type is missing
  # sum(n_adj$n_adj)
  # 
  # # For each adjusted n, we can get individual total abundances with full weights on a lscape type
  # n_adj$n_full
  # 
  # # Manual way to calculate the above
  # n_fast_lscape <- sum(lscape_speeds$Speed == "Fast")
  # fast_cam_n <- mean(count_data$count[count_data$speed == "Fast"]) * 
  #   n_fast_lscape / (cam.A * t.steps * stay_time_summary$mu_prop[stay_time_summary$Speed == "Fast"])
  # 
  # n_slow_lscape <- sum(lscape_speeds$Speed == "Slow")
  # slow_cam_n <- mean(count_data$count[count_data$speed == "Slow"]) * n_slow_lscape / 
  #   (cam.A * t.steps * stay_time_summary$mu_prop[stay_time_summary$Speed == "Slow"])
  # adjusted_tot_n_slow <- stay_time_summary$mu_prop[stay_time_summary$Speed == "Slow"] *
  #   slow_cam_n
  
  
  # Run models only if any data points were collected
  if(max(count_data$count) > 0 & sum(!is.na(TTE_data)) > 0) {
    
    ################################
    # MCMC methods
    ################################
    # Run models in parallel
    # n_cores <- parallel::detectCores() $ Check number of cores available
    my_cluster <- parallel::makeCluster(3, type = "PSOCK")
    #register cluster to be used by %dopar%
    doParallel::registerDoParallel(cl = my_cluster)
    
    D.chain <- foreach::foreach(iter=1:num_models) %dopar% {
      ###################################
      # TDST
      ###################################
      if (iter == 1) {
        # print("Fit TDST with MCMC")
        ptm <- proc.time()
        # unpack tidyr if extract has no applicable method
        # .rs.unloadPackage("tidyr")
        chain.TDST <- fit.model.mcmc.TDST.cov(
          n.iter = n.iter,
          gamma.start = log(mean(count_data$count)),
          kappa.start = rep(log(mean(stay_time_data,na.rm=T)), 3),
          gamma.prior.var = 10^6,
          kappa.prior.var = 10^6,
          gamma.tune = -1,
          kappa.tune = c(-1, -1, -1),
          cam.counts = count_data$count,
          t.staying.dat = stay_time_data,
          covariate_labels = covariate_labels,
          covariates.index = covariates.index,
          t.steps = t.steps,
          cam.A = cam.A,
          cell.A = dx*dy,
          censor = t.censor)
        proc.time() - ptm
        
        # ## Posterior summaries
        pop.ind.TDST <- which(names(chain.TDST) == "u")
        MCMC.parms.TDST.cov <- mcmcr::as.mcmc(do.call(cbind, chain.TDST[-pop.ind.TDST])[-c(1:burn.in), ])
        summary(MCMC.parms.TDST.cov)
        
        # plot(chain.TDST$tot.u[burn.in:n.iter])
        D.TDST.MCMC <- mean(chain.TDST$tot.u[burn.in:n.iter])
        SD.TDST.MCMC <- sd(chain.TDST$tot.u[burn.in:n.iter])
        
        Prop_speeds <- c(mean(chain.TDST$u[slow_inds]),
                         mean(chain.TDST$u[med_inds]),
                         mean(chain.TDST$u[fast_inds]))/
          mean(mean(chain.TDST$u[slow_inds])+
                 mean(chain.TDST$u[med_inds])+
                 mean(chain.TDST$u[fast_inds]))
        
        if(any(colMeans(chain.TDST$accept[burn.in:n.iter,])< 0.2) ||
           any(colMeans(chain.TDST$accept[burn.in:n.iter,])> 0.7)){
          warning(('TDST accept rate OOB'))
          D.TDST.MCMC <- NA
          SD.TDST.MCMC <- NA
        }
        
        D.chain <- tibble::tibble(
          Model = "TDST",
          Covariate = "Covariate",
          Est = D.TDST.MCMC,
          SD = SD.TDST.MCMC,
          Prop_speeds = list(Prop_speeds)
        )
      }
      
      ###################################
      # TDST w/ prior
      ###################################
      if (iter == 5) {
        ptm <- proc.time()
        # unpack tidyr if extract has no applicable method
        # .rs.unloadPackage("tidyr")
        chain.TDST.prior <- fit.model.mcmc.TDST.cov(
          n.iter = n.iter,
          gamma.start = log(mean(count_data$count)),
          kappa.start = kappa.prior.mu.tdst,#rep(log(mean(stay_time_data,na.rm=T)), 3),
          gamma.prior.var = 10^6,
          kappa.prior.mu = kappa.prior.mu.tdst,
          kappa.prior.var = exp(kappa.prior.var.tdst)/1000,
          gamma.tune = -1,
          kappa.tune = c(-1, -1, -1),
          cam.counts = count_data$count,
          t.staying.dat = NULL, #stay_time_data,
          covariate_labels = covariate_labels,
          covariates.index = covariates.index,
          t.steps = t.steps,
          cam.A = cam.A,
          cell.A = dx*dy,
          censor = t.censor)
        proc.time() - ptm

        # ## Posterior summaries
        pop.ind.TDST.prior <- which(names(chain.TDST.prior) == "u")
        MCMC.parms.TDST.cov.prior <- mcmcr::as.mcmc(do.call(cbind, chain.TDST.prior[-pop.ind.TDST.prior])[-c(1:burn.in), ])
        summary(MCMC.parms.TDST.cov.prior)
        
        # plot(chain.TDST.prior$tot.u[burn.in:n.iter])
        D.TDST.MCMC.prior <- mean(chain.TDST.prior$tot.u[burn.in:n.iter])
        SD.TDST.MCMC.prior <- sd(chain.TDST.prior$tot.u[burn.in:n.iter])

        Prop_speeds <- c(mean(chain.TDST.prior$u[slow_inds]),
                         mean(chain.TDST.prior$u[med_inds]),
                         mean(chain.TDST.prior$u[fast_inds]))/
          mean(mean(chain.TDST.prior$u[slow_inds])+
                 mean(chain.TDST.prior$u[med_inds])+
                 mean(chain.TDST.prior$u[fast_inds]))

        if(any(colMeans(chain.TDST.prior$accept[burn.in:n.iter,])< 0.2) ||
           any(colMeans(chain.TDST.prior$accept[burn.in:n.iter,])> 0.7)){
          warning(('TDST accept rate OOB'))
          D.TDST.MCMC.prior <- NA
          SD.TDST.MCMC.prior <- NA
        }

        D.chain <- tibble::tibble(
          Model = "TDST priors",
          Covariate = "Covariate",
          Est = D.TDST.MCMC.prior,
          SD = SD.TDST.MCMC.prior,
          Prop_speeds = list(Prop_speeds)
        )
      }
      
      ########################################
      ## PR no covariates
      ########################################
      if (iter == 2) {
        # print("Fit Poisson Regression model with MCMC, no covariates")
        ptm <- proc.time()
        # unpack tidyr if extract has no applicable method
        # .rs.unloadPackage("tidyr")
        chain.PR <-fit.model.mcmc.PR(
          n.iter = n.iter,
          gamma.start = log(mean(count_data$count)),
          gamma.prior.var = 10^6,
          gamma.tune = -1,
          cam.counts = count_data$count,
          sample_frame = tot.A,
          cam.A,
          t.steps)
        proc.time() - ptm
        
        # ## Posterior summaries
        MCMC.parms.PR <- mcmcr::as.mcmc(do.call(cbind, chain.PR)[-c(1:burn.in), ])
        summary(MCMC.parms.PR)
        
        # plot(chain.PR$tot.u[burn.in:n.iter])
        D.PR.MCMC <- mean(chain.PR$tot.u[burn.in:n.iter])
        SD.PR.MCMC <- sd(chain.PR$tot.u[burn.in:n.iter])
        
        if(mean(chain.PR$accept[burn.in:n.iter,])< 0.2 || mean(chain.PR$accept[burn.in:n.iter,])> 0.7){
          warning(('Poisson Regression accept rate OOB'))
          D.PR.MCMC <- NA
          SD.PR.MCMC <- NA
        }
        
        Prop_speeds <- c(mean(chain.PR$u[slow_inds]),
                         mean(chain.PR$u[med_inds]),
                         mean(chain.PR$u[fast_inds]))/
          mean(mean(chain.PR$u[slow_inds])+
                 mean(chain.PR$u[med_inds])+
                 mean(chain.PR$u[fast_inds]))
        
        D.chain <- tibble::tibble(
          Model = "PR",
          Covariate = "Non-Covariate",
          Est = D.PR.MCMC,
          SD = SD.PR.MCMC,
          Prop_speeds = list(Prop_speeds)
        )
      }
      
      ########################################
      ## PR habitat strat
      ########################################
      if (iter == 3) {
        # print("Fit Poisson Regression model with MCMC, no covariates")
        ptm <- proc.time()

        # .rs.unloadPackage("tidyr")
        chain.PR.habitat <- fit.model.mcmc.PR.habitat(
          n.iter = n.iter,
          gamma.start = rep(log(mean(count_data$count)), 3),
          gamma.prior.var = 10^6,
          gamma.tune = c(-1, -1, -1),
          cam.counts = count_data,
          cam.A,
          t.steps,
          habitat_summary = habitat_summary)
        proc.time() - ptm
        
        # ## Posterior summaries
        MCMC.parms.PR.habitat <- mcmcr::as.mcmc(do.call(cbind, chain.PR.habitat)[-c(1:burn.in), ])
        summary(MCMC.parms.PR.habitat)
        
        # plot(chain.PR.habitat$tot.u[burn.in:n.iter])
        D.PR.habitat.MCMC <- mean(chain.PR.habitat$tot.u[burn.in:n.iter])
        SD.PR.habitat.MCMC <- sd(chain.PR.habitat$tot.u[burn.in:n.iter])
        

        if(mean(chain.PR.habitat$accept[burn.in:n.iter,])< 0.2 || mean(chain.PR.habitat$accept[burn.in:n.iter,])> 0.7){
          warning(('Poisson Regression accept rate OOB'))
          D.PR.habitat.MCMC <- NA
          SD.PR.habitat.MCMC <- NA
        }

        Prop_speeds <- c(mean(chain.PR.habitat$u[slow_inds]),
                         mean(chain.PR.habitat$u[med_inds]),
                         mean(chain.PR.habitat$u[fast_inds]))/
          mean(mean(chain.PR.habitat$u[slow_inds])+
                 mean(chain.PR.habitat$u[med_inds])+
                 mean(chain.PR.habitat$u[fast_inds]))
        
        D.chain <- tibble::tibble(
          Model = "PR Habitat",
          Covariate = "Non-Covariate",
          Est = D.PR.habitat.MCMC,
          SD = SD.PR.habitat.MCMC,
          Prop_speeds = list(Prop_speeds)
        )
      }
      
      ########################################
      ## PR w/ covariates
      ########################################
      if (iter == 4) {
        # print("Fit Poisson Regression model with MCMC w/ covariates")
        ptm <- proc.time()
        # unpack tidyr if extract has no applicable method
        # .rs.unloadPackage("tidyr")
        chain.PR.cov <-fit.model.mcmc.PR.cov(
          n.iter = n.iter,
          gamma.start = rep(log(mean(count_data$count)), 3),
          gamma.prior.var = 10^6,
          gamma.tune = c(-1, -1, -1),
          cam.counts = count_data$count,
          covariate_labels = covariate_labels,
          covariates.index = covariates.index,
          cam.A = cam.A,
          cell.A = dx*dy)
        proc.time() - ptm
        
        # ## Posterior summaries
        pop.ind.PR <- which(names(chain.PR.cov) == "u")
        MCMC.parms.PR.cov <- mcmcr::as.mcmc(do.call(cbind, chain.PR.cov[-pop.ind.PR])[-c(1:burn.in), ])
        summary(MCMC.parms.PR.cov)
        
        # plot(chain.PR.cov$tot.u[burn.in:n.iter])
        D.PR.MCMC.cov <- mean(chain.PR.cov$tot.u[burn.in:n.iter])
        SD.PR.MCMC.cov <- sd(chain.PR.cov$tot.u[burn.in:n.iter])
        
        Prop_speeds <- c(mean(chain.PR.cov$u[slow_inds]),
                         mean(chain.PR.cov$u[med_inds]),
                         mean(chain.PR.cov$u[fast_inds]))/
          mean(mean(chain.PR.cov$u[slow_inds])+
                 mean(chain.PR.cov$u[med_inds])+
                 mean(chain.PR.cov$u[fast_inds]))
        
        if(mean(chain.PR.cov$accept[burn.in:n.iter,])< 0.2 || mean(chain.PR.cov$accept[burn.in:n.iter,])> 0.7){
          warning(('Poisson Regression accept rate OOB'))
          D.PR.MCMC.cov <- NA
          SD.PR.MCMC.cov <- NA
        }
        
        D.chain <- tibble::tibble(
          Model = "PR",
          Covariate = "Covariate",
          Est = D.PR.MCMC.cov,
          SD = SD.PR.MCMC.cov,
          Prop_speeds = list(Prop_speeds)
        )
      }
      
      D.chain <- D.chain
    }
    stopCluster(my_cluster)
    
  }
  
  D.all[(num_models * (run - 1) + 1):(num_models * run), ] <- dplyr::bind_rows(D.chain)
}


ggplot(D.all, aes(x = Model, y = Est, fill = Covariate)) +
  geom_boxplot(lwd = .1, fatten = .1) +
  # geom_boxplot(data=subset(D.all, D.all$Model %in% c("REST", "TTE", "PR")), colour = "black") +
  # geom_boxplot(data=subset(D.all, D.all$Model == "TDST"), colour=c("black","white")) +
  # geom_boxplot(data=subset(D.all, D.all$Model == "STE"), colour=c("white","black")) +
  geom_hline(yintercept=100, linetype="dashed", size=1) +
  labs(x = "Model",
       # y = paste(cam.props.label[lscape_var+1], " Landscape \n Mean Abundance")) +
       y = "Mean Abundance") +
  scale_fill_manual(values=c('grey40','Grey')) +
  theme(text = element_text(size = 20),
        legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = c(0.17, 0.84),
        legend.background = element_blank(),
        legend.spacing.y = unit(0, "mm"),
        legend.box.background = element_rect(colour = "black"))

ggplot(D.all, aes(x = Model, y = SD, fill = Covariate)) +
  geom_boxplot(lwd = .1, fatten = .1) +
  # geom_boxplot(data=subset(D.all, D.all$Model %in% c("REST", "TTE", "PR")), colour = "black") +
  # geom_boxplot(data=subset(D.all, D.all$Model == "TDST"), colour=c("black","white")) +
  # geom_boxplot(data=subset(D.all, D.all$Model == "STE"), colour=c("white","black")) +
  labs(x = "Model",
       # y = paste(cam.props.label[lscape_var+1], " Landscape \n Mean Abundance")) +
       y = "SD") +
  scale_fill_manual(values=c('grey40','Grey')) +
  theme(text = element_text(size = 20),
        legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = c(0.17, 0.84),
        legend.background = element_blank(),
        legend.spacing.y = unit(0, "mm"),
        legend.box.background = element_rect(colour = "black"))


# D.all <- D.all[-1,]

# # Remove outlier estimates
# D.all$Est[D.all$Est > 5*nind] <- NA
# D.all$SD[D.all$SD > 5*nind] <- NA
# 
# D.all$Model <- factor(D.all$Model, levels = c("TDST", "TDST priors", "REST", "TTE", "PR", "STE"))
# 
# D_summary <- D.all |>
#   dplyr::group_by(Model) |>
#   dplyr::summarise(num_NAs = sum(is.na(Est)))
# 
# 
# ####################################
# # Plot the stuff
# ####################################
# cam.props.label <- c("Camera Bias: Random",
#                      "Camera Bias: Slow",
#                      "Camera Bias: Medium",
#                      "Camera Bias: Fast")
# 
# # if(num_runs == 1) {
# #   plot_onerun_results()
# # } else{
# #   plot_multirun_means()
# #   # plot_multirun_sds()
# #   # plot_multirun_hist()
# # }
# # 
# # plot_count_data(fill = "speed")
# # plot_encounter_data(fill = "speed")
# # plot_staytime_data(fill = "speed")
# # plot_TTE_data(fill = "speed")
# 
# # # Plot ABM simulations
# # plot_ABM()
# # plot_space_use()
# # plot_ABM_stay_proportions()
# # plot_model_proportions()
# 
# # # Save Results
# # save(D.all, animalxy.list, tri_cam_list, file = paste0("Sim_results/Sim_", sim_vars$sim_names[sim_num], "_all_vars.RData"))


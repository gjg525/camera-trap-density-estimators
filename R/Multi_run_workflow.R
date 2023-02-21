# Workflow for agent-based simulations for repeated simulations

# TODO: CREATE one data variable for fitting models
# Adapt covariate functions for tri viewsheds

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
  sim_names = c("Original", "Slow_landscape", "Medium_landscape", "Fast_landscape", "Slow_cams", "Medium_cams", "Fast_cams"),
  lscape_tag = c("Random", rep("Homogeneous", 3), rep("Random", 3)),
  all_speed = c(1, 1, 2, 3, rep(1, 3)),
  cam.dist.set = c(rep(1, 4), 2, 3, 4)
)

lscape_tag <- sim_vars$lscape_tag[sim_num]
all_speed <- sim_vars$all_speed[sim_num] # 1: Slow, 2: Medium, 3: Fast

# Cam sample designs (1: random, 2-4: biased slow, medium, fast)
cam.dist.set <- sim_vars$cam.dist.set[sim_num]

# D.all <- data.frame(Model = character(),
#                     Covariate = character(),
#                     Est = double(),
#                     SD = double(),
#                     Prop_speeds = double()
# )
D.all <- data.frame(Model = NA,
                    Covariate = NA,
                    Est = NA,
                    SD = NA,
                    Prop_speeds = NA
)
num_runs <- 2

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
# cam_length <- dx*.7 # length of all viewshed sides
cam.A <- cam_length ^ 2 / 2

# MCMC parms
n.iter <- 10000
burn.in <- 5000

################################################################################
# Define movement speeds for each cell across landscape
################################################################################

# Repeat Simulations
for (run in 1:num_runs) {
  
  # # Create custom landscape (tag = "Random", "Homogeneous")
  lscape_speeds <- lscape_creator(tag = lscape_tag, speed_index = all_speed)
  
  # Calculate total area using available space
  tot.A <- (bounds[2]-bounds[1])^2 * sum(!is.na(lscape_speeds$Road))/q
  
  # indices for all covariates including intercept in spatial.covariates
  # Assumes 3 covariate types (slow, medium, fast)
  covariates.index <- c(0,rep(1,3))
  covariate_labels <- c("Slow", "Medium", "Fast")
  
  # Define covariate speeds (There is a better way to do this)
  Z <- matrix(0, nrow = q, ncol = length(covariate_labels))
  Z[lscape_speeds$Index[lscape_speeds$Speed == "Slow"], 1] <- lscape_speeds$Value[lscape_speeds$Speed == "Slow"]
  Z[lscape_speeds$Index[lscape_speeds$Speed == "Medium"], 2] <- lscape_speeds$Value[lscape_speeds$Speed == "Medium"]
  Z[lscape_speeds$Index[lscape_speeds$Speed == "Fast"], 3] <- lscape_speeds$Value[lscape_speeds$Speed == "Fast"]

  slow_inds <- which(lscape_speeds$Speed == "Slow")
  med_inds <- which(lscape_speeds$Speed == "Medium")
  fast_inds <- which(lscape_speeds$Speed == "Fast")
  
  print(paste("Run", run, "of", num_runs))

  # Create camera sample designs
  if(cam.dist.set == 1) {
    # # Random camera sampling
    cam.samps <- sample(which(!is.na(lscape_speeds$Road)), ncam, replace = F)
  } else {
    # Biased sample design (lscape_speeds)
    cam.samps <- sample_speeds(cam.dist.set)

  }

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

  # # Run agent-based model
  animalxy.all <- ABM_sim(bounds,
                          t.steps,
                          speeds = matrix(lscape_speeds$Value, q^0.5, q^0.5),
                          direction = matrix(lscape_speeds$Direction, q^0.5, q^0.5),
                          kappa = matrix(lscape_speeds$Kappa, q^0.5, q^0.5),
                          road = matrix(lscape_speeds$Road, q^0.5, q^0.5),
                          clump_sizes,
                          clump.rad)

  ################################
  # Collect data
  ################################
  '%notin%' <- Negate('%in%')
  source("./R/Collect_data.R")

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

    D.chain <- foreach::foreach(iter=1:8) %dopar% {
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
    # pop.ind.TDST <- which(names(chain.TDST) == "u")
    # MCMC.parms.TDST.cov <- mcmcr::as.mcmc(do.call(cbind, chain.TDST[-pop.ind.TDST])[-c(1:burn.in), ])
    # summary(MCMC.parms.TDST.cov)

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

    # D.all<- rbind(D.all, data.frame(
    #   Model = "TDST",
    #   Covariate = "Covariate",
    #   Est = D.TDST.MCMC,
    #   SD = SD.TDST.MCMC))

    D.chain <- tibble::tibble(
      Model = "TDST",
      Covariate = "Covariate",
      Est = D.TDST.MCMC,
      SD = SD.TDST.MCMC,
      Prop_speeds = list(Prop_speeds)
    )
    }

    ###################################
    # REST no covariates
    ###################################
    if (iter == 2) {
    # print("Fit REST with MCMC, no covariates")
    ptm <- proc.time()
    # unpack tidyr if extract has no applicable method
    # .rs.unloadPackage("tidyr")
    chain.REST <-fit.model.mcmc.REST(
      n.iter = n.iter,
      gamma.start = log(mean(encounter_data$encounter)),
      kappa.start = log(mean(stay_time_data,na.rm=T)),
      gamma.prior.var = 10^6,
      kappa.prior.var = 10^6,
      gamma.tune = c(-1, -1, -1),
      kappa.tune = c(-1, -1, -1),
      num.encounters.dat = encounter_data$encounter,
      t.staying.dat = stay_time_data,
      t.steps = t.steps,
      cam.A = cam.A,
      censor = t.censor)
    proc.time() - ptm

    # ## Posterior summaries
    MCMC.parms.REST <- mcmcr::as.mcmc(do.call(cbind, chain.REST)[-c(1:burn.in), ])
    summary(MCMC.parms.REST)

    # plot(chain.REST$tot.u[burn.in:n.iter])
    D.REST.MCMC <- mean(chain.REST$tot.u[burn.in:n.iter])
    SD.REST.MCMC <- sd(chain.REST$tot.u[burn.in:n.iter])

    if(any(colMeans(chain.REST$accept[burn.in:n.iter,])< 0.2) || any(colMeans(chain.REST$accept[burn.in:n.iter,])> 0.7)){
      warning(('REST accept rate OOB'))
      D.REST.MCMC <- NA
      SD.REST.MCMC <- NA
    }

    Prop_speeds <- c(mean(chain.REST$u[slow_inds]),
                     mean(chain.REST$u[med_inds]),
                     mean(chain.REST$u[fast_inds]))/
      mean(mean(chain.REST$u[slow_inds])+
             mean(chain.REST$u[med_inds])+
             mean(chain.REST$u[fast_inds]))
    
    # D.all<- rbind(D.all, data.frame(
    #   Model = "REST",
    #   Covariate = "Non-Covariate",
    #   Est = D.REST.MCMC,
    #   SD = SD.REST.MCMC))

    D.chain <- tibble::tibble(
      Model = "REST",
      Covariate = "Non-Covariate",
      Est = D.REST.MCMC,
      SD = SD.REST.MCMC,
      Prop_speeds = list(Prop_speeds)
    )
    }

    ###################################
    # REST w/ covariates
    ###################################
    if (iter == 3) {
    # print("Fit REST with MCMC")
    ptm <- proc.time()
    # unpack tidyr if extract has no applicable method
    # .rs.unloadPackage("tidyr")
    chain.REST.cov <- fit.model.mcmc.REST.cov(
      n.iter = n.iter,
      gamma.start = rep(log(mean(encounter_data$encounter)), 3),
      kappa.start = rep(log(mean(stay_time_data,na.rm=T)), 3),
      gamma.prior.var = 10^6,
      kappa.prior.var = 10^6,
      gamma.tune = c(-1, -1, -1),
      kappa.tune = c(-1, -1, -1),
      num.encounters.dat = encounter_data$encounter,
      t.staying.dat = stay_time_data,
      covariate_labels = covariate_labels,
      covariates.index = covariates.index,
      t.steps = t.steps,
      cam.A = cam.A,
      cell.A = dx*dy,
      censor = t.censor)
    proc.time() - ptm

    # ## Posterior summaries
    pop.ind.REST <- which(names(chain.REST.cov) == "u")
    MCMC.parms.REST.cov <- mcmcr::as.mcmc(do.call(cbind, chain.REST.cov[-pop.ind.REST])[-c(1:burn.in), ])
    summary(MCMC.parms.REST.cov)

    # plot(chain.REST$tot.u[burn.in:n.iter])
    D.REST.MCMC.cov <- mean(chain.REST.cov$tot.u[burn.in:n.iter])
    SD.REST.MCMC.cov <- sd(chain.REST.cov$tot.u[burn.in:n.iter])

    if(any(colMeans(chain.REST.cov$accept[burn.in:n.iter,])< 0.2) ||
       any(colMeans(chain.REST.cov$accept[burn.in:n.iter,])> 0.7)){
      warning(('REST accept rate OOB'))
      D.REST.MCMC.cov <- NA
      SD.REST.MCMC.cov <- NA
    }

    Prop_speeds <- c(mean(chain.REST.cov$u[slow_inds]),
                     mean(chain.REST.cov$u[med_inds]),
                     mean(chain.REST.cov$u[fast_inds]))/
      mean(mean(chain.REST.cov$u[slow_inds])+
           mean(chain.REST.cov$u[med_inds])+
           mean(chain.REST.cov$u[fast_inds]))
    
    # D.all<- rbind(D.all, data.frame(
    #   Model = "REST",
    #   Covariate = "Covariate",
    #   Est = D.REST.MCMC.cov,
    #   SD = SD.REST.MCMC.cov))

    D.chain <- tibble::tibble(
      Model = "REST",
      Covariate = "Covariate",
      Est = D.REST.MCMC.cov,
      SD = SD.REST.MCMC.cov,
      Prop_speeds = list(Prop_speeds)
    )
    }

    ########################################
    ## TTE no covariates
    ########################################
    if (iter == 4) {
    # print("Fit TTE with MCMC, no covariates")
    ptm <- proc.time()
    # unpack tidyr if extract has no applicable method
    # .rs.unloadPackage("tidyr")
    chain.TTE <-fit.model.mcmc.TTE(
      n.iter = n.iter,
      gamma.start = log(mean(TTE_data, na.rm = T)),
      kappa.start = log(mean(stay_time_data,na.rm=T)),
      gamma.prior.var = 10^6,
      kappa.prior.var = 10^6,
      gamma.tune = -1,
      kappa.tune = -1,
      TTE.dat = TTE_data,
      t.staying.dat = stay_time_data,
      censor = t.censor,
      JJ = JJ,
      cam.A = cam.A)
    proc.time() - ptm

    Prop_speeds <- c(mean(chain.TTE$u[slow_inds]),
                     mean(chain.TTE$u[med_inds]),
                     mean(chain.TTE$u[fast_inds]))/
      mean(mean(chain.TTE$u[slow_inds])+
             mean(chain.TTE$u[med_inds])+
             mean(chain.TTE$u[fast_inds]))
    
    # ## Posterior summaries
    # MCMC.parms.TTE <- mcmcr::as.mcmc(do.call(cbind, chain.TTE)[-c(1:burn.in), ])
    # summary(MCMC.parms.TTE)

    # plot(chain.TTE$tot.u[burn.in:n.iter])
    D.TTE.MCMC <- mean(chain.TTE$tot.u[burn.in:n.iter])
    SD.TTE.MCMC <- sd(chain.TTE$tot.u[burn.in:n.iter])

    if(any(colMeans(chain.TTE$accept[burn.in:n.iter,])< 0.2) || any(colMeans(chain.TTE$accept[burn.in:n.iter,])> 0.7)){
      warning(('TTE accept rate OOB'))
      D.TTE.MCMC <- NA
      SD.TTE.MCMC <- NA
    }

    # D.all<- rbind(D.all, data.frame(
    #   Model = "TTE",
    #   Covariate = "Non-Covariate",
    #   Est = D.TTE.MCMC,
    #   SD = SD.TTE.MCMC))

    D.chain <- tibble::tibble(
      Model = "TTE",
      Covariate = "Non-Covariate",
      Est = D.TTE.MCMC,
      SD = SD.TTE.MCMC,
      Prop_speeds = list(Prop_speeds)
    )

    }

    ########################################
    ## TTE w/ covariates
    ########################################
    if (iter == 5) {
    # print("Fit TTE with MCMC w/ covariates")
    ptm <- proc.time()
    # unpack tidyr if extract has no applicable method
    # .rs.unloadPackage("tidyr")
    chain.TTE.cov <-fit.model.mcmc.TTE.cov(
      n.iter = n.iter,
      gamma.start = rep(log(mean(TTE_data, na.rm = T)), 3),
      kappa.start = rep(log(mean(stay_time_data,na.rm=T)), 3),
      gamma.prior.var = 10^6,
      kappa.prior.var = 10^6,
      gamma.tune = c(-1, -1, -1),
      kappa.tune = c(-1, -1, -1),
      TTE.dat = TTE_data,
      t.staying.dat = stay_time_data,
      covariate_labels = covariate_labels,
      covariates.index = covariates.index,
      censor = t.censor,
      JJ = JJ,
      cam.A = cam.A,
      cell.A = dx*dy)
    proc.time() - ptm

    # ## Posterior summaries
    # pop.ind.TTE <- which(names(chain.TTE.cov) == "u")
    # MCMC.parms.TTE.cov <- mcmcr::as.mcmc(do.call(cbind, chain.TTE.cov[-pop.ind.TTE])[-c(1:burn.in), ])
    # summary(MCMC.parms.TTE.cov)

    # plot(chain.TTE.cov$tot.u[burn.in:n.iter])
    D.TTE.MCMC.cov <- mean(chain.TTE.cov$tot.u[burn.in:n.iter])
    SD.TTE.MCMC.cov <- sd(chain.TTE.cov$tot.u[burn.in:n.iter])

    Prop_speeds <- c(mean(chain.TTE.cov$u[slow_inds]),
                     mean(chain.TTE.cov$u[med_inds]),
                     mean(chain.TTE.cov$u[fast_inds]))/
      mean(mean(chain.TTE.cov$u[slow_inds])+
             mean(chain.TTE.cov$u[med_inds])+
             mean(chain.TTE.cov$u[fast_inds]))
    
    if(any(colMeans(chain.TTE.cov$accept[burn.in:n.iter,])< 0.2) ||
       any(colMeans(chain.TTE.cov$accept[burn.in:n.iter,])> 0.7)){
      warning(('TTE accept rate OOB'))
      D.TTE.MCMC.cov <- NA
      SD.TTE.MCMC.cov <- NA
    }

    # D.all<- rbind(D.all, data.frame(
    #   Model = "TTE",
    #   Covariate = "Covariate",
    #   Est = D.TTE.MCMC.cov,
    #   SD = SD.TTE.MCMC.cov))

    D.chain <- tibble::tibble(
      Model = "TTE",
      Covariate = "Covariate",
      Est = D.TTE.MCMC.cov,
      SD = SD.TTE.MCMC.cov,
      Prop_speeds = list(Prop_speeds)
    )
    }

    ########################################
    ## MCT no covariates
    ########################################
    if (iter == 6) {
    # print("Fit Mean Count model with MCMC, no covariates")
    ptm <- proc.time()
    # unpack tidyr if extract has no applicable method
    # .rs.unloadPackage("tidyr")
    chain.MCT <-fit.model.mcmc.MCT(
      n.iter = n.iter,
      gamma.start = log(mean(count_data$count)),
      gamma.prior.var = 10^6,
      gamma.tune = -1,
      cam.counts = count_data$count)
    proc.time() - ptm

    # ## Posterior summaries
    # MCMC.parms.MCT <- mcmcr::as.mcmc(do.call(cbind, chain.MCT)[-c(1:burn.in), ])
    # summary(MCMC.parms.MCT)

    # plot(chain.MCT$tot.u[burn.in:n.iter])
    D.MCT.MCMC <- mean(chain.MCT$tot.u[burn.in:n.iter])
    SD.MCT.MCMC <- sd(chain.MCT$tot.u[burn.in:n.iter])

    Prop_speeds <- c(mean(chain.MCT$u[slow_inds]),
                     mean(chain.MCT$u[med_inds]),
                     mean(chain.MCT$u[fast_inds]))/
      mean(mean(chain.MCT$u[slow_inds])+
             mean(chain.MCT$u[med_inds])+
             mean(chain.MCT$u[fast_inds]))
    
    if(mean(chain.MCT$accept[burn.in:n.iter,])< 0.2 || mean(chain.MCT$accept[burn.in:n.iter,])> 0.7){
      warning(('Mean Count accept rate OOB'))
      D.MCT.MCMC <- NA
      SD.MCT.MCMC <- NA
    }

    # D.all<- rbind(D.all, data.frame(
    #   Model = "MCT",
    #   Covariate = "Non-Covariate",
    #   Est = D.MCT.MCMC,
    #   SD = SD.MCT.MCMC))

    D.chain <- tibble::tibble(
      Model = "MCT",
      Covariate = "Non-Covariate",
      Est = D.MCT.MCMC,
      SD = SD.MCT.MCMC,
      Prop_speeds = list(Prop_speeds)
    )
    }

    ########################################
    ## MCT w/ covariates
    ########################################
    if (iter == 7) {
    # print("Fit Mean Count model with MCMC w/ covariates")
    ptm <- proc.time()
    # unpack tidyr if extract has no applicable method
    # .rs.unloadPackage("tidyr")
    chain.MCT.cov <-fit.model.mcmc.MCT.cov(
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
    # pop.ind.MCT <- which(names(chain.MCT.cov) == "u")
    # MCMC.parms.MCT.cov <- mcmcr::as.mcmc(do.call(cbind, chain.MCT.cov[-pop.ind.MCT])[-c(1:burn.in), ])
    # summary(MCMC.parms.MCT.cov)

    # plot(chain.MCT$tot.u[burn.in:n.iter])
    D.MCT.MCMC.cov <- mean(chain.MCT.cov$tot.u[burn.in:n.iter])
    SD.MCT.MCMC.cov <- sd(chain.MCT.cov$tot.u[burn.in:n.iter])

    Prop_speeds <- c(mean(chain.MCT.cov$u[slow_inds]),
                     mean(chain.MCT.cov$u[med_inds]),
                     mean(chain.MCT.cov$u[fast_inds]))/
      mean(mean(chain.MCT.cov$u[slow_inds])+
             mean(chain.MCT.cov$u[med_inds])+
             mean(chain.MCT.cov$u[fast_inds]))
    
    if(mean(chain.MCT.cov$accept[burn.in:n.iter,])< 0.2 || mean(chain.MCT.cov$accept[burn.in:n.iter,])> 0.7){
      warning(('Mean Count accept rate OOB'))
      D.MCT.MCMC.cov <- NA
      SD.MCT.MCMC.cov <- NA
    }

    # D.all<- rbind(D.all, data.frame(
    #   Model = "MCT",
    #   Covariate = "Covariate",
    #   Est = D.MCT.MCMC.cov,
    #   SD = SD.MCT.MCMC.cov))

    D.chain <- tibble::tibble(
      Model = "MCT",
      Covariate = "Covariate",
      Est = D.MCT.MCMC.cov,
      SD = SD.MCT.MCMC.cov,
      Prop_speeds = list(Prop_speeds)
    )
    }

  ########################################
  ## STE no covariates
  ########################################
  if (iter == 8) {
  # print("Fit STE with MCMC, no covariates")
  ptm <- proc.time()
  # unpack tidyr if extract has no applicable method
  # .rs.unloadPackage("tidyr")
  chain.STE <-fit.model.mcmc.STE(
    n.iter = n.iter,
    gamma.start = 1/mean(STE_data$STE, na.rm = T),
    gamma.prior.var = 10^6,
    gamma.tune = -1,
    STE.dat = STE_data$STE,
    censor = ncam*cam.A,
    cam.A = cam.A)
  proc.time() - ptm

  # ## Posterior summaries
  # # plot(chain.STE$tot.u[burn.in:n.iter])
  # MCMC.parms.STE <- mcmcr::as.mcmc(do.call(cbind, chain.STE)[-c(1:burn.in), ])
  # summary(MCMC.parms.STE)

  # plot(chain.STE$tot.u[burn.in:n.iter])
  D.STE.MCMC <- mean(chain.STE$tot.u[burn.in:n.iter])
  SD.STE.MCMC <- sd(chain.STE$tot.u[burn.in:n.iter])

  Prop_speeds <- c(mean(chain.STE$u[slow_inds]),
                   mean(chain.STE$u[med_inds]),
                   mean(chain.STE$u[fast_inds]))/
    mean(mean(chain.STE$u[slow_inds])+
           mean(chain.STE$u[med_inds])+
           mean(chain.STE$u[fast_inds]))
  
  if(mean(chain.STE$accept[burn.in:n.iter])< 0.2 || mean(chain.STE$accept[burn.in:n.iter,])> 0.7){
    warning(('STE accept rate OOB'))
    D.STE.MCMC <- NA
    SD.STE.MCMC <- NA
  }

  # D.all<- rbind(D.all, data.frame(
  #   Model = "STE",
  #   Covariate = "Non-Covariate",
  #   Est = D.STE.MCMC,
  #   SD = SD.STE.MCMC))

  D.chain <- tibble::tibble(
    Model = "STE",
    Covariate = "Non-Covariate",
    Est = D.STE.MCMC,
    SD = SD.STE.MCMC,
    Prop_speeds = list(Prop_speeds)
  )
  }

    D.chain <- D.chain
    }
  } # End

  # Stop cluster
  stopCluster(my_cluster)

  # Is there a better way to select each element in the list?
  for (ch in 1:8) {
    chain <- D.chain[[ch]]
    D.all <- D.all |>
      dplyr::add_row(Model = chain$Model,
                     Covariate = chain$Covariate,
                     Est = chain$Est,
                     SD = chain$SD,
                     Prop_speeds = chain$Prop_speeds
      )
  }
}

# Remove outlier estimates
D.all$Est[D.all$Est > 5*nind] <- NA
D.all$SD[D.all$SD > 5*nind] <- NA

D.all$Model <- factor(D.all$Model, levels = c("TDST", "REST", "TTE", "MCT", "STE"))

D_summary <- D.all |>
  dplyr::group_by(Model) |>
  dplyr::summarise(num_NAs = sum(is.na(Est)))

####################################
# Plot the stuff
####################################
cam.props.label <- c("Camera Bias: Random",
                     "Camera Bias: Slow",
                     "Camera Bias: Medium",
                     "Camera Bias: Fast")

if(num_runs == 1) {
  plot_onerun_results()
} else{
  plot_multirun_means()
  # plot_multirun_sds()
  # plot_multirun_hist()
}

# number of counts across whole landscape for each covariate type
u.abm.all = matrix(0, nrow = q^0.5, ncol = q^0.5)
for(xx in 1:nrow(animalxy.all)) {
  x.round <- ceiling(animalxy.all$x[xx]*(q^0.5/max(bounds)))
  y.round <- ceiling(animalxy.all$y[xx]*(q^0.5/max(bounds)))
  
  # Transpose for converting matrix to raster
  u.abm.all[x.round,y.round] <- u.abm.all[x.round,y.round]+1
}

slow.counts <- animal_data %>% filter(lscape_type == "Slow") %>% pull(t_spent)
med.counts <- animal_data %>% filter(lscape_type == "Medium") %>% pull(t_spent)
fast.counts <- animal_data %>% filter(lscape_type == "Fast") %>% pull(t_spent)

Prop_all <- D.all %>% 
  dplyr::filter(Covariate == "Covariate") %>% 
  unnest_wider(Prop_speeds, names_sep="_") %>% 
  rename(Slow = Prop_speeds_1,
         Medium = Prop_speeds_2,
         Fast = Prop_speeds_3) %>% 
  select(Model, Slow, Medium, Fast) %>% 
  pivot_longer(!Model, names_to = "Speed", values_to = "Proportions") %>% 
  group_by(Model, Speed) %>% 
  summarise(Means = mean(Proportions),
            SDs = sd(Proportions),
            .groups = 'drop') %>% 
  add_row(Model = "ABM",
          Speed = c("Slow", "Medium", "Fast"),
          Means = c(slow.counts, med.counts, fast.counts))
  

  ggplot(Prop_all, aes(x=Speed, y = Means, fill = Model)) +
        geom_bar(stat="identity", color="black", position=position_dodge()) +
        geom_errorbar(aes(ymin=Means-SDs, ymax=Means+SDs), width=.2,
                      position=position_dodge(.9), size = .7) +
        # scale_y_continuous(limits=c(0, max(Prop_all$Means) +.01), expand = c(0, 0)) +
        labs(x = "Landscape Type",
             y = "Relative Distributions") +
        scale_fill_manual(values = c("grey40",fig_colors[1:4])) +
        theme(text = element_text(size = 20),
              legend.title=element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black"),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              legend.position = c(0.87, 0.7),
              legend.background = element_blank(),
              legend.spacing.y = unit(0, "mm"),
              legend.box.background = element_rect(colour = "black"))

# plot_count_data(fill = "speed")
# plot_encounter_data(fill = "speed")
# plot_staytime_data(fill = "speed")
# plot_TTE_data(fill = "speed")

# # Plot ABM simulations
# plot_ABM()
# plot_space_use()

# # Save Results
# write.csv(D.all, paste0("Sim_results/Sim_", sim_vars$sim_names[sim_num], ".csv"))


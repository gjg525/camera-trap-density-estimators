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

file_path <- "C:/Users/guen.grosklos/Google Drive/Missoula_postdoc/Images/"
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
#                     Prop_speeds = list())
D.all <- data.frame(Model = NA,
                    Covariate = NA,
                    Est = NA,
                    SD = NA,
                    Prop_speeds = NA
)
num_runs <- 1

# Define number of clumps
num.clumps <- 100

# Define clump sizes for every clump
clump_sizes <- rep(1,num.clumps)
nind <- sum(clump_sizes)

# Landscape parms
q <- 30^2             # Number grid cells
bounds <- c(0, q^0.5) # Sampling area boundaries
t.steps <- 1000        # Number of time steps
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
ncam <- q
cam_length <- dx*.3 # length of all viewshed sides
cam.A <- cam_length ^ 2 / 2 # Assumes equilateral triangle viewsheds

################################################################################
# Define movement speeds for each cell across landscape
################################################################################

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
  
  ################################
  # Collect data
  ################################
  '%notin%' <- Negate('%in%')
  source("./R/Collect_data.R")
  
}

####################################
# Plot the stuff
####################################
plot_count_data(fill = "speed")
ggsave(
  "Count_dat_histogram.png",
  plot = last_plot(),
  path = file_path,
  scale = 1,
  width = 6,
  height = 3,
  # units = c("in", "cm", "mm", "px"),
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

plot_encounter_data(fill = "speed")
ggsave(
  "Encounter_dat_histogram.png",
  plot = last_plot(),
  path = file_path,
  scale = 1,
  width = 6,
  height = 3,
  # units = c("in", "cm", "mm", "px"),
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

plot_staytime_data(fill = "speed")
ggsave(
  "Stay_time_dat_histogram.png",
  plot = last_plot(),
  path = file_path,
  scale = 1,
  width = 6,
  height = 3,
  # units = c("in", "cm", "mm", "px"),
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

plot_TTE_data(fill = "speed")
ggsave(
  "TTE_dat_histogram.png",
  plot = last_plot(),
  path = file_path,
  scale = 1,
  width = 6,
  height = 3,
  # units = c("in", "cm", "mm", "px"),
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)


# # Plot ABM simulations
# plot_ABM()
# plot_space_use()
# plot_ABM_stay_proportions()
# plot_model_proportions()

# # Save Results
# write.csv(D.all, paste0("Sim_results/Sim_", sim_vars$sim_names[sim_num], ".csv"))
# saveRDS(D.all, paste0("Sim_results/Sim_", sim_vars$sim_names[sim_num], "_slower_speeds.rds"))


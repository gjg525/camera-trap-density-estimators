# Workflow for agent-based simulations

library(truncnorm)
library(Rfast)
library(raster)
library(tidyverse)
library(secr)
library(AHMbook)
library(R2jags)

source("./R/utils.R")
source("./R/plot_funs.R")
source("./R/ABM_sim.R")
source("./R/MLE_functions.R")
source("./R/MCMC_functions.R")

################################################################################
# Initializations
################################################################################
D.all <- c()

# Define number of clumps
num.clumps <- 100
# Define clump sizes for every clump
clump_sizes <- rep(1,num.clumps)
nind <- sum(clump_sizes)

# Landscape parms
q <- 100^2         # Number grid cells
bounds <- c(0, q^0.5) # Sampling area boundaries
# bounds <- c(0, 1) # Sampling area boundaries
t.steps <- 500        # Number of time steps
dt <- 1               # Time step size
# Grid cell lengths
dx <- (bounds[2]-bounds[1])/q^0.5
dy <- (bounds[2]-bounds[1])/q^0.5

# Random walk parms
# 0 is uncorrelated random walk, inf is ideal gas model (5 is good correlation)
default_kappa <- 0
clump.rad <- dx/2 # Tightness of clumping behavior

# Number cameras
ncam <-  200
# cam.A <- ((bounds[2]-bounds[1])/q^0.5)^2 # For now, camera area is same as cell area

# MCMC parms
n.iter <- 10000
burn.in <- 5000

################################################################################
# Define movement speeds for each cell across landscape
################################################################################
# # Create custom landscape
source("./R/Create_landscape.R")

# Calculate total area using available space
tot.A <- (bounds[2]-bounds[1])^2 * sum(!is.na(lscape_speeds$Road))/q

# Random camera sampling
cam.samps <- sample(which(!is.na(lscape_speeds$Road)), ncam, replace = F)
# Biased camera sampling
# cam.samps <- c(sample(which(lscape_speeds$Road == "On Trail"), ncam*0.7, replace = F),
#                sample(which(lscape_speeds$Road == "Off Trail"), ncam*0.3, replace = F)
#                )

# Create triangle camera viewsheds
cam_length <- dx/5
tri_cam_samps <- data.frame(cam_ID = 1:ncam,
                            lscape_index = cam.samps,
                            x = ((cam.samps-1) %% q^0.5 + 1)*dx - dx,
                            y = ceiling(cam.samps/q^0.5)*dx - dx) |>
  group_by(cam_ID) |>
  summarise(lscape_index = lscape_index,
            x = runif(1, x + cam_length*1.2, x + dx - cam_length*1.2) + c(0, -0.5*cam_length, 0.5*cam_length),
            y = runif(1, y + cam_length*1.2, y + dx - cam_length*1.2) + c(0, cam_length, cam_length),
            vertex = c(1, 2, 3),
            cam_area = calc_tri_area(x, y),
            .groups = 'drop'
  )

# Triangle viewshed cameras
dt.triangle <- data.frame(group = rep(1:ncam, each = 3),
                          polygon.x = tri_cam_samps$x,
                          polygon.y = tri_cam_samps$y)

# For now, cam area is the same for all cameras
cam.A <- cam_length^2 / 2

# Run agent-based model
animalxy.all <- ABM_sim(bounds,
                        t.steps,
                        speeds = matrix(lscape_speeds$Value, q^0.5, q^0.5),
                        direction = matrix(lscape_speeds$Direction, q^0.5, q^0.5),
                        kappa = matrix(lscape_speeds$Kappa, q^0.5, q^0.5),
                        road = matrix(lscape_speeds$Road, q^0.5, q^0.5),
                        clump_sizes,
                        clump.rad)
animalxy.all$lscape_type <- lscape_speeds$Speed[animalxy.all$XY_inds]
animalxy.all$road <- lscape_speeds$Road[animalxy.all$XY_inds]

# Plot ABM simulations
# plot_ABM()
# plot_space_use()

################################
# Collect data
################################
# TTE and staying time censors
JJ <- 20  # Occasion length
t.censor <- JJ
num.occ <- t.steps/JJ

'%notin%' <- Negate('%in%')
source("./R/Collect_data.R")

# # Plot data histograms
# plot_count_data(fill = "road")
# plot_encounter_data(fill = "road")
# plot_staytime_data(fill = "road")
# plot_TTE_data(fill = "road")
# plot_space_data(fill = "road")

# Check if any data points were collected
if(max(count_data$count) > 0) {

################################
# MLE methods
################################
################################
# REST no covariates
################################
# # print("Fit REST with MLE")
# REST.start <- c(log(mean(encounter_data$encounter)),log(mean(stay_time_data$t_stay,na.rm=T)))
# opt.REST <- optim(REST.start,
#                   REST.fn,
#                   num.encounters.dat = encounter_data$encounter,
#                   t.staying.dat = stay_time_data$t_stay,
#                   censor = t.censor,
#                   t.steps = t.steps,
#                   cam.A = cam.A,
#                   control = list(fnscale = -1, maxit = 5000),
#                   hessian = T)
# REST.u <- exp(opt.REST$par[1])
# REST.Ts <- exp(opt.REST$par[2])
# D.REST.MLE <- tot.A*REST.u

################################
# TTE no covariates
################################
# # print("Fit TTE with MLE")
# TTE.start <- c(log(mean(TTE_data$TTE, na.rm = T)),
#                log(mean(stay_time_data$t_stay,na.rm=T)))
# opt.TTE <- optim(TTE.start,
#                  TTE.fn,
#                  TTE.dat = TTE_data$TTE,
#                  t.staying.dat = stay_time_data$t_stay,
#                  censor = t.censor,
#                  cam.A = cam.A,
#                  control = list(fnscale = -1, maxit = 5000),
#                  hessian = T)
#
# TTE.u <- exp(opt.TTE$par[1])
# TTE.Ts <- exp(opt.TTE$par[2])
# D.TTE.MLE <- TTE.u*tot.A

################################
# MCT no covariates
################################
# # print("Fit Mean Count with MLE")
# MCT.start <- log(mean(count_data$count))
# opt.MCT <- optim(MCT.start,
#                  MCT.fn,
#                  cam.counts = count_data$count,
#                  control = list(fnscale = -1, maxit = 5000),
#                  hessian = T,
#                  method = "Brent", lower = -10, upper = 10)
# MCT.d <- exp(opt.MCT$par)
# D.MCT.MLE <- tot.A*MCT.d/(cam.A*t.steps)


####################################
# STE no covariates
####################################
# # print("Fit STE with MLE")
# STE.start <- log(1/mean(STE_data$STE, na.rm = T))
# opt.STE <- optim(STE.start,
#                  STE.fn,
#                  STE.dat = STE_data$STE,
#                  censor = ncam*cam.A,
#                  control = list(fnscale = -1, maxit = 5000),
#                  hessian = T,
#                  method = "Brent",
#                  lower = -10,
#                  upper = 10)
#
# u.STE <- exp(opt.STE$par)
# D.STE.MLE <- u.STE*tot.A
# # D.STE.MLE <- u.STE*tot.A*mean_clump


################################
# MCMC methods
################################
########################################
## EEDE w/ spatial covariates
########################################
# # print("Fit EEDE with MCMC")
# ptm <- proc.time()
# # unpack tidyr if extract has no applicable method
# # .rs.unloadPackage("tidyr")
# chain.EEDE.cov <-fit.model.mcmc.EEDE.cov(
#   n.iter = n.iter,
#   gamma.start = gamma.EEDE.start,
#   kappa.start = kappa.EEDE.start,
#   gamma.prior.var = gamma.EEDE.prior.var,
#   kappa.prior.var = kappa.EEDE.prior.var,
#   gamma.tune = gamma.EEDE.tune,
#   kappa.tune = kappa.EEDE.tune,
#   cam.counts = cam.counts.sum,
#   t.staying.dat = t.staying.dat,
#   spatial.covariates = spatial.covariates,
#   covariates.index = covariates.index,
#   cam.A = cam.A,
#   censor = t.censor)
# proc.time() - ptm
#
# ## Posterior summaries
# pop.ind.EEDE.cov <- which(names(chain.EEDE.cov) == "u")
# MCMC.parms.EEDE.cov <- as.mcmc(do.call(cbind, chain.EEDE.cov[-pop.ind.EEDE.cov])[-c(1:burn.in), ])
# summary(MCMC.parms.EEDE.cov)
#
# if(any(colMeans(chain.EEDE.cov$accept[burn.in:n.iter,])<0.3) || any(colMeans(chain.EEDE.cov$accept[burn.in:n.iter,])>0.6)){
#   warning(('EEDE accept rate OOB'))
# }
#
# # # plot(chain.EEDE.cov$tot.u[burn.in:n.iter])
# D.EEDE <- exp(mean(chain.EEDE.cov$gamma[burn.in:n.iter]))/t.steps
# phi.EEDE <- exp(Z%*%colMeans(chain.EEDE.cov$kappa[burn.in:n.iter,]))
# u.EEDE <- D.EEDE*phi.EEDE/sum(phi.EEDE)
#
# D.EEDE.cov <- mean(chain.EEDE.cov$tot.u[burn.in:n.iter])
# SD.EEDE.cov <- sd(chain.EEDE.cov$tot.u[burn.in:n.iter])

###################################
# REST no covariates
###################################
# print("Fit REST with MCMC, no covariates")
ptm <- proc.time()
# unpack tidyr if extract has no applicable method
# .rs.unloadPackage("tidyr")
chain.REST <-fit.model.mcmc.REST(
  n.iter = n.iter,
  gamma.start = log(mean(encounter_data$encounter)),
  kappa.start = log(mean(stay_time_data$t_stay,na.rm=T)),
  gamma.prior.var = 10^6,
  kappa.prior.var = 10^6,
  gamma.tune = -1,
  kappa.tune = -1,
  num.encounters.dat = encounter_data$encounter,
  t.staying.dat = stay_time_data$t_stay,
  t.steps = t.steps,
  cam.A = cam.A,
  censor = t.censor)
proc.time() - ptm

## Posterior summaries
MCMC.parms.REST <- as.mcmc(do.call(cbind, chain.REST)[-c(1:burn.in), ])
summary(MCMC.parms.REST)

# plot(chain.REST$tot.u[burn.in:n.iter])
D.REST.MCMC <- mean(chain.REST$tot.u[burn.in:n.iter])
SD.REST.MCMC <- sd(chain.REST$tot.u[burn.in:n.iter])

if(any(colMeans(chain.REST$accept[burn.in:n.iter,])<0.3) || any(colMeans(chain.REST$accept[burn.in:n.iter,])>0.6)){
  warning(('REST accept rate OOB'))
}


###################################
# REST w/ spatial covariates
###################################
# # print("Fit REST with MCMC")
# ptm <- proc.time()
# # unpack tidyr if extract has no applicable method
# # .rs.unloadPackage("tidyr")
# chain.REST.cov <-fit.model.mcmc.REST.cov(
#   n.iter = n.iter,
#   gamma.start = gamma.REST.start,
#   kappa.start = kappa.REST.start,
#   gamma.prior.var = gamma.REST.prior.var,
#   kappa.prior.var = kappa.REST.prior.var,
#   gamma.tune = gamma.REST.tune,
#   kappa.tune = kappa.REST.tune,
#   num.encounters.dat = num.encounters.dat,
#   t.staying.dat = t.staying.dat,
#   spatial.covariates = spatial.covariates,
#   covariates.index = covariates.index,
#   t.steps = t.steps,
#   cam.A = cam.A,
#   censor = t.censor)
# proc.time() - ptm
#
# ## Posterior summaries
# pop.ind.REST.cov <- which(names(chain.REST.cov) == "u")
# MCMC.parms.REST.cov <- as.mcmc(do.call(cbind, chain.REST.cov[-pop.ind.REST.cov])[-c(1:burn.in), ])
# summary(MCMC.parms.REST.cov)
#
# # plot(chain.REST.cov$tot.u[burn.in:n.iter])
# theta.REST <- exp(Z%*%colMeans(chain.REST.cov$gamma[burn.in:n.iter,]))
# psi.REST <- exp(Z%*%colMeans(chain.REST.cov$kappa[burn.in:n.iter,]))
# u.REST <- theta.REST
#
# D.REST.MCMC.cov <- mean(chain.REST.cov$tot.u[burn.in:n.iter])
# SD.REST.MCMC.cov <- sd(chain.REST.cov$tot.u[burn.in:n.iter])
#
# if(any(colMeans(chain.REST.cov$accept[burn.in:n.iter,])<0.3) || any(colMeans(chain.REST.cov$accept[burn.in:n.iter,])>0.6)){
#   warning(('REST accept rate OOB'))
# }

########################################
## TTE no covariates
########################################
# print("Fit TTE with MCMC, no covariates")
ptm <- proc.time()
# unpack tidyr if extract has no applicable method
# .rs.unloadPackage("tidyr")
chain.TTE <-fit.model.mcmc.TTE(
  n.iter = n.iter,
  gamma.start = log(mean(TTE_data$TTE, na.rm = T)),
  kappa.start = log(mean(stay_time_data$t_stay,na.rm=T)),
  gamma.prior.var = 10^6,
  kappa.prior.var = 10^6,
  gamma.tune = -1,
  kappa.tune = -1,
  TTE.dat = TTE_data$TTE,
  t.staying.dat = stay_time_data$t_stay,
  censor = t.censor,
  JJ = JJ,
  cam.A = cam.A)
proc.time() - ptm

## Posterior summaries
MCMC.parms.TTE <- as.mcmc(do.call(cbind, chain.TTE)[-c(1:burn.in), ])
summary(MCMC.parms.TTE)

if(any(colMeans(chain.TTE$accept[burn.in:n.iter,])<0.3) || any(colMeans(chain.TTE$accept[burn.in:n.iter,])>0.6)){
  warning(('TTE accept rate OOB'))
}

# plot(chain.TTE$tot.u[burn.in:n.iter])
D.TTE.MCMC <- mean(chain.TTE$tot.u[burn.in:n.iter])
SD.TTE.MCMC <- sd(chain.TTE$tot.u[burn.in:n.iter])


########################################
## TTE w/ spatial covariates
########################################
# # print("Fit TTE with MCMC")
# ptm <- proc.time()
# # unpack tidyr if extract has no applicable method
# # .rs.unloadPackage("tidyr")
# chain.TTE.cov <-fit.model.mcmc.TTE.cov(
#   n.iter = n.iter,
#   gamma.start = gamma.TTE.start,
#   kappa.start = kappa.TTE.start,
#   gamma.prior.var = gamma.TTE.prior.var,
#   kappa.prior.var = kappa.TTE.prior.var,
#   gamma.tune = gamma.TTE.tune,
#   kappa.tune = kappa.TTE.tune,
#   TTE.dat = TTE.dat,
#   t.staying.dat = t.staying.dat,
#   censor = t.censor,
#   spatial.covariates = spatial.covariates,
#   covariates.index = covariates.index,
#   JJ = JJ,
#   cam.A = cam.A)
# proc.time() - ptm
#
# ## Posterior summaries
# pop.ind.TTE.cov <- which(names(chain.TTE.cov) == "u")
# MCMC.parms.TTE.cov <- as.mcmc(do.call(cbind, chain.TTE.cov[-pop.ind.TTE.cov])[-c(1:burn.in), ])
# summary(MCMC.parms.TTE.cov)
#
# if(any(colMeans(chain.TTE.cov$accept[burn.in:n.iter,])<0.3) || any(colMeans(chain.TTE.cov$accept[burn.in:n.iter,])>0.6)){
#   warning(('TTE accept rate OOB'))
# }
#
# # plot(chain.TTE.cov$tot.u[burn.in:n.iter])
# u.TTE <- exp(Z%*%colMeans(chain.TTE.cov$gamma[burn.in:n.iter,]))
# psi.TTE <- exp(Z%*%colMeans(chain.TTE.cov$kappa[burn.in:n.iter,]))
#
# D.TTE.MCMC.cov <- mean(chain.TTE.cov$tot.u[burn.in:n.iter])
# SD.TTE.MCMC.cov <- sd(chain.TTE.cov$tot.u[burn.in:n.iter])


########################################
## MCT no covariates
########################################
# print("Fit Mean Count model with MCMC, no covariates")
ptm <- proc.time()
# unpack tidyr if extract has no applicable method
# .rs.unloadPackage("tidyr")
chain.MCT <- fit.model.mcmc.MCT(
  n.iter = n.iter,
  gamma.start = mean(count_data$count),
  gamma.prior.var = 10^6,
  gamma.tune = -1,
  cam.counts = count_data$count)
proc.time() - ptm

## Posterior summaries
MCMC.parms.MCT <- as.mcmc(do.call(cbind, chain.MCT)[-c(1:burn.in), ])
summary(MCMC.parms.MCT)

if(mean(chain.MCT$accept[burn.in:n.iter,])<0.3 || mean(chain.MCT$accept[burn.in:n.iter,])>0.6){
  warning(('Mean Count accept rate OOB'))
}

# plot(chain.MCT$tot.u)
D.MCT.MCMC <- mean(chain.MCT$tot.u[burn.in:n.iter])
SD.MCT.MCMC <- sd(chain.MCT$tot.u[burn.in:n.iter])


########################################
## MCT w/ spatial covariates
########################################
# # print("Fit Mean Count model with MCMC")
# ptm <- proc.time()
# # unpack tidyr if extract has no applicable method
# # .rs.unloadPackage("tidyr")
# chain.MCT.cov <-fit.model.mcmc.MCT.cov(
#   n.iter = n.iter,
#   gamma.start = gamma.MCT.start,
#   gamma.prior.var = gamma.MCT.prior.var,
#   gamma.tune = gamma.MCT.tune,
#   cam.counts = cam.counts.sum,
#   spatial.covariates = spatial.covariates,
#   covariates.index = covariates.index)
# proc.time() - ptm
#
# ## Posterior summaries
# pop.ind.MCT.cov <- which(names(chain.MCT.cov) == "u")
# MCMC.parms.MCT.cov <- as.mcmc(do.call(cbind, chain.MCT.cov[-pop.ind.MCT.cov])[-c(1:burn.in), ])
# summary(MCMC.parms.MCT.cov)
#
# if(any(colMeans(chain.MCT.cov$accept[burn.in:n.iter,])<0.3) || any(colMeans(chain.MCT.cov$accept[burn.in:n.iter,])>0.6)){
#   warning(('Mean Count accept rate OOB'))
# }
#
# # plot(chain.MCT.cov$tot.u[(n.iter-1000):n.iter])
# u.MCT <- exp(Z%*%colMeans(chain.MCT.cov$gamma[burn.in:n.iter,]))
#
# D.MCT.MCMC.cov <- mean(chain.MCT.cov$tot.u[burn.in:n.iter])
# SD.MCT.MCMC.cov <- sd(chain.MCT.cov$tot.u[burn.in:n.iter])

########################################
## STE no covariates
########################################
# print("Fit STE with MCMC, no covariates")
ptm <- proc.time()
# unpack tidyr if extract has no applicable method
# .rs.unloadPackage("tidyr")
chain.STE <-fit.model.mcmc.STE(
  n.iter = n.iter,
  gamma.start = mean(STE_data$STE, na.rm = T),
  gamma.prior.var = 10^6,
  gamma.tune = -1,
  STE.dat = STE_data$STE,
  censor = ncam*cam.A,
  cam.A = cam.A)
proc.time() - ptm

## Posterior summaries
# plot(chain.STE$tot.u[burn.in:n.iter])
MCMC.parms.STE <- as.mcmc(do.call(cbind, chain.STE)[-c(1:burn.in), ])
summary(MCMC.parms.STE)

if(mean(chain.STE$accept[burn.in:n.iter])<0.3 || mean(chain.STE$accept[burn.in:n.iter,])>0.6){
  warning(('STE accept rate OOB'))
}

D.STE.MCMC <- mean(chain.STE$tot.u[burn.in:n.iter])
SD.STE.MCMC <- sd(chain.STE$tot.u[burn.in:n.iter])

D.STE.Adjust <- D.STE.MCMC*mean_clump

####################################
# SECR no covariates
####################################
# chain.SECR <-fit.model.mcmc.SECR(as.matrix(SECR_data$dat),
#                                  psi.start = 0.5,
#                                  sigma.start = 1.5,
#                                  lambda_0.start = 0.5,
#                                  S.tune = .01,
#                                  sigma.tune = -1,
#                                  lambda_0.tune = -1.5,
#                                  X = as.matrix(SECR_data$X),
#                                  M = M,
#                                  x = bounds,
#                                  y = bounds,
#                                  n.iter = n.iter)
#
# # # Royle original mcmc sampler
# # SCR_OG <- SCR0pois(y = as.matrix(SECR_data$dat),
# #            X = as.matrix(SECR_data$X),
# #            M = M,
# #            xl = bounds[1],
# #            xu = bounds[2],
# #            yl = bounds[1],
# #            yu = bounds[2],
# #            delta = c(0.1, 0.1, .01),
# #            niter = n.iter)
#
# ## Posterior summaries
# # plot(chain.SECR$tot.u[burn.in:n.iter])
# # plot(chain.SECR$lambda_0[burn.in:n.iter])
# # plot(chain.SECR$sigma[burn.in:n.iter])
# # plot(chain.SECR$psi[burn.in:n.iter])
# # plot(chain.SECR$lambda_0[burn.in:n.iter], chain.SECR$sigma[burn.in:n.iter])
# MCMC.parms.SECR <- as.mcmc(do.call(cbind, chain.SECR)[-c(1:burn.in), ])
# summary(MCMC.parms.SECR)
#
# if(mean(chain.SECR$accept[burn.in:n.iter])<0.3 || mean(chain.SECR$accept[burn.in:n.iter,])>0.6){
#   warning(('SECR accept rate OOB'))
# }
#
# D.SECR.MCMC <- mean(chain.SECR$tot.u[burn.in:n.iter])
# SD.SECR.MCMC <- sd(chain.SECR$tot.u[burn.in:n.iter])
D.SECR.MCMC <- NA
SD.SECR.MCMC <- NA


####################################
# Plot the stuff
####################################
D.all <- data.frame(Model = c("REST", "TTE", "MCT", "STE", "SECR"),
                    Est = c(D.REST.MCMC, D.TTE.MCMC, D.MCT.MCMC, D.STE.MCMC, D.SECR.MCMC),
                    SD = c(SD.REST.MCMC, SD.TTE.MCMC, SD.MCT.MCMC, SD.STE.MCMC, SD.SECR.MCMC))
D.all$Model <- factor(D.all$Model, levels = c("REST", "TTE", "MCT", "STE", "SECR"))

} else{
  warning("No data points")
  D.all <- data.frame(Model = c("REST", "TTE", "MCT", "STE", "SECR"),
                      Est = rep(NA, 5),
                      SD = rep(NA, 5))
  D.all$Model <- factor(D.all$Model, levels = c("REST", "TTE", "MCT", "STE", "SECR"))
}

plot_onerun_results()

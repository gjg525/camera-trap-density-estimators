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
  sim_names = c("Original", "Slow_landscape", "Medium_landscape", "Fast_landscape", "Slow_cams", "Medium_cams", "Fast_cams"),
  lscape_tag = c("Random", rep("Homogeneous", 3), rep("Random", 3)),
  all_speed = c(1, 1, 2, 3, rep(1, 3)),
  cam.dist.set = c(rep(1, 4), 2, 3, 4)
)

load(file = paste0("Sim_results/Sim_", sim_vars$sim_names[sim_num], "_all_vars.RData"))

lscape_tag <- sim_vars$lscape_tag[sim_num]
all_speed <- sim_vars$all_speed[sim_num] # 1: Slow, 2: Medium, 3: Fast

# Cam sample designs (1: random, 2-4: biased slow, medium, fast)
cam.dist.set <- sim_vars$cam.dist.set[sim_num]

cam.samps <- unique(tri_cam_samps$lscape_index)

num_runs <- 1000

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
n.iter <- 20000
burn.in <- 10000

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


################################
# Collect data
################################
'%notin%' <- Negate('%in%')
source("./R/Collect_data.R")

# Quick calculation of stay time priors
stay_time_summary <- stay_time_raw %>% 
  dplyr::group_by(speed) %>% 
  dplyr::summarise(mu = mean(log(t_stay / cam.A * dx * dy)),
                   sd = sd(log(t_stay / cam.A * dx * dy))) %>% 
  dplyr::mutate(speed = factor(speed, levels = c("Slow", "Medium", "Fast"))) %>% 
  dplyr::arrange(speed)

kappa.prior.mu.tdst <- stay_time_summary %>% 
  dplyr::pull(mu)
kappa.prior.var.tdst <- stay_time_summary %>% 
  dplyr::pull(sd)

# Run Original models
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

plot(chain.TDST$tot.u[burn.in:n.iter])
D.TDST.MCMC <- mean(chain.TDST$tot.u[burn.in:n.iter])
SD.TDST.MCMC <- sd(chain.TDST$tot.u[burn.in:n.iter])

D.chain <- tibble::tibble(
  Model = "TDST no priors",
  Est = D.TDST.MCMC,
  SD = SD.TDST.MCMC,
  Kappa_mean = list(colMeans(chain.TDST$kappa[burn.in:n.iter,])),
  Kappa_sd = list(apply(chain.TDST$kappa[burn.in:n.iter,], 2, sd))
)


chain.TDST.prior <- fit.model.mcmc.TDST.cov(
  n.iter = n.iter,
  gamma.start = log(mean(count_data$count)),
  kappa.start = rep(log(mean(stay_time_data,na.rm=T)), 3),
  gamma.prior.var = 10^6,
  kappa.prior.mu = kappa.prior.mu.tdst,
  kappa.prior.var = kappa.prior.var.tdst,
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

plot(chain.TDST.prior$tot.u[burn.in:n.iter])
D.TDST.prior.MCMC <- mean(chain.TDST.prior$tot.u[burn.in:n.iter])
SD.TDST.prior.MCMC <- sd(chain.TDST.prior$tot.u[burn.in:n.iter])

D.chain <- dplyr::bind_rows(D.chain,
                            tibble::tibble(
                              Model = "TDST priors",
                              Est = D.TDST.prior.MCMC,
                              SD = SD.TDST.prior.MCMC,
                              Kappa_mean = list(colMeans(chain.TDST.prior$kappa[burn.in:n.iter,])),
                              Kappa_sd = list(apply(chain.TDST.prior$kappa[burn.in:n.iter,], 2, sd))
                            )
)


# REMOVE STAY TIME DATA INFORMATION
fit.model.mcmc.TDST.cov <- function(n.iter,
                                    gamma.start,
                                    kappa.start,
                                    gamma.prior.var,
                                    kappa.prior.mu = 0,
                                    kappa.prior.var,
                                    gamma.tune,
                                    kappa.tune,
                                    cam.counts,
                                    t.staying.dat,
                                    covariate_labels,
                                    covariates.index,
                                    cam.A,
                                    cell.A,
                                    censor,
                                    t.steps = t.steps){
  
  # Variables that will be saved
  gamma <- matrix(,n.iter+1,1)
  kappa <- matrix(,n.iter+1,sum(covariates.index==1))
  tot.u <- matrix(,n.iter+1,1)
  accept <- matrix(, n.iter+1, 1 + sum(covariates.index==1))
  gamma[1] <- gamma.start
  colnames(gamma) <- "gamma"
  kappa[1,] <- kappa.start
  colnames(kappa) <- paste0("kappa.", covariate_labels)
  colnames(tot.u) <- "Total estimate"
  colnames(accept) <- c("accept.rate.gamma",
                        paste0("accept.rate.kappa.", covariate_labels))
  
  # Account for censored times
  t.staying.dat.all <- t.staying.dat
  t.staying.dat.all[t.staying.dat.all>=censor] <- NA
  t.staying.dat.censor <- t.staying.dat
  t.staying.dat.censor[t.staying.dat.censor<censor] <- NA
  
  # Initialize with landscape-scale covariates
  beta <- exp(gamma[1])
  phi <- exp(Z%*%kappa[1,])
  u <- beta*phi/sum(phi)
  
  tune.check <- 100
  batch_n <- 0
  
  # Begin MCMC loop
  # prog_bar <- txtProgressBar(min = 0,max = n.iter,style = 3,width = 50,char = "=")
  for(i in 1:n.iter){
    # setTxtProgressBar(prog_bar, i)
    #Sample gamma
    gamma.star <- rnorm(1,gamma[i],exp(2*gamma.tune))
    beta.star <- exp(gamma.star)
    u.star <- beta.star*phi/sum(phi)
    
    # repeat estimated parms for fitting
    u.star.cams <- u.star[cam.samps] * cam.A / cell.A
    u.cams <- u[cam.samps] * cam.A / cell.A
    
    if(all(u.star.cams>0)& all(!is.infinite(u.star.cams))){
      mh1 <- sum(dpois(cam.counts,u.star.cams,log=TRUE),na.rm=TRUE) +
        sum(dnorm(gamma.star,0,gamma.prior.var^0.5,log=TRUE))
      mh2 <- sum(dpois(cam.counts,u.cams,log=TRUE),na.rm=TRUE) +
        sum(dnorm(gamma[i,],0,gamma.prior.var^0.5,log=TRUE))
      mh <- exp(mh1-mh2)
      
      if(mh>runif(1)){gamma[i+1] <- gamma.star;
      accept[i+1,1] <- 1;
      u <- u.star;
      beta <- beta.star}
      else{
        gamma[i+1] <- gamma[i];
        accept[i+1,1] <- 0}
    } else{
      gamma[i+1] <- gamma[i];
      accept[i+1,1] <- 0}
    
    beta <- exp(gamma[i+1])
    
    #Sample kappa
    kappa.temp <- kappa[i,]
    for(kk in 1:sum(covariates.index==1)){
      kappa.star <- kappa[i,]
      kappa.star[kk] <- rnorm(1,kappa[i,kk],exp(2*kappa.tune[kk]))
      phi.star <- exp(Z%*%kappa.star)
      u.star <- beta*phi.star/sum(phi.star)
      
      # # repeat estimated parms for fitting
      phi.star.cams <- phi.star[cam.samps] * cam.A / cell.A
      phi.star.cam.rep <- matrix(rep(phi.star.cams, dim(t.staying.dat.all)[2]),
                                 nrow = ncam,ncol = dim(t.staying.dat.all)[2])
      phi.all.cams <- phi[cam.samps] * cam.A / cell.A
      phi.all.cam.rep <- matrix(rep(phi.all.cams,dim(t.staying.dat.all)[2]),
                                nrow = ncam,ncol = dim(t.staying.dat.all)[2])
      
      
      if(all(1/phi.star.cams>0) & all(!is.infinite(1/phi.star.cams))){
        mh1 <-sum(dnorm(kappa.star,kappa.prior.mu,kappa.prior.var^0.5,log=TRUE))
        mh2 <-sum(dnorm(kappa[i,],kappa.prior.mu,kappa.prior.var^0.5,log=TRUE))
        mh <- exp(mh1-mh2)
        # if(all(phi.star.cams>0)){
        #   mh1 <- sum(dgamma(t.staying.dat.all, phi.star.cam.rep,log=TRUE),na.rm=TRUE) +
        #     sum(pgamma(t.staying.dat.censor, phi.star.cam.rep, lower.tail = F, log = TRUE),na.rm=TRUE) +
        #     sum(dnorm(kappa.star,0,kappa.prior.var^0.5,log=TRUE))
        #   mh2 <- sum(dgamma(t.staying.dat.all, phi.all.cam.rep,log=TRUE),na.rm=TRUE) +
        #     sum(pgamma(t.staying.dat.censor, phi.all.cam.rep, lower.tail = F, log = TRUE),na.rm=TRUE) +
        #     sum(dnorm(kappa[i,],0,kappa.prior.var^0.5,log=TRUE))
        #   mh <- exp(mh1-mh2)
        
        if(mh>runif(1)){
          kappa[i,] <- kappa.star;
          accept[i+1,1+kk] <- 1;
          u <-u.star;
          phi <- phi.star} else{
            kappa[i,] <- kappa[i,];
            accept[i+1,1+kk] <- 0}
      } else{
        kappa[i,] <- kappa[i,];
        accept[i+1,1+kk] <- 0}
    }
    kappa[i+1,] <- kappa[i,];
    phi <- exp(Z%*%kappa[i+1,])
    
    tot.u[i+1] <- sum(u)/t.steps
    
    # Update tuning parms
    if(i%%tune.check == 0){
      batch_n <- batch_n+1
      delta_n <- batch_n^-1
      # delta_n <- min(0.01,batch_n^-1)
      accept.gamma.check <- mean(accept[(i-tune.check+1):i,1],na.rm=T)
      gamma.tune[which(accept.gamma.check>0.44)] <- gamma.tune[which(accept.gamma.check>0.44)] + delta_n
      gamma.tune[which(accept.gamma.check<=0.44)] <- gamma.tune[which(accept.gamma.check<=0.44)] - delta_n
      accept.kappa.check <- colMeans(accept[(i-tune.check+1):i,2:(1+sum(covariates.index==1))],na.rm=T)
      kappa.tune[which(accept.kappa.check>0.44)] <- kappa.tune[which(accept.kappa.check>0.44)] + delta_n
      kappa.tune[which(accept.kappa.check<=0.44)] <- kappa.tune[which(accept.kappa.check<=0.44)] - delta_n
    }
  }
  # print("MCMC complete")
  
  list(accept = accept,gamma = gamma,kappa = kappa,tot.u = tot.u, u = u)
}

chain.TDST.nodata <- fit.model.mcmc.TDST.cov(
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

plot(chain.TDST.nodata$tot.u[burn.in:n.iter])
ct <- chain.TDST.nodata$tot.u[burn.in:n.iter]
ct <- ct[!is.na(ct)]
ct <- ct[!is.infinite(ct)]
D.TDST.nodata.MCMC <- mean(ct)
SD.TDST.nodata.MCMC <- sd(ct)

D.chain <- dplyr::bind_rows(D.chain,
                            tibble::tibble(
                              Model = "TDST no priors no data",
                              Est = D.TDST.nodata.MCMC,
                              SD = SD.TDST.nodata.MCMC,
                              Kappa_mean = list(colMeans(chain.TDST.nodata$kappa[burn.in:n.iter,])),
                              Kappa_sd = list(apply(chain.TDST.nodata$kappa[burn.in:n.iter,], 2, sd))
                            )
)


chain.TDST.prior.nodata <- fit.model.mcmc.TDST.cov(
  n.iter = n.iter,
  gamma.start = log(mean(count_data$count)),
  kappa.start = rep(log(mean(stay_time_data,na.rm=T)), 3),
  gamma.prior.var = 10^6,
  kappa.prior.mu = kappa.prior.mu.tdst,
  kappa.prior.var = kappa.prior.var.tdst,
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

plot(chain.TDST.prior.nodata$tot.u[burn.in:n.iter])
D.TDST.nodata.prior.MCMC <- mean(chain.TDST.prior.nodata$tot.u[burn.in:n.iter])
SD.TDST.nodata.prior.MCMC <- sd(chain.TDST.prior.nodata$tot.u[burn.in:n.iter])

D.chain <- dplyr::bind_rows(D.chain,
                            tibble::tibble(
                              Model = "TDST priors no data",
                              Est = D.TDST.nodata.prior.MCMC,
                              SD = SD.TDST.nodata.prior.MCMC,
                              Kappa_mean = list(colMeans(chain.TDST.prior.nodata$kappa[burn.in:n.iter,])),
                              Kappa_sd = list(apply(chain.TDST.prior.nodata$kappa[burn.in:n.iter,], 2, sd))
                            )
)

D.chain %>% 
  ggplot2::ggplot(ggplot2::aes(x = Model, y = Est)) +
  ggplot2::geom_point() +
  ggplot2::geom_errorbar(ggplot2::aes(ymin=Est-SD, ymax=Est+SD), width=.2, 
                         position=position_dodge(0.05))



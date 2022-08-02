# # Main script

library(tidyverse)
# library(pbapply) 
library(msm) 
library(MASS)
library(Directional)
library(raster)
# library(fields)
library(RColorBrewer)
# library(rasterVis)
library(truncnorm)
# library(mvtnorm)
# library(gridExtra)
library(fBasics)
library(coda)
# library(Rlab)
# library(fitdistrplus)
# library(patchwork)
library(crayon)
# library(GMCM)
library(beepr)
options(error = beep)

gc() # take out the trash
set.seed(1)

# save files for plotting
fig_dir <- "C:/Users/guen.grosklos/Google Drive/Missoula_postdoc/Code/All_models/"

# Define functions for MCMC and MLE samplers
source("./Main/MCMC_functions.R")
source("./Main/MLE_functions.R")

###########################
# Initialize
###########################
# Define camera sample variations
# Cam sample design (1: random, 2: 80% slow, 3: 80% medium, 4: 80% fast)
cs.all <- c(1)
cam.dist.labels <- c("random","slow","med","fast")
cam.dist.labels.caps <- c("Random","Slow","Medium","Fast")

# Define landscape variations
# 1: all slow, 2: all medium, 3: all fast, 4: equal slow, medium, fast 5: 80% fast
lv.all <- c(3)
lv.labels <- c("_slow_lscape_all","_med_lscape_all","_fast_lscape_all","","_fast_lscape")

# Variable parms
num.clumps <- 100
clump.size <- rep(1,num.clumps)
nind<-sum(clump.size)
num.runs <- 100 # number of repeated runs

# Correlated random walk parameter
# 0 is uncorrelated random walk, inf is ideal gas model (5 is good correlation)
corr.walk.kappa <- 5

# # Number cameras
ncam <-  20

# # Legend labels
leg1<-c("EEDE", "REST", "TTE", "MCT", "STE")
leg.props<-c("EEDE (MLE)", "EEDE (MCMC)","REST (MLE)", "REST (MCMC)", 
             "TTE (MLE)", "TTE (MCMC)", "MCT (MLE)", "MCT (MCMC)")

# Agent based model parms
q <- 30^2   # number grid cells
bounds <- c(0, q^0.5)
t.steps <- 200
dt <- 1 # time step size
dx <- (bounds[2]-bounds[1])/q^0.5 # space step size in x
dy <- (bounds[2]-bounds[1])/q^0.5 # space step size in y 
clump.rad <- dx/2 # tightness of clumping behavior

tot.A <- (bounds[2]-bounds[1])^2    # total area
cam.A <- ((bounds[2]-bounds[1])/q^0.5)^2 # For now, camera area is same as cell area

# Fitting parms for MCMC
n.iter <- 5000
burn.in <- 2000

# moving speed bounds for slow, medium, fast speeds
# speed.bounds <- dx/rbind(c(5.1,4.9),c(3.1,2.9),c(1.1,0.9))
speed.bounds <- dx/rbind(c(10.1,9.9),c(2.6,2.4),c(0.6,0.4))

# indices for all covariates including intercept in spatial.covariates
# Assumes 3 covariate types (slow, medium, fast)
covariates.index <- c(0,rep(1,3))

# Data collection parms
# Occasion length for TTE
JJ <- 10  
# TTE and staying time censor
t.censor <- JJ
num.occ <- t.steps/JJ
# STE censor
STE.censor <- ncam*cam.A

################################
# Loop through camera sample designs
################################
for (cam.dist.set in cs.all){
  for (lscape_var in lv.all){
  # Labels for saving figures
  means_label <- paste(cam.dist.labels[cam.dist.set],"_cams_means",lv.labels[lscape_var], sep = "")
  props_label <- paste(cam.dist.labels[cam.dist.set],"_cams_props",lv.labels[lscape_var], sep = "")
  
  ################################
  # Initialize multi-run results
  ################################
  D.all.MLE <- matrix(NA, nrow = num.runs, ncol = 5)
  SD.all.MLE <- matrix(NA, nrow = num.runs, ncol = 5)
  D.all.MLE.cov <- matrix(NA, nrow = num.runs, ncol = 5)
  SD.all.MLE.cov <- matrix(NA, nrow = num.runs, ncol = 5)
  
  D.all.MCMC <- matrix(NA, nrow = num.runs, ncol = 5)
  SD.all.MCMC <- matrix(NA, nrow = num.runs, ncol = 5)
  D.all.MCMC.cov <- matrix(NA, nrow = num.runs, ncol = 5)
  SD.all.MCMC.cov <- matrix(NA, nrow = num.runs, ncol = 5)
  
  D.all <- matrix(NA, nrow = num.runs, ncol = 20)
  SD.all <- matrix(NA, nrow = num.runs, ncol = 20)
  
  D.IS <- c()
  
  abm.props <- matrix(NA, nrow = num.runs, ncol = 3)
  EEDE.props.MLE <- matrix(NA, nrow = num.runs, ncol = 3)
  REST.props.MLE <- matrix(NA, nrow = num.runs, ncol = 3)
  TTE.props.MLE <- matrix(NA, nrow = num.runs, ncol = 3)
  MCT.props.MLE <- matrix(NA, nrow = num.runs, ncol = 3)
  EEDE.props <- matrix(NA, nrow = num.runs, ncol = 3)
  REST.props <- matrix(NA, nrow = num.runs, ncol = 3)
  TTE.props <- matrix(NA, nrow = num.runs, ncol = 3)
  MCT.props <- matrix(NA, nrow = num.runs, ncol = 3)
  
  # Run multiple iterations of code
  for (rr in 1:num.runs) {
    cat(green("Run", rr, "of", num.runs, "\n"))
    
    #########################################
    # Initialize landscape
    #########################################
    # Define slow, medium, and fast movement speed cells
    z.slow <- matrix(0,q^0.5,q^0.5)
    z.med <- matrix(0,q^0.5,q^0.5)
    z.fast <- matrix(runif(q,speed.bounds[3,1],speed.bounds[3,2]),q^0.5,q^0.5)
    
    # Initialize landscape with equal slow, medium, and fast cells
    num.slow.inds <- round(q/3)
    num.med.inds <- round(q/3)
    
    # Initialize landscape with mostly fast cells
    if (lscape_var == 5){
    num.slow.inds <- q/10
    num.med.inds <- q/10
    }
    
    # Calculate the number of fast cells
    num.fast.inds <- q-num.slow.inds-num.med.inds
    sample.inds <- sample(1:q,num.slow.inds+num.med.inds,replace = F)
    fast.inds <- 1:q
    fast.inds<-fast.inds[-sample.inds]
    slow.inds <- sample.inds[1:num.slow.inds]
    med.inds <- sample.inds[(num.slow.inds+1):(num.slow.inds+num.med.inds)]
    z.fast[sample.inds] <- 0
    z.slow[slow.inds] <- runif(num.slow.inds,speed.bounds[1,1],speed.bounds[1,2])
    z.med[med.inds] <- runif(num.med.inds,speed.bounds[2,1],speed.bounds[2,2])
    
    # Initialize homogeneous landscapes
    if (lscape_var == 1) {
      z.med[med.inds] <- runif(num.med.inds,speed.bounds[1,1],speed.bounds[1,2])
      z.fast[fast.inds] <- runif(num.fast.inds,speed.bounds[1,1],speed.bounds[1,2])
    } 
    else if (lscape_var == 2) {
      z.slow[slow.inds] <- runif(num.slow.inds,speed.bounds[2,1],speed.bounds[2,2])
      z.fast[fast.inds] <- runif(num.fast.inds,speed.bounds[2,1],speed.bounds[2,2])
    }
    else if (lscape_var == 3) {
      z.slow[slow.inds] <- runif(num.slow.inds,speed.bounds[3,1],speed.bounds[3,2])
      z.med[med.inds] <- runif(num.med.inds,speed.bounds[3,1],speed.bounds[3,2])
    }
    
    ## Create rasterstack object for covariates
    spatial.covariates <- raster(, nrows = q^0.5, ncols = q^0.5, xmn = 0, xmx = 1,
                                 ymn = 0, ymx = 1, crs = NA)
    spatial.covariates$cell <- c(1:q)
    spatial.covariates$slow <- c(z.slow)/dx
    spatial.covariates$med <- c(z.med)/dx
    spatial.covariates$fast <- c(z.fast)/dx
    
    Z <- as.matrix(spatial.covariates[[which(covariates.index==1)]])
    
    # Define movement speeds for each cell
    v.abm <- matrix(z.slow+z.med+z.fast,q^0.5,q^0.5)
    
    # Define staying time
    t.stay.abm <- 4/(dx*v.abm)
    
    # Create cumulative distribution to distribute population according to stay time
    t.stay.scale <- t.stay.abm/sum(t.stay.abm)
    t.stay.cdf <- rep(0,q+1)
    for (ii in 1:q){
      t.stay.cdf[ii+1] <- sum(t.stay.scale[1:ii]) 
    }
    
    #########################################
    # Random walks for nind animals
    #########################################
    # Run agent-based model
    source("./Main/ABM_sims.R")
    
    if(cam.dist.set == 1){
      # # 1-cell sized cameras, randomly selected
      cam.samps <- sample(1:q,ncam, replace=FALSE)
      # indices for each covariate type
      cam.slow <- which(cam.samps %in% slow.inds)
      cam.med <- which(cam.samps %in% med.inds)
      cam.fast <- which(cam.samps %in% fast.inds)
      cam.props <- c(length(cam.slow),length(cam.med),length(cam.fast))/ncam
      cam.props.rounds <- round(cam.props*100)/100
      
      
      cam.counts <- u.abm.all[cam.samps] # num counts at each camera for each time
      # Make sure all landscape types are sampled
      sample_check <- c(sum(cam.counts[cam.slow]), 
                        sum(cam.counts[cam.med]), 
                        sum(cam.counts[cam.fast]))
      while(any(sample_check == 0)) {
        cam.samps <- sample(1:q,ncam, replace=FALSE)
        # indices for each covariate type
        cam.slow <- which(cam.samps %in% slow.inds)
        cam.med <- which(cam.samps %in% med.inds)
        cam.fast <- which(cam.samps %in% fast.inds)
        cam.props <- c(length(cam.slow),length(cam.med),length(cam.fast))/ncam
        cam.props.rounds <- round(cam.props*100)/100
        
        cam.counts <- u.abm.all[cam.samps] # num counts at each camera for each time
        sample_check <- c(sum(cam.counts[cam.slow]), 
                          sum(cam.counts[cam.med]), 
                          sum(cam.counts[cam.fast]))
      }
      }else{
      # # proportion sampled slow, med, fast
      ps <- ncam*c(0.1,0.1,0.1)
      ps[cam.dist.set-1] <- ncam*0.8 
      cam.samps <- c(sample(slow.inds,ps[1],replace=F),sample(med.inds,ps[2],replace=F),
                     sample(fast.inds,ps[3],replace=F))
      # indices for each covariate type
      cam.slow <- which(cam.samps %in% slow.inds)
      cam.med <- which(cam.samps %in% med.inds)
      cam.fast <- which(cam.samps %in% fast.inds)
      cam.props <- c(length(cam.slow),length(cam.med),length(cam.fast))/ncam
      cam.props.rounds <- round(cam.props*100)/100
      
      cam.counts <- u.abm.all[cam.samps] # num counts at each camera for each time
      # Make sure all landscape types are sampled
      sample_check <- c(sum(cam.counts[cam.slow]), 
                        sum(cam.counts[cam.med]), 
                        sum(cam.counts[cam.fast]))
      while(any(sample_check == 0)) {
        cam.samps <- c(sample(slow.inds,ps[1],replace=F),sample(med.inds,ps[2],replace=F),
                       sample(fast.inds,ps[3],replace=F))
        # indices for each covariate type
        cam.slow <- which(cam.samps %in% slow.inds)
        cam.med <- which(cam.samps %in% med.inds)
        cam.fast <- which(cam.samps %in% fast.inds)
        cam.props <- c(length(cam.slow),length(cam.med),length(cam.fast))/ncam
        cam.props.rounds <- round(cam.props*100)/100
        
        cam.counts <- u.abm.all[cam.samps] # num counts at each camera for each time
        sample_check <- c(sum(cam.counts[cam.slow]), 
                          sum(cam.counts[cam.med]), 
                          sum(cam.counts[cam.fast]))
      }
      }
    
    # Make sure cam.samps adds up correctly
    if(length(cam.samps)!=ncam){
      stop("camera samples not equal to number cameras")
    }
    
    # Camera statistics
    cam.counts.sum <- rowSums(cam.counts) # num counts at each camera 
    
    # number of counts from cameras for each covariate type
    cam.counts.slow <- cam.counts[cam.slow,]
    cam.counts.med <- cam.counts[cam.med,]
    cam.counts.fast <- cam.counts[cam.fast,]
    cam.counts.props <- c(sum(cam.counts.slow),sum(cam.counts.med),sum(cam.counts.fast))/sum(cam.counts)
    
    # number of counts across whole landscape for each covariate type
    slow.counts <- sum(u.abm.all[slow.inds])/t.steps
    med.counts <- sum(u.abm.all[med.inds])/t.steps
    fast.counts <- sum(u.abm.all[fast.inds])/t.steps
    
    # Check that abm distributions match predicted distributions based on staying time
    stay.time.distribution <- c(mean(t.stay.abm[slow.inds]),
                                mean(t.stay.abm[med.inds]),
                                mean(t.stay.abm[fast.inds]))
    stay.time.distribution.scale <-stay.time.distribution/sum(stay.time.distribution)
    abm.distribution <- c(slow.counts,med.counts,fast.counts)/c(num.slow.inds,num.med.inds,num.fast.inds)
    abm.distribution.scale <- abm.distribution/sum(abm.distribution)
    
    ################################
    # Collect data
    ################################
    source("./Main/collect_data.R")
  
    # Initialize start values for covariate fits based on data
    Zcam <- Z[cam.samps,]
    
    if (all(cam.counts[cam.slow,]==0) || all(cam.counts[cam.med,]==0) || all(cam.counts[cam.fast,]==0)) {
      num.enc.cov.start <- rep(log(mean(num.encounters.dat))/mean(Zcam),3)
      t.staying.cov.start <- rep(log(mean(t.staying.dat,na.rm=T))/mean(Zcam),3)
      TTE.dat.cov.start <- rep(log(mean(TTE.dat,na.rm=T)/t.steps)/mean(Zcam),3)
      cam.counts.cov.start <- rep(log(mean(cam.counts))/mean(Zcam),3)
    } else {
    num.enc.cov.start <- c((log(mean(num.encounters.dat[cam.slow])))/mean(Zcam[cam.slow,1]),
                           (log(mean(num.encounters.dat[cam.med])))/mean(Zcam[cam.med,2]),
                           (log(mean(num.encounters.dat[cam.fast])))/mean(Zcam[cam.fast,3]))
    t.staying.cov.start <- c((log(mean(t.staying.dat[cam.slow,],na.rm=T)))/mean(Zcam[cam.slow,1]),
                             (log(mean(t.staying.dat[cam.med,],na.rm=T)))/mean(Zcam[cam.med,2]),
                             (log(mean(t.staying.dat[cam.fast,],na.rm=T)))/mean(Zcam[cam.fast,3]))
    TTE.dat.cov.start <- c((log(mean(TTE.dat[cam.slow,],na.rm=T)/t.steps))/mean(Zcam[cam.slow,1]),
                           (log(mean(TTE.dat[cam.med,],na.rm=T)/t.steps))/mean(Zcam[cam.med,2]),
                           (log(mean(TTE.dat[cam.fast,],na.rm=T)/t.steps))/mean(Zcam[cam.fast,3]))
    cam.counts.cov.start <- c((log(mean(cam.counts[cam.slow,])))/mean(Zcam[cam.slow,1]),
                              (log(mean(cam.counts[cam.med,])))/mean(Zcam[cam.med,2]),
                              (log(mean(cam.counts[cam.fast,])))/mean(Zcam[cam.fast,3]))
    }
    
    
    gamma.EEDE.start <- log(mean(cam.counts))-mean(Zcam%*%c(t.staying.cov.start))
    kappa.EEDE.start <- c(t.staying.cov.start)
    gamma.EEDE.prior.var <- 10^6
    kappa.EEDE.prior.var <- 10^6
    gamma.EEDE.tune <- -1
    kappa.EEDE.tune <- c(-1,-1,-1)
  
    gamma.REST.start <- -c(num.enc.cov.start)
    kappa.REST.start <- c(t.staying.cov.start)
    gamma.REST.prior.var <- 10^6
    kappa.REST.prior.var <- 10^6
    gamma.REST.tune <- c(-1,-1,-1)
    kappa.REST.tune <-c(-1,-1,-1)
    
    gamma.TTE.start <- c(TTE.dat.cov.start)
    kappa.TTE.start <- c(t.staying.cov.start)
    gamma.TTE.prior.var <- 10^6
    kappa.TTE.prior.var <- 10^6
    gamma.TTE.tune <- c(-1,-1,-1)
    kappa.TTE.tune <- c(-1,-1,-1)
    
    gamma.MCT.start <-c(cam.counts.cov.start)
    gamma.MCT.prior.var <- 10^6
    gamma.MCT.tune <- c(-1,-1,-1)
    
    STE.prior.var <- 10^6
    STE.tune <- -1.5
    
  ################################
  # MLE methods
  ################################
    ################################
    # EEDE w/ spatial covariates
    ################################
    # print("Fit EEDE with MLE covariates")
    EEDE.cov.start <- c(gamma.EEDE.start,kappa.EEDE.start)
    opt.EEDE.cov <- optim(EEDE.cov.start,
                         EEDE.cov.fn,
                         cam.counts = cam.counts.sum,
                         t.staying.dat = t.staying.dat,
                         censor = t.censor,
                         Z = Z,
                         cam.samps = cam.samps,
                         t.steps = dim(cam.counts)[2],
                         control = list(fnscale = -1, maxit = 5000),
                         hessian = T)
    # Calculate densities
    beta.EEDE.cov <- exp(opt.EEDE.cov$par[1])
    phi.EEDE.cov <- exp(Z%*%opt.EEDE.cov$par[2:length(EEDE.cov.start)])
    u.EEDE.cov <- beta.EEDE.cov*phi.EEDE.cov/sum(phi.EEDE.cov)
    D.EEDE.MLE.cov <- sum(u.EEDE.cov)/t.steps
    
    # Variance with msm::deltamethod
    varB <- -ginv(opt.EEDE.cov$hessian)
    form <- sprintf("~ %f * exp(x1)", q)
    SD.EEDE.MLE.cov <- deltamethod(g = as.formula(form), mean = opt.EEDE.cov$par, cov = varB, ses = T)
    
    if(opt.EEDE.cov$convergence != 0){
      warning(('Convergence unsuccessful for EEDE MLE w/ covariates'))
    }
    
    ################################
    # REST no covariates
    ################################
    # print("Fit REST with MLE")
    REST.start <- c(log(mean(num.encounters.dat)),log(mean(t.staying.dat,na.rm=T)))
    opt.REST <- optim(REST.start,
                      REST.fn,
                      num.encounters.dat = num.encounters.dat,
                      t.staying.dat = t.staying.dat,
                      censor = t.censor,
                      t.steps = t.steps,
                      cam.A = cam.A,
                      control = list(fnscale = -1, maxit = 5000),
                      hessian = T)
    REST.u <- exp(opt.REST$par[1])
    REST.Ts <- exp(opt.REST$par[2])
    D.REST.MLE <- tot.A*REST.u
    
    # Variance with msm::deltamethod
    varB <- -ginv(opt.REST$hessian)
    form <- sprintf("~ %f * exp(x1)", tot.A)
    SD.REST.MLE <- deltamethod(g = as.formula(form), mean = opt.REST$par, cov = varB, ses = T)
    
    if(opt.REST$convergence != 0){
      warning(('Convergence unsuccessful for REST MLE no covariates'))
    }
    
    ################################
    # REST w/ spatial covariates
    ################################
    # print("Fit REST with MLE covariates")
    REST.cov.start <- c(gamma.REST.start,kappa.REST.start)
    opt.REST.cov <- optim(REST.cov.start,
                      REST.cov.fn,
                      num.encounters.dat = num.encounters.dat,
                      t.staying.dat = t.staying.dat,
                      censor = t.censor,
                      Z = Z,
                      cov.inds = covariates.index,
                      t.steps = t.steps,
                      control = list(fnscale = -1, maxit = 5000),
                      hessian = T)
    # Calculate densities
    u.REST.cov <- exp(Z%*%opt.REST.cov$par[1:sum(covariates.index)])
    phi.REST.cov <- exp(Z%*%opt.REST.cov$par[(sum(covariates.index)+1):(2*sum(covariates.index))])
    D.REST.MLE.cov <- sum(u.REST.cov)
    
    # Variance with msm::deltamethod
    varB <- -ginv(opt.REST.cov$hessian)
    form <- sprintf("~ %f * exp(x1)", q)
    SD.REST.MLE.cov <- deltamethod(g = as.formula(form), mean = opt.REST.cov$par, cov = varB, ses = T)
    
    if(opt.REST.cov$convergence != 0){
      warning(warning('Convergence unsuccessful for REST MLE w/ covariates'))
    }
    
    ################################
    # TTE no covariates
    ################################
    # print("Fit TTE with MLE")
    TTE.start <- c(log(mean(TTE.dat, na.rm = T)),
                   log(mean(t.staying.dat,na.rm=T)))
    opt.TTE <- optim(TTE.start,
                     TTE.fn,
                     TTE.dat = TTE.dat,
                     t.staying.dat = t.staying.dat,
                     censor = t.censor,
                     cam.A = cam.A,
                     control = list(fnscale = -1, maxit = 5000),
                     hessian = T)
    
    TTE.u <- exp(opt.TTE$par[1])
    TTE.Ts <- exp(opt.TTE$par[2])
    D.TTE.MLE <- TTE.u*tot.A
    
    # Variance with msm::deltamethod
    varB <- -ginv(opt.TTE$hessian)
    form <- sprintf("~ %f * exp(x1)", tot.A/cam.A)
    SD.TTE.MLE <- deltamethod(g = as.formula(form), mean = opt.TTE$par, cov = varB, ses = T)
    
    if(opt.TTE$convergence != 0){
      warning(warning('Convergence unsuccessful for TTE MLE no covariates'))
    }
    
    ################################
    # TTE w/ spatial covariates
    ################################
    # print("Fit TTE with MLE covariates")
    TTE.cov.start <- c(gamma.TTE.start,kappa.TTE.start)
    opt.TTE.cov <- optim(TTE.cov.start,
                         TTE.cov.fn,
                         TTE.dat = TTE.dat,
                         t.staying.dat = t.staying.dat,
                         censor = t.censor,
                         Z = Z,
                         cam.samps = cam.samps,
                         num.occ = num.occ,
                         cov.inds = covariates.index,
                         control = list(fnscale = -1, maxit = 5000),
                         hessian = T)
    # Calculate densities
    u.TTE.cov <- exp(Z%*%opt.TTE.cov$par[1:sum(covariates.index)])
    phi.TTE.cov <- exp(Z%*%opt.TTE.cov$par[(sum(covariates.index)+1):(2*sum(covariates.index))])
    D.TTE.MLE.cov <- sum(u.TTE.cov)
    
    # Variance with msm::deltamethod
    varB <- -ginv(opt.TTE.cov$hessian)
    form <- sprintf("~ %f * exp(x1)", q)
    SD.TTE.MLE.cov <- deltamethod(g = as.formula(form), mean = opt.TTE.cov$par, cov = varB, ses = T)
    
    if(opt.TTE.cov$convergence != 0){
      warning(('Convergence unsuccessful for TTE MLE w/ covariates'))
    }
    
    ################################
    # MCT no covariates 
    ################################
    # print("Fit Mean Count with MLE")
    MCT.start <- log(mean(cam.counts))
    opt.MCT <- optim(MCT.start,
                     MCT.fn,
                     cam.counts = cam.counts.sum,
                     control = list(fnscale = -1, maxit = 5000),
                     hessian = T, 
                     method = "Brent", lower = -10, upper = 10)
    MCT.d <- exp(opt.MCT$par)
    D.MCT.MLE <- tot.A*MCT.d/(cam.A*t.steps)
    
    # Variance with msm::deltamethod
    varB <- -ginv(opt.MCT$hessian)
    form <- sprintf("~ %f * exp(x1)", tot.A/cam.A)
    SD.MCT.MLE <- deltamethod(g = as.formula(form), mean = opt.MCT$par, cov = varB, ses = T)
    
    if(opt.MCT$convergence != 0){
      warning(('Convergence unsuccessful for MCT MLE no covariates'))
    }
    
    ################################
    # MCT w/ spatial covariates
    ################################
    # print("Fit Mean Count with MLE covariates")
    MCT.cov.start <- gamma.MCT.start
    opt.MCT.cov <- optim(MCT.cov.start,
                         MCT.cov.fn,
                         cam.counts = cam.counts.sum,
                         Z = Z,
                         cam.samps = cam.samps,
                         cov.inds = covariates.index,
                         control = list(fnscale = -1, maxit = 5000),
                         hessian = T)
    # Calculate densities
    u.MCT.cov <- exp(Z%*%opt.MCT.cov$par[1:sum(covariates.index)])
    D.MCT.MLE.cov <- sum(u.MCT.cov)/t.steps
    
    # Variance with msm::deltamethod
    varB <- -ginv(opt.MCT.cov$hessian)
    form <- sprintf("~ %f * exp(x1)", q)
    SD.MCT.MLE.cov <- deltamethod(g = as.formula(form), 
                                    mean = opt.MCT.cov$par, 
                                    cov = varB, ses = T)
    
    if(opt.MCT.cov$convergence != 0){
      warning(('Convergence unsuccessful for MCT MLE w/ covariates'))
    }
    
    ####################################
    # STE no covariates
    ####################################
    # print("Fit STE with MLE")
    STE.start <- log(1/mean(STE.dat, na.rm = T))
    opt.STE <- optim(STE.start,
                     STE.fn,
                     STE.dat = STE.dat,
                     censor = STE.censor,
                     control = list(fnscale = -1, maxit = 5000),
                     hessian = T,
                     method = "Brent",
                     lower = -q,
                     upper = q)
    
    u.STE <- exp(opt.STE$par)
    D.STE.MLE <- u.STE*tot.A
    # Variance with msm::deltamethod
    varB <- -ginv(opt.STE$hessian)
    form <- sprintf("~ %f * exp(x1)", tot.A/cam.A)
    SD.STE.MLE <- deltamethod(g = as.formula(form), mean = opt.STE$par, cov = varB, ses = T)
    
    if(opt.STE$convergence != 0){
      warning(('Convergence unsuccessful for STE MLE no covariates'))
    }
    
  ####################################
  # MCMC methods
  ####################################
    ########################################
    ## EEDE w/ spatial covariates
    ########################################
    # print("Fit EEDE with MCMC")
    ptm <- proc.time()
    # unpack tidyr if extract has no applicable method
    # .rs.unloadPackage("tidyr")
    chain.EEDE.cov <-fit.model.mcmc.EEDE.cov(
      n.iter = n.iter,
      gamma.start = gamma.EEDE.start,
      kappa.start = kappa.EEDE.start,
      gamma.prior.var = gamma.EEDE.prior.var,
      kappa.prior.var = kappa.EEDE.prior.var,
      gamma.tune = gamma.EEDE.tune,
      kappa.tune = kappa.EEDE.tune,
      cam.counts = cam.counts.sum,
      t.staying.dat = t.staying.dat,
      spatial.covariates = spatial.covariates,
      covariates.index = covariates.index,
      cam.A = cam.A,
      censor = t.censor)
    proc.time() - ptm
    
    ## Posterior summaries
    pop.ind.EEDE.cov <- which(names(chain.EEDE.cov) == "u")
    MCMC.parms.EEDE.cov <- as.mcmc(do.call(cbind, chain.EEDE.cov[-pop.ind.EEDE.cov])[-c(1:burn.in), ])
    summary(MCMC.parms.EEDE.cov)
    
    if(any(colMeans(chain.EEDE.cov$accept[burn.in:n.iter,])<0.3) || any(colMeans(chain.EEDE.cov$accept[burn.in:n.iter,])>0.6)){
      warning(('EEDE accept rate OOB'))
    }
    
    # # plot(chain.EEDE.cov$tot.u[burn.in:n.iter])
    beta.EEDE <- exp(mean(chain.EEDE.cov$gamma[burn.in:n.iter]))
    phi.EEDE <- exp(Z%*%colMeans(chain.EEDE.cov$kappa[burn.in:n.iter,]))
    u.EEDE <- beta.EEDE*phi.EEDE/sum(phi.EEDE)
    
    D.EEDE.cov <- mean(chain.EEDE.cov$tot.u[burn.in:n.iter])
    SD.EEDE.cov <- sd(chain.EEDE.cov$tot.u[burn.in:n.iter])
  
    ###################################
    # REST no covariates
    ###################################
    # print("Fit REST with MCMC, no covariates")
    ptm <- proc.time()
    # unpack tidyr if extract has no applicable method
    # .rs.unloadPackage("tidyr")
    chain.REST <-fit.model.mcmc.REST(
      n.iter = n.iter,
      gamma.start = log(mean(num.encounters.dat)),
      kappa.start = log(mean(t.staying.dat,na.rm=T)),
      gamma.prior.var = gamma.REST.prior.var,
      kappa.prior.var = kappa.REST.prior.var,
      gamma.tune = gamma.REST.tune[1],
      kappa.tune = kappa.REST.tune[1],
      num.encounters.dat = num.encounters.dat,
      t.staying.dat = t.staying.dat,
      spatial.covariates = spatial.covariates,
      covariates.index = covariates.index,
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
    # print("Fit REST with MCMC")
    ptm <- proc.time()
    # unpack tidyr if extract has no applicable method
    # .rs.unloadPackage("tidyr")
    chain.REST.cov <-fit.model.mcmc.REST.cov(
      n.iter = n.iter,
      gamma.start = gamma.REST.start,
      kappa.start = kappa.REST.start,
      gamma.prior.var = gamma.REST.prior.var,
      kappa.prior.var = kappa.REST.prior.var,
      gamma.tune = gamma.REST.tune,
      kappa.tune = kappa.REST.tune,
      num.encounters.dat = num.encounters.dat,
      t.staying.dat = t.staying.dat,
      spatial.covariates = spatial.covariates,
      covariates.index = covariates.index,
      t.steps = t.steps,
      cam.A = cam.A,
      censor = t.censor)
    proc.time() - ptm
    
    ## Posterior summaries
    pop.ind.REST.cov <- which(names(chain.REST.cov) == "u")
    MCMC.parms.REST.cov <- as.mcmc(do.call(cbind, chain.REST.cov[-pop.ind.REST.cov])[-c(1:burn.in), ])
    summary(MCMC.parms.REST.cov)
    
    # plot(chain.REST.cov$tot.u[burn.in:n.iter])
    theta.REST <- exp(Z%*%colMeans(chain.REST.cov$gamma[burn.in:n.iter,]))
    psi.REST <- exp(Z%*%colMeans(chain.REST.cov$kappa[burn.in:n.iter,]))
    u.REST <- theta.REST
    
    D.REST.MCMC.cov <- mean(chain.REST.cov$tot.u[burn.in:n.iter])
    SD.REST.MCMC.cov <- sd(chain.REST.cov$tot.u[burn.in:n.iter])
    
    if(any(colMeans(chain.REST.cov$accept[burn.in:n.iter,])<0.3) || any(colMeans(chain.REST.cov$accept[burn.in:n.iter,])>0.6)){
      warning(('REST accept rate OOB'))
    }
    
    ########################################
    ## TTE no covariates
    ########################################
    # print("Fit TTE with MCMC, no covariates")
    ptm <- proc.time()
    # unpack tidyr if extract has no applicable method
    # .rs.unloadPackage("tidyr")
    chain.TTE <-fit.model.mcmc.TTE(
      n.iter = n.iter,
      gamma.start = log(mean(TTE.dat, na.rm = T)),
      kappa.start = log(mean(t.staying.dat,na.rm=T)),
      gamma.prior.var = gamma.TTE.prior.var,
      kappa.prior.var = kappa.TTE.prior.var,
      gamma.tune = gamma.TTE.tune[1],
      kappa.tune = kappa.TTE.tune[1],
      TTE.dat = TTE.dat,
      t.staying.dat = t.staying.dat,
      censor = t.censor,
      spatial.covariates = spatial.covariates,
      covariates.index = covariates.index,
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
    # print("Fit TTE with MCMC")
    ptm <- proc.time()
    # unpack tidyr if extract has no applicable method
    # .rs.unloadPackage("tidyr")
    chain.TTE.cov <-fit.model.mcmc.TTE.cov(
      n.iter = n.iter,
      gamma.start = gamma.TTE.start,
      kappa.start = kappa.TTE.start,
      gamma.prior.var = gamma.TTE.prior.var,
      kappa.prior.var = kappa.TTE.prior.var,
      gamma.tune = gamma.TTE.tune,
      kappa.tune = kappa.TTE.tune,
      TTE.dat = TTE.dat,
      t.staying.dat = t.staying.dat,
      censor = t.censor,
      spatial.covariates = spatial.covariates,
      covariates.index = covariates.index,
      JJ = JJ,
      cam.A = cam.A)
    proc.time() - ptm
    
    ## Posterior summaries
    pop.ind.TTE.cov <- which(names(chain.TTE.cov) == "u")
    MCMC.parms.TTE.cov <- as.mcmc(do.call(cbind, chain.TTE.cov[-pop.ind.TTE.cov])[-c(1:burn.in), ])
    summary(MCMC.parms.TTE.cov)
    
    if(any(colMeans(chain.TTE.cov$accept[burn.in:n.iter,])<0.3) || any(colMeans(chain.TTE.cov$accept[burn.in:n.iter,])>0.6)){
      warning(('TTE accept rate OOB'))
    }
    
    # plot(chain.TTE.cov$tot.u[burn.in:n.iter])
    u.TTE <- exp(Z%*%colMeans(chain.TTE.cov$gamma[burn.in:n.iter,]))
    psi.TTE <- exp(Z%*%colMeans(chain.TTE.cov$kappa[burn.in:n.iter,]))

    D.TTE.MCMC.cov <- mean(chain.TTE.cov$tot.u[burn.in:n.iter])
    SD.TTE.MCMC.cov <- sd(chain.TTE.cov$tot.u[burn.in:n.iter])
    
    
    ########################################
    ## MCT no covariates
    ########################################
    # print("Fit Mean Count model with MCMC, no covariates")
    ptm <- proc.time()
    # unpack tidyr if extract has no applicable method
    # .rs.unloadPackage("tidyr")
    chain.MCT <-fit.model.mcmc.MCT(
      n.iter = n.iter,
      gamma.start = mean(cam.counts),
      gamma.prior.var = gamma.MCT.prior.var,
      gamma.tune = gamma.MCT.tune[1],
      cam.counts = cam.counts.sum,
      spatial.covariates = spatial.covariates,
      covariates.index = covariates.index)
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
    # print("Fit Mean Count model with MCMC")
    ptm <- proc.time()
    # unpack tidyr if extract has no applicable method
    # .rs.unloadPackage("tidyr")
    chain.MCT.cov <-fit.model.mcmc.MCT.cov(
      n.iter = n.iter,
      gamma.start = gamma.MCT.start,
      gamma.prior.var = gamma.MCT.prior.var,
      gamma.tune = gamma.MCT.tune,
      cam.counts = cam.counts.sum,
      spatial.covariates = spatial.covariates,
      covariates.index = covariates.index)
    proc.time() - ptm
    
    ## Posterior summaries
    pop.ind.MCT.cov <- which(names(chain.MCT.cov) == "u")
    MCMC.parms.MCT.cov <- as.mcmc(do.call(cbind, chain.MCT.cov[-pop.ind.MCT.cov])[-c(1:burn.in), ])
    summary(MCMC.parms.MCT.cov)
    
    if(any(colMeans(chain.MCT.cov$accept[burn.in:n.iter,])<0.3) || any(colMeans(chain.MCT.cov$accept[burn.in:n.iter,])>0.6)){
      warning(('Mean Count accept rate OOB'))
    }
    
    # plot(chain.MCT.cov$tot.u)
    u.MCT <- exp(Z%*%colMeans(chain.MCT.cov$gamma[burn.in:n.iter,]))

    D.MCT.MCMC.cov <- mean(chain.MCT.cov$tot.u[burn.in:n.iter])
    SD.MCT.MCMC.cov <- sd(chain.MCT.cov$tot.u[burn.in:n.iter])
    
    ########################################
    ## STE no covariates
    ########################################
    # print("Fit STE with MCMC, no covariates")
    ptm <- proc.time()
    # unpack tidyr if extract has no applicable method
    # .rs.unloadPackage("tidyr")
    chain.STE <-fit.model.mcmc.STE(
      n.iter = n.iter,
      gamma.start = 1/mean(STE.dat, na.rm = T),
      gamma.prior.var = STE.prior.var,
      gamma.tune = STE.tune,
      STE.dat = STE.dat,
      censor = STE.censor,
      cam.A = cam.A,
      spatial.covariates = spatial.covariates,
      covariates.index = covariates.index)
    proc.time() - ptm
    
    ## Posterior summaries
    # plot(chain.STE$tot.u[burn.in:n.iter])
    MCMC.parms.STE <- as.mcmc(do.call(cbind, chain.STE)[-c(1:burn.in), ])
    summary(MCMC.parms.STE)
    
    if(mean(chain.STE$accept[burn.in:n.iter])<0.3 || mean(chain.STE$accept[burn.in:n.iter,])>0.6){
      warning(('STE accept rate OOB'))
    }
    
    # plot(chain.STE$tot.u[burn.in:n.iter])
    D.STE.MCMC <- mean(chain.STE$tot.u[burn.in:n.iter])
    SD.STE.MCMC <- sd(chain.STE$tot.u[burn.in:n.iter])
    
    
    ####################################
    # IS (Note: same mean as Mean Count w/o covariates)
    ####################################
    D.IS[rr] <- mean(cam.counts,na.rm=T)*q
      
  ####################################
  # Accumulate all results
  ####################################
  D.all.MLE[rr,] <- c(NA,D.REST.MLE,D.TTE.MLE,D.MCT.MLE,D.STE.MLE)
  SD.all.MLE[rr,] <- c(NA,SD.REST.MLE,SD.TTE.MLE,SD.MCT.MLE,SD.STE.MLE)
  D.all.MLE.cov[rr,] <- c(D.EEDE.MLE.cov,D.REST.MLE.cov,D.TTE.MLE.cov,D.MCT.MLE.cov,NA)
  SD.all.MLE.cov[rr,] <- c(SD.EEDE.MLE.cov,SD.REST.MLE.cov,SD.TTE.MLE.cov,SD.MCT.MLE.cov,NA)
  
  D.all.MCMC[rr,] <- c(NA,D.REST.MCMC,D.TTE.MCMC,D.MCT.MCMC,D.STE.MCMC)
  SD.all.MCMC[rr,] <- c(NA,SD.REST.MCMC,SD.TTE.MCMC,SD.MCT.MCMC,SD.STE.MCMC)
  D.all.MCMC.cov[rr,] <- c(D.EEDE.cov,D.REST.MCMC.cov,D.TTE.MCMC.cov,D.MCT.MCMC.cov,NA)
  SD.all.MCMC.cov[rr,] <- c(SD.EEDE.cov,SD.REST.MCMC.cov,SD.TTE.MCMC.cov,SD.MCT.MCMC.cov,NA)
  
  D.all[rr,] <- c(D.all.MLE[rr,],D.all.MLE.cov[rr,],D.all.MCMC[rr,],D.all.MCMC.cov[rr,])
  SD.all[rr,] <- c(SD.all.MLE[rr,],SD.all.MLE.cov[rr,],SD.all.MCMC[rr,],SD.all.MCMC.cov[rr,])
  colnames(D.all) <- rep(c("EEDE","REST","TTE","MCT","STE"),4)
  colnames(SD.all) <- rep(c("EEDE","REST","TTE","MCT","STE"),4)
  
  # Calculate proportions where individuals are likely to be found
  abm.props[rr,] <- c(slow.counts,med.counts,fast.counts)/nind
  EEDE.props.MLE[rr,] <- c(sum(u.EEDE.cov[slow.inds]),sum(u.EEDE.cov[med.inds]),sum(u.EEDE.cov[fast.inds]))/sum(u.EEDE.cov[])
  REST.props.MLE[rr,] <- c(sum(u.REST.cov[slow.inds]),sum(u.REST.cov[med.inds]),sum(u.REST.cov[fast.inds]))/sum(u.REST.cov)
  TTE.props.MLE[rr,] <- c(sum(u.TTE.cov[slow.inds]),sum(u.TTE.cov[med.inds]),sum(u.TTE.cov[fast.inds]))/sum(u.TTE.cov)
  MCT.props.MLE[rr,] <- c(sum(u.MCT.cov[slow.inds]),sum(u.MCT.cov[med.inds]),sum(u.MCT.cov[fast.inds]))/sum(u.MCT.cov)
  EEDE.props[rr,] <- c(sum(u.EEDE[slow.inds]),sum(u.EEDE[med.inds]),sum(u.EEDE[fast.inds]))/sum(u.EEDE[])
  REST.props[rr,] <- c(sum(u.REST[slow.inds]),sum(u.REST[med.inds]),sum(u.REST[fast.inds]))/sum(u.REST[])
  TTE.props[rr,] <- c(sum(u.TTE[slow.inds]),sum(u.TTE[med.inds]),sum(u.TTE[fast.inds]))/sum(u.TTE[])
  MCT.props[rr,] <- c(sum(u.MCT[slow.inds]),sum(u.MCT[med.inds]),sum(u.MCT[fast.inds]))/sum(u.MCT[])
  
  }
  
  ####################################
  # Calculate Summaries
  ####################################
  D.all.Means <- colMeans(D.all)
  # colSds for some reason does not like NA vals
  D.all2 <- D.all
  D.all2[is.na(D.all2)] <- 0
  D.all.Sds <- colSds(D.all2)
  D.all.Sds[D.all.Sds == 0] <- NA
  
  D.all.MLE.Means <- colMeans(D.all.MLE)
  # colSds for some reason does not like NA vals
  D.all2.MLE <- D.all.MLE
  D.all2.MLE[is.na(D.all2.MLE)] <- 0
  D.all.MLE.Sds <- colSds(D.all2.MLE)
  D.all.MLE.Sds[D.all.MLE.Sds == 0] <- NA
  
  D.all.MLE.cov.Means <- colMeans(D.all.MLE.cov)
  # colSds for some reason does not like NA vals
  D.all2.MLE.cov<- D.all.MLE.cov
  D.all2.MLE.cov[is.na(D.all2.MLE.cov)] <- 0
  D.all.MLE.cov.Sds <- colSds(D.all2.MLE.cov)
  D.all.MLE.cov.Sds[D.all.MLE.cov.Sds == 0] <- NA
  
  D.all.MCMC.Means <- colMeans(D.all.MCMC)
  # colSds for some reason does not like NA vals
  D.all2.MCMC <- D.all.MCMC
  D.all2.MCMC[is.na(D.all2.MCMC)] <- 0
  D.all.MCMC.Sds <- colSds(D.all2.MCMC)
  D.all.MCMC.Sds[D.all.MCMC.Sds == 0] <- NA
  
  D.all.MCMC.cov.Means <- colMeans(D.all.MCMC.cov)
  # colSds for some reason does not like NA vals
  D.all2.MCMC.cov <- D.all.MCMC.cov
  D.all2.MCMC.cov[is.na(D.all2.MCMC.cov)] <- 0
  D.all.MCMC.cov.Sds <- colSds(D.all2.MCMC.cov)
  D.all.MCMC.cov.Sds[D.all.MCMC.cov.Sds == 0] <- NA
  
  D.all.Means.mat <- matrix(D.all.Means, nrow = 4, ncol = 5, byrow = T) 
  D.all.Sds.mat <- matrix(D.all.Sds, nrow = 4, ncol = 5, byrow = T) 
  rownames(D.all.Means.mat) <- c("MLE (no covariates)","MLE (covariates)","MCMC (no covariates)","MCMC (covariates)")
  colnames(D.all.Means.mat) <- c("EEDE","REST","TTE","MCT","STE")
  rownames(D.all.Sds.mat) <- c("MLE (no covariates)","MLE (covariates)","MCMC (no covariates)","MCMC (covariates)")
  colnames(D.all.Sds.mat) <- c("EEDE","REST","TTE","MCT","STE")
  
  abm.props.Means <- colMeans(abm.props)
  EEDE.props.MLE.Means <- colMeans(EEDE.props.MLE)
  REST.props.MLE.Means <- colMeans(REST.props.MLE)
  TTE.props.MLE.Means <- colMeans(TTE.props.MLE)
  MCT.props.MLE.Means <- colMeans(MCT.props.MLE)
  EEDE.props.Means <- colMeans(EEDE.props)
  REST.props.Means <- colMeans(REST.props)
  TTE.props.Means <- colMeans(TTE.props)
  MCT.props.Means <- colMeans(MCT.props)
  
  abm.props.Sds <- colSds(abm.props)
  EEDE.props.MLE.Sds <- colSds(EEDE.props.MLE)
  REST.props.MLE.Sds <- colSds(REST.props.MLE)
  TTE.props.MLE.Sds <- colSds(TTE.props.MLE)
  MCT.props.MLE.Sds <- colSds(MCT.props.MLE)
  EEDE.props.Sds <- colSds(EEDE.props)
  REST.props.Sds <- colSds(REST.props)
  TTE.props.Sds <- colSds(TTE.props)
  MCT.props.Sds <- colSds(MCT.props)
  
  # all.props.Means <- rbind(abm.props.Means,EEDE.props.MLE.Means,REST.props.MLE.Means,
  #                          TTE.props.MLE.Means,MCT.props.MLE.Means,EEDE.props.Means,
  #                          REST.props.Means,TTE.props.Means,MCT.props.Means)
  # all.props.Sds <- rbind(abm.props.Sds,EEDE.props.MLE.Sds,REST.props.MLE.Sds,
  #                        TTE.props.MLE.Sds,MCT.props.MLE.Sds,EEDE.props.Sds,
  #                        REST.props.Sds,TTE.props.Sds,MCT.props.Sds)
  all.props.Means <- rbind(abm.props.Means,EEDE.props.MLE.Means,EEDE.props.Means,
                           REST.props.MLE.Means,REST.props.Means,
                           TTE.props.MLE.Means,TTE.props.Means,
                           MCT.props.MLE.Means,MCT.props.Means)
  all.props.Sds <- rbind(abm.props.Sds,EEDE.props.MLE.Sds,EEDE.props.Sds,
                         REST.props.MLE.Sds,REST.props.Sds,
                         TTE.props.MLE.Sds,TTE.props.Sds,
                         MCT.props.MLE.Sds,MCT.props.Sds)
  rownames(all.props.Means) <- c("ABM (Truth)",leg.props)
  colnames(all.props.Means) <- c("Slow","Medium","fast")
  rownames(all.props.Sds) <- c("ABM (Truth)",leg.props)
  colnames(all.props.Sds) <- c("Slow","Medium","fast")
  
  ####################################
  # Plot Results
  ####################################
  
  # # Save .csv files
  write.csv(D.all.Means.mat,paste(fig_dir,"sim_data/",means_label,".csv", sep = ""))
  write.csv(D.all.Sds.mat,paste(fig_dir,"sim_data/",means_label,"_sds.csv", sep = ""))
  write.csv(all.props.Means,paste(fig_dir,"sim_data/",props_label,".csv", sep = ""))
  write.csv(all.props.Sds,paste(fig_dir,"sim_data/",props_label,"_sds.csv", sep = ""))
  
  setEPS()
  postscript(paste(fig_dir,"figs/",means_label,".eps", sep = ""),width=8,height=5)
  op <- par(mar=c(5, 6, 4, 2) + 0.1)
  cam.props.label <- paste("Camera Bias: ", cam.dist.labels.caps[cam.dist.set], sep = "")
  plot(seq(0.7,4.7,by=1), 
       D.all.MLE.Means, ylim=c(min(nind,D.all.Means - D.all.Sds,na.rm=T),
       max(nind,D.all.Means + D.all.Sds,na.rm=T)), 
       xlim=c(0.5,5.5),xlab="Method", 
       ylab=paste(cam.props.label,"\n Mean Estimates"), pch=16, cex=1.8, xaxt = "n", cex.lab = 1.5, cex.axis = 1.3)
  # ylab="Mean Estimates", pch=16, cex=1.8, xaxt = "n", cex.lab = 1.5, cex.axis = 1.3)
points(seq(0.9,4.9,by=1),D.all.MLE.cov.Means, col="black", pch=1,cex=1.8)
  points(seq(1.1,5.1,by=1),D.all.MCMC.Means, col="blue", pch=16,cex=1.8)
  points(seq(1.3,5.3,by=1),D.all.MCMC.cov.Means, col="blue", pch=1,cex=1.8)
  arrows(x0=seq(0.7,4.7,by=1), y0=D.all.MLE.Means-D.all.MLE.Sds, 
         x1=seq(0.7,4.7,by=1), y1=D.all.MLE.Means+D.all.MLE.Sds, 
         code=3, angle=90, length=0.1, col="black", lwd=2)
  arrows(x0=seq(0.9,4.9,by=1), y0=D.all.MLE.cov.Means-D.all.MLE.cov.Sds, 
         x1=seq(0.9,4.9,by=1), y1=D.all.MLE.cov.Means+D.all.MLE.cov.Sds, 
         code=3, angle=90, length=0.1, col="black", lwd=2)
  arrows(x0=seq(1.1,5.1,by=1), y0=D.all.MCMC.Means-D.all.MCMC.Sds, 
         x1=seq(1.1,5.1,by=1), y1=D.all.MCMC.Means+D.all.MCMC.Sds, 
         code=3, angle=90, length=0.1, col="blue", lwd=2)
  arrows(x0=seq(1.3,5.3,by=1), y0=D.all.MCMC.cov.Means-D.all.MCMC.cov.Sds, 
         x1=seq(1.3,5.3,by=1), y1=D.all.MCMC.cov.Means+D.all.MCMC.cov.Sds, 
         code=3, angle=90, length=0.1, col="blue", lwd=2)
  lines(c(0.2,6.8), c(nind,nind), type = "l", lty = 2, lwd = 3.5)
  axis(1, at = c(1:5), labels = leg1, cex.axis = 1.3)
  legend("topright", c("MLE (no covariates)","MLE (covariates)","MCMC (no covariates)","MCMC (covariates)"),
         pch=c(16,1,16,1), col=c("black","black","blue","blue"))
  par(op)
  dev.off()
  
  setEPS()
  postscript(paste(fig_dir,"figs/",props_label,".eps", sep = ""),width=8,height=5)
  bar.p <- barplot(all.props.Means, 
                   names = c("Slow", "Medium", "Fast"), 
                   beside = T, 
                   legend.text=T,
                   # col = c("black",brewer.pal(4, "Set2"),brewer.pal(4, "Set2")),
                   col = c("black",brewer.pal(8, "Paired")),
                   # border = c("black",brewer.pal(8, "Paired")),
                   # alpha = c(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5),
                   # col=c("black","blue","red","green","violet","blue","red","green","violet"),
                   density=c(300,300,10,300,10,300,10,300,10),
                   angle=c(0,0,45,0,45,0,45,0,45),
                   ylim=c(0,1) , ylab="Proportion Occupied",xlab = "Movement Speeds",
                   # args.legend = list(x = "topright",
                   #                    inset = c(0.38, 0)),
                   cex=1.3, 
                   cex.lab = 1.5, 
                   cex.axis = 1.3)
  arrows(bar.p,all.props.Means+all.props.Sds, 
         bar.p, 
         all.props.Means-all.props.Sds, 
         angle=90, 
         code=3, 
         length=0.05)
  dev.off()
  
  
  # # Compare proportional staying time distributions to abm distributions 
  # comp.props <- rbind(abm.distribution.scale,stay.time.distribution.scale)
  # colnames(comp.props) <- c("Slow","Medium","Fast")
  # barplot(comp.props,legend = c("ABM Abundances","Staying Time"),
  #         args.legend = list(x = "topright"),
  #         beside=TRUE,
  #         xlab = "Landscape Type",
  #         # ylab = "Relative Quantities")
  #         ylab = "Proportional Quantities",
  #         cex=1.3, cex.lab = 1.5, cex.axis = 1.3)
  
  # # Compare absolute staying time distributions to abm distributions 
  # DF <- data.frame(Speed = c("Slow","Medium","Fast"), 
  #                  ABM = abm.distribution.scale, 
  #                  Stay_time = t.stay.props)
  # DFlong <- DF %>% 
  #   pivot_longer(cols = -Speed,names_to = "Type") %>%
  #   mutate(scaled_value=ifelse(Type=="ABM",value,value/10))
  # ggplot(DFlong,aes(x=Speed, y = scaled_value,fill= Type)) + 
  #   geom_col(position="dodge") + 
  #   scale_y_continuous(sec.axis = sec_axis(~ . * 10, name = "Staying Time"))+
  #   labs(y="ABM",x="Movement Speeds")
}
}
beep(sound=3)

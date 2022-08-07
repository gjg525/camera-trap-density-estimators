################################
# EEDE w/ covariates
################################
EEDE.cov.fn <- function(cam.counts, 
                        t.staying.dat, 
                        censor, 
                        param, 
                        Z, 
                        cam.samps,
                        t.steps){
  # Initialize with landscape-scale covariates
  beta <- exp(param[1])
  phi <- exp(Z%*%param[2:length(param)])
  
  u.all <- beta*phi/sum(phi)
  
  # repeat estimated parms for fitting
  u.all.cams <- u.all[cam.samps]
  phi.cams <- phi[cam.samps]
  phi.cam.rep <- matrix(rep(phi.cams,dim(t.staying.dat)[2]),
                        nrow = ncam,
                        ncol = dim(t.staying.dat)[2])
  
  # Account for censored times
  t.stay <- t.staying.dat
  t.stay[t.stay>=censor] <- NA
  stay.censor <- t.staying.dat
  stay.censor[stay.censor<censor] <- NA
  
  logL <- sum(dpois(cam.counts,u.all.cams,log=TRUE)) +
    sum(dexp(t.stay, 1/phi.cam.rep, log = TRUE),na.rm=T) +
    sum(pexp(stay.censor, 1/phi.cam.rep, lower.tail = F, log = TRUE),na.rm=T)
  return(logL)
}
################################
# REST no covariates
################################
REST.fn <- function(num.encounters.dat, 
                    t.staying.dat, 
                    censor,
                    t.steps,
                    cam.A,
                    param){
  # param: beta parameter for lambda
  u <- exp(param[1])
  phi <- exp(param[2])
  
  Y <- u*dt*t.steps*cam.A/phi

  # Account for censored times
  t.stay <- t.staying.dat
  t.stay[t.stay >= censor] <- NA
  stay.censor <- t.staying.dat
  stay.censor[stay.censor<censor] <- NA
  
  # logL <- sum(dpois(num.encounters.dat,Y,log=TRUE)) +
  logL <- sum(dpois(num.encounters.dat,Y,log=TRUE)) +
    sum(dexp(t.stay, 1/phi, log = TRUE),na.rm=T) +
    sum(pexp(stay.censor, 1/phi, lower.tail = F, log = TRUE),na.rm=T)
  # sum(dgamma(t.stay, phi, log = TRUE),na.rm=T) +
  # sum(pgamma(stay.censor, phi, lower.tail = F, log = TRUE),na.rm=T)
  return(logL)
}


################################
# REST w/ covariates
################################
REST.cov.fn <- function(num.encounters.dat, 
                        t.staying.dat, 
                        censor, 
                        param, 
                        Z, 
                        cov.inds,
                        t.steps){
  # Initialize with landscape-scale covariates
  num.covs <- sum(cov.inds)
  u <- exp(Z%*%param[1:num.covs])
  phi <- exp(Z%*%param[(num.covs+1):length(param)])

  eta <- u*dt*t.steps/phi
  
  # Account for censored times
  t.stay <- t.staying.dat
  t.stay[t.stay>=censor] <- NA
  t.stay.censor <- t.staying.dat
  t.stay.censor[t.stay.censor<censor] <- NA
  
  # match data dimensions
  eta.cams <- eta[cam.samps]
  phi.cams <- phi[cam.samps]
  phi.cam.rep <- matrix(rep(phi.cams,dim(t.staying.dat)[2]),nrow = ncam,ncol = dim(t.staying.dat)[2])
  
  logL <- sum(dpois(num.encounters.dat,eta.cams,log=TRUE)) +
    sum(dexp(t.stay, 1/phi.cam.rep, log = TRUE),na.rm=T) +
    sum(pexp(t.stay.censor, 1/phi.cam.rep, lower.tail = F, log = TRUE),na.rm=T)
    # sum(dgamma(t.stay, phi.cam.rep, log = TRUE),na.rm=T) +
    # sum(pgamma(t.stay.censor, phi.cam.rep, lower.tail = F, log = TRUE),na.rm=T)
  return(logL)
}


################################
# TTE no covariates
################################
TTE.fn <- function(TTE.dat,
                   t.staying.dat,
                   censor,
                   cam.A,
                   param){
  u <- exp(param[1])
  phi <- exp(param[2])
  eta <- u*cam.A/phi

  # Calculate censored time to event data
  TTE.dat[TTE.dat == censor] <- NA
  TTE.censor <- TTE.dat
  TTE.censor[] <- 0
  TTE.censor[which(is.na(TTE.dat))] <- censor

  # Calculate censored staying time data
  t.stay <- t.staying.dat
  t.stay[t.stay>=censor] <- NA
  stay.censor <- t.staying.dat
  stay.censor[stay.censor<censor] <- NA

  logL <- sum(dexp(TTE.dat, eta, log = TRUE),na.rm=T) +
    sum(pexp(TTE.censor, eta, lower.tail = F, log = TRUE)) +
    # sum(dgamma(t.stay, phi, log = TRUE),na.rm=T) +
    # sum(pgamma(stay.censor, phi, lower.tail = F, log = TRUE),na.rm=T)
    sum(dexp(t.stay, 1/phi, log = TRUE),na.rm=T) +
    sum(pexp(stay.censor, 1/phi, lower.tail = F, log = TRUE),na.rm=T)
}


################################
# TTE w/ covariates
################################
TTE.cov.fn <- function(TTE.dat, 
                       t.staying.dat, 
                       censor, 
                       param, 
                       Z, 
                       cam.samps, 
                       num.occ, 
                       cov.inds){
  # Initialize with landscape-scale covariates
  num.covs <- sum(cov.inds)
  u <- exp(Z%*%param[1:num.covs])
  phi <- exp(Z%*%param[(num.covs+1):length(param)])
  eta <- u/phi
  
  # repeat estimated parms to match data dimensions
  eta.cams <- eta[cam.samps]
  eta.cam.rep <- matrix(rep(eta.cams,num.occ),nrow = ncam,ncol = num.occ)
  phi.cams <- phi[cam.samps]
  phi.cam.rep <- matrix(rep(phi.cams,dim(t.staying.dat)[2]),nrow = ncam,ncol = dim(t.staying.dat)[2])
  
  # Define censored time to events
  TTE.dat[TTE.dat == censor] <- NA
  TTE.censor <- TTE.dat
  TTE.censor[] <- 0
  TTE.censor[which(is.na(TTE.dat))] <- censor
  
  # Account for censored staying times
  t.stay <- t.staying.dat
  t.stay[t.stay>=censor] <- NA
  stay.censor <- t.staying.dat
  stay.censor[stay.censor<censor] <- NA
  
  logL <- sum(dexp(TTE.dat, eta.cam.rep, log = TRUE),na.rm=TRUE) +
    sum(pexp(TTE.censor, eta.cam.rep, lower.tail = F, log = TRUE)) +
    sum(dexp(t.stay, 1/phi.cam.rep, log = TRUE),na.rm=T) +
    sum(pexp(stay.censor, 1/phi.cam.rep, lower.tail = F, log = TRUE),na.rm=T)
  # sum(dgamma(t.stay, phi.cam.rep, log = TRUE),na.rm=T) +
  # sum(pgamma(stay.censor, phi.cam.rep, lower.tail = F, log = TRUE),na.rm=T)
  return(logL)
}


################################
# MCT no covariates
################################
MCT.fn <- function(cam.counts, param){
  # param: beta parameter for lambda
  d <- exp(param[1])
  
  logL <- sum(dpois(cam.counts,d,log=TRUE))
  return(logL)
}


################################
# MCT w/ covariates
################################
MCT.cov.fn <- function(cam.counts, param, Z, cam.samps, cov.inds){
  # Initialize with landscape-scale covariates
  num.covs <- sum(cov.inds)
  d <- exp(Z%*%param[1:num.covs])

  # repeat estimated parms for fitting
  d.cams <- d[cam.samps]
  d.cam.rep <- matrix(rep(d.cams,t.steps),nrow = ncam,ncol = t.steps)

  logL <- sum(dpois(cam.counts,d.cam.rep,log=TRUE)) 
  return(logL)
}






####################################
# STE no covariates
####################################
STE.fn <- function(STE.dat, censor, param){
  # param: beta parameter for lambda
  u <- exp(param)
  
  # Define censored space to event
  STE.dat[STE.dat == censor] <- NA
  STE.censor <- STE.dat
  STE.censor[] <- 0
  STE.censor[which(is.na(STE.dat))] <- censor
  
  logL <- sum(dexp(STE.dat, u, log = TRUE),na.rm=TRUE) +
    sum(pexp(STE.censor, u, lower.tail = F, log = TRUE),na.rm=TRUE)
  return(logL)
}




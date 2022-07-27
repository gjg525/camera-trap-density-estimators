## NOTE tidyr may cause errors in extract function
library(tidyverse)
library(pbapply) # For progress bar in bootstrap
library(msm) # Delta method from Fewster
library(MASS)
library(Directional)
library(raster)
library(fields)
library(RColorBrewer)
library(rasterVis)
library(truncnorm)
library(mvtnorm)
library(gridExtra)
library(fBasics)
library(coda)
library(Rlab)
library(fitdistrplus)
library(patchwork)
library(crayon)
library(beepr)
options(error = beep)

############################################
# Define functions
############################################

########################################
## Collect REST and TTE data
########################################
collect.TTE.data <- function(cam.counts,JJ,ncam,num.occ,t.steps) {
  num.encounters.dat <- matrix(NA,ncam,1)
  TTE.dat <- matrix(NA,ncam,num.occ)
  t.staying.dat <- matrix(NA,ncam,1)  # Extend matrix as number of encounters increase
  # # Calculate number of first encounters and staying time from simulation (REST, TTE)
  for(cc in 1:ncam){
    t.stay.jj <- c()
    num.encounters.temp <- c()
    
    # Determine if a new individual has entered or exited the camera
    all.counts <- cam.counts[cc,]  # number of counts in sample occasion
    count.diff <- all.counts[2:t.steps]- all.counts[1:(t.steps-1)] 
    
    # find enter and exit times for each individual
    t.enter <- c()
    t.exit <- c()
    if(max(count.diff)>0){
      for(ent.id in 1:max(count.diff)){
        t.enter <- c(t.enter,rep(which(count.diff == ent.id),ent.id)) # individual enters camera
      }
    }
    if(min(count.diff)<0){
      for(ex.id in -1:min(count.diff)){
        t.exit <- c(t.exit, rep(which(count.diff == ex.id),abs(ex.id))) # individual leaves camera
      }
    }
    # Account for individuals already in camera at start and finish times
    if(all.counts[1]>0){
      t.enter <- c(rep(0,all.counts[1]),t.enter)
    }
    if(all.counts[t.steps]>0){
      t.exit <- c(t.exit,rep(t.steps,all.counts[t.steps]))
    }
    t.enter <- sort(t.enter)
    t.exit <- sort(t.exit)
    
    # Number of encounters (REST)
    # num.encounters.temp <- c(num.encounters.temp,length(t.enter))
    t.stay.jj <- c(t.stay.jj,t.exit-t.enter)
    
    # Calculate number of encounters for each camera (REST)
    num.encounters.dat[cc] <- length(t.enter)
    
    for(jj in 1:num.occ){
      occ.inds <- (JJ*(jj-1)+1):(JJ*jj)      # sampling occasion indices
      occ.counts <- cam.counts[cc,occ.inds]  # number of counts in sample occasion
      
      # Calculate TTE
      if(all(occ.counts==0)){
        TTE.dat[cc,jj] <- NA
      }else{
        first.enc <- which(occ.counts>0)
        TTE.dat[cc,jj] <- first.enc[1]
      }
    }
    # Calculate staying times for each encounter (REST and TTE)
    if(sum(num.encounters.dat[cc])>0){
      # Extend matrix by number of new encounters
      if(sum(num.encounters.dat[cc])>dim(t.staying.dat)[2]){
        t.staying.dat <- cbind(t.staying.dat, matrix(NA,ncam,sum(num.encounters.dat[cc,])-1))
      }
      t.staying.dat[cc,1:sum(num.encounters.dat[cc])] <- t.stay.jj
    }
  }
  list(num.encounters.dat,TTE.dat,t.staying.dat)
}

########################################
## Collect STE data
########################################
collect.STE.data <- function(cam.counts,CC,ncam,cam.occ,t.steps) {
  STE.dat <- matrix(0,cam.occ,t.steps)
  samp.inds <- matrix(0,ncam,t.steps)   # randomized camera indices

  # Calculate space to encounters for IS and STE
  for(tt in 1:t.steps){
    # Randomize camera samples
    samp.inds[,tt] <- sample(1:ncam)
    for(cc in 1:cam.occ){
      STE.counts <- cam.counts[samp.inds[,tt],tt]
      if(all(STE.counts==0)){
        STE.dat[cc,tt] <- NA
      }else{
        first.enc <- which(STE.counts>0)
        STE.dat[cc,tt] <- first.enc[1]
      }
    }
  }
  list(STE.dat,samp.inds)
}   

########################################
## Collect REST and TTE data
########################################
num.encounters.dat <- matrix(NA,ncam,1)
TTE.dat <- matrix(NA,ncam,num.occ)
t.staying.dat <- matrix(NA,ncam,1)  # Extend matrix as number of encounters increase

# # Calculate number of first encounters and staying time from simulation (REST, TTE)
for(cc in 1:ncam){
  t.stay.jj <- c()
  num.encounters.temp <- c()
  
  cam_ind <- cam.samps[cc]
  
  cam_dat_temp <- filter(animalxy.all, XY_inds == 11) |> 
                  group_by(ID) |> 
                  summarise(inds = min(ii):max(ii),
                            i_diff = animalxy.all$ii[(inds[1]+1):inds[2]] - 
                            animalxy.all$ii[inds[1]:(inds[2]-1)],
                            t_diff = animalxy.all$t[(inds[1]+1):inds[2]] - 
                              animalxy.all$t[inds[1]:(inds[2]-1)],
                            .groups = 'drop')
  cam_dat_temp <- group_by(animalxy.all, ID) |> 
    summarise(t_diff =  t[2:nrow(animalxy.all)] - 
                t[1:(nrow(animalxy.all)-1)],
              XY_inds = XY_inds[1:(nrow(animalxy.all)-1)],
              .groups = 'drop'
                )
  
  num.encounters.dat[cc] <- sum(cam_dat_temp$i_diff>1) + length(unique(cam_dat_temp$ID))
  
  # Determine if a new individual has entered or exited the camera
  all.counts <- cam.counts[cc,]  # number of counts in sample occasion
  count.diff <- all.counts[2:t.steps]- all.counts[1:(t.steps-1)] 
  
  # Find enter and exit times for each individual
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
  
  # Staying time for all encounters
  t.stay.jj <- c(t.stay.jj,t.exit-t.enter)
  
  # Calculate number of encounters for each camera (REST)
  num.encounters.dat[cc] <- length(t.enter)
  
  # Collect TTE
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

########################################
## Collect STE data
########################################
  
STE.dat <- matrix(0,1,t.steps)
samp.inds <- matrix(0,ncam,t.steps)   # randomized camera indices

# Calculate space to encounters for STE
for(tt in 1:t.steps){
  # Randomize camera samples
  samp.inds[,tt] <- sample(1:ncam)
  STE.counts <- cam.counts[samp.inds[,tt],tt]
  if(all(STE.counts==0)){
    STE.dat[tt] <- NA
  }else{
    first.enc <- which(STE.counts>0)
    STE.dat[tt] <- first.enc[1]*cam.A
  }
}

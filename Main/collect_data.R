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
  
  cam_dat_temp <- group_by(animalxy.all, ID) %>%  
              summarise(t_diff =  t[2:length(t)] - t[1:(length(t)-1)],
              XY_inds = XY_inds[1:(length(XY_inds)-1)],
              i_diff = ii[2:length(ii)] - ii[1:(length(ii)-1)],
              ii = ii[1:(length(ii)-1)],
              .groups = 'drop'
                ) 
  if(sum(cam_dat_temp$XY_inds == cam_ind) == 0) {
    num.encounters.dat[cc] <- 0
    TTE.dat[cc, ] <- NA
  } else {
    t_staying_temp <- dplyr::filter(cam_dat_temp, XY_inds == cam_ind) %>% 
      group_by(ID) %>% 
      mutate(ii_diff = ifelse(length(ii) == 1,
                               1,
                               c(1,ii[2:length(ii)] - ii[1:(length(ii)-1)])),
             encounter = cumsum(c(1, abs(ii[-length(ii)] - ii[-1]) > 1)))
    
    num.encounters.dat[cc] <- nrow(group_by(t_staying_temp, ID) %>% 
                                     summarise(mm = unique(encounter),
                                   .groups = 'drop'))
  
    t_staying <- group_by(t_staying_temp, ID, encounter) %>% 
      summarise(t_stay = sum(t_diff),
                .groups = 'drop')
      
    # Extend matrix by number of new encounters
    if(sum(num.encounters.dat[cc])>dim(t.staying.dat)[2]){
      t.staying.dat <- cbind(t.staying.dat, matrix(NA,ncam,sum(num.encounters.dat[cc,])-1))
    }
    t.staying.dat[cc,1:sum(num.encounters.dat[cc])] <- t_staying$t_stay
    
    # Calculate TTE
    for(jj in 1:num.occ){
      occ.times <- t.steps*c((jj-1)/num.occ, jj/num.occ)
      TTE_temp <- dplyr::filter(animalxy.all, 
                         XY_inds == cam_ind, 
                         t > occ.times[1], 
                         t<= occ.times[2]) 
    
    if(nrow(TTE_temp) == 0){
      TTE.dat[cc,jj] <- NA
    }else{
      TTE.dat[cc,jj] <- min(TTE_temp$t) - (jj-1)*JJ
    }
    }
  }
}
TTE.dat[TTE.dat == 1] <- 0

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

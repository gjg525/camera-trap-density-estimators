# Metropolis-within-Gibbs sampler
########################################
# MCMC for TDST w/ covariates
########################################
fit.model.mcmc.TDST.cov <- function(n.iter,
                                    gamma.start,
                                    kappa.start,
                                    gamma.prior.var,
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
        mh1 <- sum(dexp(t.staying.dat.all, 1/phi.star.cam.rep,log=TRUE),na.rm=TRUE) +
          sum(pexp(t.staying.dat.censor, 1/phi.star.cam.rep, lower.tail = F, log = TRUE),na.rm=TRUE) +
          sum(dnorm(kappa.star,0,kappa.prior.var^0.5,log=TRUE))
        mh2 <- sum(dexp(t.staying.dat.all, 1/phi.all.cam.rep,log=TRUE),na.rm=TRUE) +
          sum(pexp(t.staying.dat.censor, 1/phi.all.cam.rep, lower.tail = F, log = TRUE),na.rm=TRUE) +
          sum(dnorm(kappa[i,],0,kappa.prior.var^0.5,log=TRUE))
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


########################################
# MCMC for REST no covariates
########################################
fit.model.mcmc.REST <- function(n.iter,
                                gamma.start,
                                kappa.start,
                                gamma.prior.var,
                                kappa.prior.var,
                                gamma.tune,
                                kappa.tune,
                                num.encounters.dat,
                                t.staying.dat,
                                cam.A,
                                t.steps,
                                censor){

  # Variables that will be saved
  gamma <- matrix(,n.iter+1,1)
  kappa <- matrix(,n.iter+1,1)
  tot.u <- matrix(,n.iter+1,1)
  accept <- matrix(,n.iter+1,2)

  gamma[1] <- gamma.start
  colnames(gamma) <- "gamma"
  kappa[1] <- kappa.start
  colnames(kappa) <- "kappa"
  colnames(tot.u) <- "Total estimate"
  colnames(accept) <- c("accept.rate.gamma","accept.rate.kappa")

  # Account for censored times
  t.staying.dat.all <- t.staying.dat
  t.staying.dat.all[t.staying.dat.all>censor] <- NA
  t.staying.dat.censor <- t.staying.dat
  t.staying.dat.censor[t.staying.dat.censor<=censor] <- NA

  u <- exp(gamma[1])
  phi <- exp(kappa[1])
  eta <- u*dt*t.steps*cam.A/phi
  tot.u <- tot.A*u

  tune.check <- 100
  batch_n <- 0

  # Begin MCMC loop
  # prog_bar <- txtProgressBar(min = 0,max = n.iter,style = 3,width = 50,char = "=")
  for(i in 1:n.iter){
    # setTxtProgressBar(prog_bar, i)
    #Sample gamma
    gamma.star <- rnorm(1,gamma[i],exp(2*gamma.tune))
    u.star <- exp(gamma.star)
    eta.star <- u.star*dt*t.steps*cam.A/phi

    if(all(eta.star>0)){
      mh1 <- sum(dpois(num.encounters.dat,eta.star,log=TRUE),na.rm=TRUE) +
        sum(dnorm(gamma.star,0,gamma.prior.var^0.5,log=TRUE))
      mh2 <- sum(dpois(num.encounters.dat,eta,log=TRUE),na.rm=TRUE) +
        sum(dnorm(gamma[i],0,gamma.prior.var^0.5,log=TRUE))
      mh <- exp(mh1-mh2)

      if(mh>runif(1)){
        gamma[i+1] <- gamma.star;
        accept[i+1,1] <- 1;
        eta <- eta.star}
      else{
        gamma[i+1] <- gamma[i];
        accept[i+1,1] <- 0}
    } else{
      gamma[i+1] <- gamma[i];
      accept[i+1,1] <- 0}

    #Sample kappa
    kappa.star <- rnorm(1,kappa[i],exp(2*kappa.tune))
    phi.star <- exp(kappa.star)

    if(1/phi.star>0 & !is.infinite(1/phi.star)){
      mh1 <- sum(dexp(t.staying.dat.all, 1/phi.star,log=TRUE),na.rm=TRUE) +
        sum(pexp(t.staying.dat.censor, 1/phi.star, lower.tail = F, log = TRUE),na.rm=TRUE) +
        sum(dnorm(kappa.star,0,kappa.prior.var^0.5,log=TRUE))
      mh2 <- sum(dexp(t.staying.dat.all, 1/phi,log=TRUE),na.rm=TRUE) +
        sum(pexp(t.staying.dat.censor, 1/phi, lower.tail = F, log = TRUE),na.rm=TRUE) +
        sum(dnorm(kappa[i],0,kappa.prior.var^0.5,log=TRUE))
      mh <- exp(mh1-mh2)

      if(mh>runif(1)){
        kappa[i+1] <- kappa.star;
        accept[i+1,2] <- 1;
        phi <- phi.star} else{
          kappa[i+1] <- kappa[i];
          accept[i+1,2] <- 0}
    } else{
      kappa[i+1] <- kappa[i];
      accept[i+1,2] <- 0}

    u <- exp(gamma[i+1])
    phi <- exp(kappa[i+1])
    eta <- u*dt*t.steps*cam.A/phi
    tot.u[i+1] <- tot.A*u

    # Update tuning parms
    if(i%%tune.check == 0){
      batch_n <- batch_n+1
      delta_n <- batch_n^-1
      # delta_n <- min(0.01,batch_n^-1)
      accept.gamma.check <- mean(accept[(i-tune.check+1):i,1],na.rm=T)
      gamma.tune[accept.gamma.check>0.44] <- gamma.tune[accept.gamma.check>0.44] + delta_n
      gamma.tune[accept.gamma.check<=0.44] <- gamma.tune[accept.gamma.check<=0.44] - delta_n
      accept.kappa.check <- mean(accept[(i-tune.check+1):i,2],na.rm=T)
      kappa.tune[accept.kappa.check>0.44] <- kappa.tune[accept.kappa.check>0.44] + delta_n
      kappa.tune[accept.kappa.check<=0.44] <- kappa.tune[accept.kappa.check<=0.44] - delta_n
    }

  }
  # print("MCMC complete")

  list(accept = accept,gamma = gamma,kappa = kappa,tot.u = tot.u)
}


########################################
# MCMC for REST w/ covariates
########################################
fit.model.mcmc.REST.cov <- function(n.iter,
                                    gamma.start,
                                    kappa.start,
                                    gamma.prior.var,
                                    kappa.prior.var,
                                    gamma.tune,
                                    kappa.tune,
                                    num.encounters.dat,
                                    t.staying.dat,
                                    covariate_labels,
                                    covariates.index,
                                    cam.A,
                                    cell.A,
                                    t.steps,
                                    censor){

  # Variables that will be saved
  gamma <- matrix(,n.iter+1,sum(covariates.index==1))
  kappa <- matrix(,n.iter+1,sum(covariates.index==1))
  tot.u <- matrix(,n.iter+1,1)
  accept <- matrix(,n.iter+1,2*sum(covariates.index==1))

  gamma[1,] <- gamma.start
  colnames(gamma) <- paste0("gamma.", covariate_labels)
  kappa[1,] <- kappa.start
  colnames(kappa) <- paste0("kappa.", covariate_labels)
  colnames(tot.u) <- "Total estimate"
  colnames(accept) <- c(paste0("accept.rate.gamma.", covariate_labels),
                        paste0("accept.rate.kappa.", covariate_labels))

  # Account for censored times
  t.staying.dat.all <- t.staying.dat
  t.staying.dat.all[t.staying.dat.all>censor] <- NA
  t.staying.dat.censor <- t.staying.dat
  t.staying.dat.censor[t.staying.dat.censor<=censor] <- NA

  # Initialize with landscape-scale covariates
  u <- exp(Z%*%gamma[1,])
  phi <- exp(Z%*%kappa[1,])
  eta <- u*t.steps/phi*cell.A
  tot.u[1] <- sum(u)*cell.A

  tune.check <- 100
  batch_n <- 0

  # Begin MCMC loop
  # prog_bar <- txtProgressBar(min = 0,max = n.iter,style = 3,width = 50,char = "=")
  for(i in 1:n.iter){
    # setTxtProgressBar(prog_bar, i)
    #Sample gamma
    for(gg in 1:sum(covariates.index==1)){
      gamma.star <- gamma[i,]
      gamma.star[gg] <- rnorm(1,gamma[i,gg],exp(2*gamma.tune[gg]))
      u.star <- exp(Z%*%gamma.star)
      eta.star <- u.star*t.steps/phi*cell.A

      eta.star.cams <- eta.star[cam.samps]
      eta.all.cams <- eta[cam.samps]

      if(all(eta.star.cams>0)){
        mh1 <- sum(dpois(num.encounters.dat,eta.star.cams,log=TRUE),na.rm=TRUE) +
          sum(dnorm(gamma.star,0,gamma.prior.var^0.5,log=TRUE))
        mh2 <- sum(dpois(num.encounters.dat,eta.all.cams,log=TRUE),na.rm=TRUE) +
          sum(dnorm(gamma[i,],0,gamma.prior.var^0.5,log=TRUE))
        mh <- exp(mh1-mh2)

        if(mh>runif(1) && !is.nan(mh)){
          gamma[i,] <- gamma.star;
          accept[i+1,gg] <- 1;
          eta <- eta.star}
        else{
          gamma[i,] <- gamma[i,];
          accept[i+1,gg] <- 0}
      } else{
        gamma[i,] <- gamma[i,];
        accept[i+1,gg] <- 0}
    }
    gamma[i+1,] <- gamma[i,];
    u <- exp(Z%*%gamma[i+1,])
    eta <- u*t.steps/phi*cell.A

    #Sample kappa
    for(kk in 1:sum(covariates.index==1)){
      kappa.star <- kappa[i,]
      kappa.star[kk] <- rnorm(1,kappa[i,kk],exp(2*kappa.tune[kk]))
      phi.star <- exp(Z%*%kappa.star)

      # repeat estimated parms for fitting
      phi.star.cams <- phi.star[cam.samps] * cam.A / cell.A
      phi.star.cam.rep <- matrix(rep(phi.star.cams,dim(t.staying.dat.all)[2]),
                                 nrow = ncam,ncol = dim(t.staying.dat.all)[2])
      phi.all.cams <- phi[cam.samps] * cam.A / cell.A
      phi.all.cam.rep <- matrix(rep(phi.all.cams,dim(t.staying.dat.all)[2]),
                                nrow = ncam,ncol = dim(t.staying.dat.all)[2])

      # if(all(phi.star.cams>0)){
      #   mh1 <- sum(dgamma(t.staying.dat.all,phi.star.cam.rep,log=TRUE),na.rm=TRUE) +
      #     sum(pgamma(t.staying.dat.censor, phi.star.cam.rep, lower.tail = F, log = TRUE),na.rm=TRUE) +
      #     sum(dnorm(kappa.star,0,kappa.prior.var^0.5,log=TRUE))
      #   mh2 <- sum(dgamma(t.staying.dat.all,phi.all.cam.rep,log=TRUE),na.rm=TRUE) +
      #     sum(pgamma(t.staying.dat.censor, phi.all.cam.rep, lower.tail = F, log = TRUE),na.rm=TRUE) +
      #     sum(dnorm(kappa[i,],0,kappa.prior.var^0.5,log=TRUE))
      if(all(1/phi.star.cams>0) & all(!is.infinite(1/phi.star.cams))){
        mh1 <- sum(dexp(t.staying.dat.all, 1/phi.star.cam.rep,log=TRUE),na.rm=TRUE) +
          sum(pexp(t.staying.dat.censor, 1/phi.star.cam.rep, lower.tail = F, log = TRUE),na.rm=TRUE) +
          sum(dnorm(kappa.star,0,kappa.prior.var^0.5,log=TRUE))
        mh2 <- sum(dexp(t.staying.dat.all, 1/phi.all.cam.rep,log=TRUE),na.rm=TRUE) +
          sum(pexp(t.staying.dat.censor, 1/phi.all.cam.rep, lower.tail = F, log = TRUE),na.rm=TRUE) +
          sum(dnorm(kappa[i,],0,kappa.prior.var^0.5,log=TRUE))
        mh <- exp(mh1-mh2)

        if(mh>runif(1)){
          kappa[i,] <- kappa.star;
          accept[i+1,sum(covariates.index==1)+kk] <- 1;
          phi <- phi.star}
        else{
          kappa[i,] <- kappa[i,];
          accept[i+1,sum(covariates.index==1)+kk] <- 0}
      } else{
        kappa[i,] <- kappa[i,];
        accept[i+1,sum(covariates.index==1)+kk] <- 0}
    }
    kappa[i+1,] <- kappa[i,];
    phi <- exp(Z%*%kappa[i+1,])

    tot.u[i+1] <- sum(u)*cell.A

    # Update tuning parms
    if(i%%tune.check == 0){
      batch_n <- batch_n+1
      delta_n <- batch_n^-1
      # delta_n <- min(0.01,batch_n^-1)
      accept.gamma.check <- colMeans(accept[(i-tune.check+1):i,1:sum(covariates.index==1)],na.rm=T)
      gamma.tune[which(accept.gamma.check>0.44)] <- gamma.tune[which(accept.gamma.check>0.44)] + delta_n
      gamma.tune[which(accept.gamma.check<=0.44)] <- gamma.tune[which(accept.gamma.check<=0.44)] - delta_n
      accept.kappa.check <- colMeans(accept[(i-tune.check+1):i,(sum(covariates.index==1)+1):(2*sum(covariates.index==1))],na.rm=T)
      kappa.tune[which(accept.kappa.check>0.44)] <- kappa.tune[which(accept.kappa.check>0.44)] + delta_n
      kappa.tune[which(accept.kappa.check<=0.44)] <- kappa.tune[which(accept.kappa.check<=0.44)] - delta_n
    }

  }
  # print("MCMC complete")

  list(accept = accept,gamma = gamma,kappa = kappa,tot.u = tot.u,u = u)
}


########################################
# MCMC for TTE no covariates
########################################
fit.model.mcmc.TTE <- function(n.iter,
                               gamma.start,
                               kappa.start,
                               gamma.prior.var,
                               kappa.prior.var,
                               gamma.tune,
                               kappa.tune,
                               TTE.dat,
                               t.staying.dat,
                               cam.A,
                               JJ,
                               censor){

  # Variables that will be saved
  gamma <- matrix(,n.iter+1,1)
  kappa <- matrix(,n.iter+1,1)
  tot.u <- matrix(,n.iter+1,1)
  accept <- matrix(,n.iter+1,2)

  gamma[1] <- gamma.start
  colnames(gamma) <- "gamma"
  kappa[1] <- kappa.start
  colnames(kappa) <- "kappa"
  colnames(tot.u) <- "Total estimate"
  colnames(accept) <- c("accept.rate.gamma","accept.rate.kappa")

  # Define censored time to events
  TTE.dat[TTE.dat == censor] <- NA
  data.censor <- TTE.dat
  data.censor[] <- 0
  data.censor[which(is.na(TTE.dat))] <- censor

  # Account for censored times
  t.staying.dat.all <- t.staying.dat
  t.staying.dat.all[t.staying.dat.all>censor] <- NA
  t.staying.dat.censor <- t.staying.dat
  t.staying.dat.censor[t.staying.dat.censor<=censor] <- NA

  u <- exp(gamma[1])
  phi <- exp(kappa[1])
  tot.u <- u*tot.A
  eta <- u*cam.A/phi

  tune.check <- 100
  batch_n <- 0

  # Begin MCMC loop
  # prog_bar <- txtProgressBar(min = 0,max = n.iter,style = 3,width = 50,char = "=")
  for(i in 1:n.iter){
    # setTxtProgressBar(prog_bar, i)
    #Sample gamma
    gamma.star <- rnorm(1,gamma[i],exp(2*gamma.tune))
    u.star <- exp(gamma.star)
    eta.star <- u.star*cam.A/phi

    if(all(eta.star>0)){
      mh1 <- sum(dexp(TTE.dat,eta.star,log=TRUE),na.rm=TRUE) +
        sum(pexp(data.censor, eta.star, lower.tail = F, log = TRUE)) +
        sum(dnorm(gamma.star,0,gamma.prior.var^0.5,log=TRUE))
      mh2 <- sum(dexp(TTE.dat,eta,log=TRUE),na.rm=TRUE) +
        sum(pexp(data.censor, eta, lower.tail = F, log = TRUE)) +
        sum(dnorm(gamma[i],0,gamma.prior.var^0.5,log=TRUE))
      mh <- exp(mh1-mh2)

      if(mh>runif(1)){gamma[i+1] <- gamma.star;
      accept[i+1,1] <- 1;
      eta <- eta.star} else{
        gamma[i+1] <- gamma[i];
        accept[i+1,1] <- 0}
    } else{
      gamma[i+1] <- gamma[i];
      accept[i+1,1] <- 0}

    #Sample kappa
    kappa.star <- rnorm(1,kappa[i],exp(2*kappa.tune))
    phi.star <- exp(kappa.star)

    if(1/phi.star>0 & !is.infinite(1/phi.star)){
      mh1 <- sum(dexp(t.staying.dat.all, 1/phi.star,log=TRUE),na.rm=TRUE) +
        sum(pexp(t.staying.dat.censor, 1/phi.star, lower.tail = F, log = TRUE),na.rm=TRUE) +
        sum(dnorm(kappa.star,0,kappa.prior.var^0.5,log=TRUE))
      mh2 <- sum(dexp(t.staying.dat.all, 1/phi,log=TRUE),na.rm=TRUE) +
        sum(pexp(t.staying.dat.censor, 1/phi, lower.tail = F, log = TRUE),na.rm=TRUE) +
        sum(dnorm(kappa[i],0,kappa.prior.var^0.5,log=TRUE))
      mh <- exp(mh1-mh2)

      if(mh>runif(1)){kappa[i+1] <- kappa.star;
      accept[i+1,2] <- 1;
      phi <- phi.star} else{
        kappa[i+1] <- kappa[i];
        accept[i+1,2] <- 0}
    } else{
      kappa[i+1] <- kappa[i];
      accept[i+1,2] <- 0}

    u <- exp(gamma[i+1])
    phi <- exp(kappa[i+1])
    tot.u[i+1] <-u*tot.A
    eta <- u*cam.A/phi

    # Update tuning parms
    if(i%%tune.check == 0){
      batch_n <- batch_n+1
      delta_n <- batch_n^-1
      # delta_n <- min(0.01,batch_n^-1)
      accept.gamma.check <- mean(accept[(i-tune.check+1):i,1],na.rm=T)
      gamma.tune[accept.gamma.check>0.44] <- gamma.tune[accept.gamma.check>0.44] + delta_n
      gamma.tune[accept.gamma.check<=0.44] <- gamma.tune[accept.gamma.check<=0.44] - delta_n
      accept.kappa.check <- mean(accept[(i-tune.check+1):i,2],na.rm=T)
      kappa.tune[accept.kappa.check>0.44] <- kappa.tune[accept.kappa.check>0.44] + delta_n
      kappa.tune[accept.kappa.check<=0.44] <- kappa.tune[accept.kappa.check<=0.44] - delta_n
    }

  }
  # print("MCMC complete")

  list(accept = accept,gamma = gamma,kappa = kappa,tot.u = tot.u)
}


########################################
# MCMC for TTE w/ covariates
########################################
fit.model.mcmc.TTE.cov <- function(n.iter,
                                   gamma.start,
                                   kappa.start,
                                   gamma.prior.var,
                                   kappa.prior.var,
                                   gamma.tune,
                                   kappa.tune,
                                   TTE.dat,
                                   t.staying.dat,
                                   censor,
                                   covariate_labels,
                                   covariates.index,
                                   cam.A,
                                   cell.A,
                                   JJ){

  # Variables that will be saved
  gamma <- matrix(,n.iter+1,sum(covariates.index==1))
  kappa <- matrix(,n.iter+1,sum(covariates.index==1))
  tot.u <- matrix(,n.iter+1,1)
  accept <- matrix(,n.iter+1,2*sum(covariates.index==1))
  gamma[1,] <- gamma.start
  colnames(gamma) <- paste0("gamma.",covariate_labels)
  kappa[1,] <- kappa.start
  colnames(kappa) <- paste0("kappa.",covariate_labels)
  colnames(tot.u) <- "Total estimate"
  colnames(accept) <- c(paste0("accept.rate.gamma.",covariate_labels),
                        paste0("accept.rate.kappa.",covariate_labels))

  # Define censored time to events
  TTE.dat[TTE.dat == censor] <- NA
  data.censor <- TTE.dat
  data.censor[] <- 0
  data.censor[which(is.na(TTE.dat))] <- censor

  # Account for censored staying times
  t.staying.dat.all <- t.staying.dat
  t.staying.dat.all[t.staying.dat.all>censor] <- NA
  t.staying.dat.censor <- t.staying.dat
  t.staying.dat.censor[t.staying.dat.censor<=censor] <- NA

  # Initialize with landscape-scale covariates
  u <- exp(Z%*%gamma[1,])
  phi <- exp(Z%*%kappa[1,])
  tot.u[1] <- sum(u)*cell.A
  eta <- u/phi*cell.A

  tune.check <- 100
  batch_n <- 0

  # Begin MCMC loop
  # prog_bar <- txtProgressBar(min = 0,max = n.iter,style = 3,width = 50,char = "=")
  for(i in 1:n.iter){
    # setTxtProgressBar(prog_bar, i)
    #Sample gamma
    for(gg in 1:sum(covariates.index==1)){
      gamma.star <- gamma[i,]
      gamma.star[gg] <- rnorm(1,gamma[i,gg],exp(2*gamma.tune[gg]))
      u.star <- exp(Z%*%gamma.star)
      eta.star <- u.star/phi*cell.A

      # repeat estimated parms for fitting
      eta.star.cams <- eta.star[cam.samps]
      eta.star.cam.rep <- matrix(rep(eta.star.cams,num.occ),nrow = ncam,ncol = num.occ)
      eta.all.cams <- eta[cam.samps]
      eta.all.cam.rep <- matrix(rep(eta.all.cams,num.occ),nrow = ncam,ncol = num.occ)

      if(all(c(eta.star.cam.rep,eta.all.cam.rep)>1e-10)){
        mh1 <- sum(dexp(TTE.dat,eta.star.cam.rep,log=TRUE),na.rm=TRUE) +
          sum(pexp(data.censor, eta.star.cam.rep, lower.tail = F, log = TRUE)) +
          sum(dnorm(gamma.star,0,gamma.prior.var^0.5,log=TRUE))
        mh2 <- sum(dexp(TTE.dat,eta.all.cam.rep,log=TRUE),na.rm=TRUE) +
          sum(pexp(data.censor, eta.all.cam.rep, lower.tail = F, log = TRUE)) +
          sum(dnorm(gamma[i,],0,gamma.prior.var^0.5,log=TRUE))
        mh <- exp(mh1-mh2)

        if(mh>runif(1)){
          gamma[i,] <- gamma.star;
          accept[i+1,gg] <- 1;
          u<- u.star;
          eta <- eta.star} else{
            gamma[i,] <- gamma[i,];
            accept[i+1,gg] <- 0}
      } else{
        gamma[i,] <- gamma[i,];
        accept[i+1,gg] <- 0}
    }
    gamma[i+1,] <- gamma[i,];
    u <- exp(Z%*%gamma[i+1,])
    eta <- u/phi*cell.A

    #Sample kappa
    for(kk in 1:sum(covariates.index==1)){
      kappa.star <- kappa[i,]
      kappa.star[kk] <- rnorm(1,kappa[i,kk],exp(2*kappa.tune[kk]))
      phi.star <- exp(Z%*%kappa.star)

      # repeat estimated parms for fitting
      phi.star.cams <- phi.star[cam.samps] * cam.A / cell.A
      phi.star.cam.rep <- matrix(rep(phi.star.cams,dim(t.staying.dat.all)[2]),nrow = ncam,ncol = dim(t.staying.dat.all)[2])
      phi.all.cams <- phi[cam.samps] * cam.A / cell.A
      phi.all.cam.rep <- matrix(rep(phi.all.cams,dim(t.staying.dat.all)[2]),nrow = ncam,ncol = dim(t.staying.dat.all)[2])

      # if(all(phi.star.cams>0)){
      # mh1 <- sum(dgamma(t.staying.dat.all,phi.star.cam.rep,log=TRUE),na.rm=TRUE) +
      #   sum(pgamma(t.staying.dat.censor, phi.star.cam.rep, lower.tail = F, log = TRUE),na.rm=TRUE) +
      #   sum(dnorm(kappa.star,0,kappa.prior.var^0.5,log=TRUE))
      # mh2 <- sum(dgamma(t.staying.dat.all,phi.all.cam.rep,log=TRUE),na.rm=TRUE) +
      #   sum(pgamma(t.staying.dat.censor, phi.all.cam.rep, lower.tail = F, log = TRUE),na.rm=TRUE) +
      #   sum(dnorm(kappa[i,],0,kappa.prior.var^0.5,log=TRUE))
      # mh <- exp(mh1-mh2)
      if(all(1/phi.star.cams>0) & all(!is.infinite(1/phi.star.cams))){
        mh1 <- sum(dexp(t.staying.dat.all, 1/phi.star.cam.rep,log=TRUE),na.rm=TRUE) +
          sum(pexp(t.staying.dat.censor, 1/phi.star.cam.rep, lower.tail = F, log = TRUE),na.rm=TRUE) +
          sum(dnorm(kappa.star,0,kappa.prior.var^0.5,log=TRUE))
        mh2 <- sum(dexp(t.staying.dat.all, 1/phi.all.cam.rep,log=TRUE),na.rm=TRUE) +
          sum(pexp(t.staying.dat.censor, 1/phi.all.cam.rep, lower.tail = F, log = TRUE),na.rm=TRUE) +
          sum(dnorm(kappa[i,],0,kappa.prior.var^0.5,log=TRUE))
        mh <- exp(mh1-mh2)

        if(is.infinite(mh)){
          hh<-0
          # message("Warning: Inf vals")
        }
        if(is.na(mh)){
          message("ERROR: na vals")
        }

        if(mh>runif(1)){
          kappa[i,] <- kappa.star;
          accept[i+1,sum(covariates.index==1)+kk] <- 1;
          phi <- phi.star} else{
            kappa[i,] <- kappa[i,];
            accept[i+1,sum(covariates.index==1)+kk] <- 0}
      } else{
        kappa[i,] <- kappa[i,];
        accept[i+1,sum(covariates.index==1)+kk] <- 0}
    }
    kappa[i+1,] <- kappa[i,];
    phi <- exp(Z%*%kappa[i+1,])

    tot.u[i+1] <- sum(u)*cell.A

    # Update tuning parms
    if(i%%tune.check == 0){
      batch_n <- batch_n+1
      delta_n <- batch_n^-1
      # delta_n <- min(0.01,batch_n^-1)
      accept.gamma.check <- colMeans(accept[(i-tune.check+1):i,1:sum(covariates.index==1)],na.rm=T)
      gamma.tune[which(accept.gamma.check>0.44)] <- gamma.tune[which(accept.gamma.check>0.44)] + delta_n
      gamma.tune[which(accept.gamma.check<=0.44)] <- gamma.tune[which(accept.gamma.check<=0.44)] - delta_n
      accept.kappa.check <- colMeans(accept[(i-tune.check+1):i,(sum(covariates.index==1)+1):(2*sum(covariates.index==1))],na.rm=T)
      kappa.tune[which(accept.kappa.check>0.44)] <- kappa.tune[which(accept.kappa.check>0.44)] + delta_n
      kappa.tune[which(accept.kappa.check<=0.44)] <- kappa.tune[which(accept.kappa.check<=0.44)] - delta_n
    }
  }
  # print("MCMC complete")

  list(accept = accept,gamma = gamma,kappa = kappa,tot.u = tot.u,u = u)
}



########################################
# MCMC for MCT no covariates
########################################
fit.model.mcmc.MCT <- function(n.iter,
                               gamma.start,
                               gamma.prior.var,
                               gamma.tune,
                               cam.counts){

  # Variables that will be saved
  gamma <- matrix(,n.iter+1,1)
  tot.u <- matrix(,n.iter+1,1)
  accept <- matrix(,n.iter+1,1)

  gamma[1] <- gamma.start
  colnames(gamma) <- "gamma"
  colnames(tot.u) <- "Total estimate"
  colnames(accept) <- "accept.rate.gamma"

  d <- exp(gamma[1])
  tot.u <- d*tot.A/(cam.A*t.steps)

  tune.check <- 100
  batch_n <- 0

  # Begin MCMC loop
  # prog_bar <- txtProgressBar(min = 0,max = n.iter,style = 3,width = 50,char = "=")
  for(i in 1:n.iter){
    # setTxtProgressBar(prog_bar, i)
    #Sample gamma
    gamma.star <- rnorm(1,gamma[i],exp(2*gamma.tune))
    d.star <- exp(gamma.star)

    if(all(d.star>0)){
      mh1 <- sum(dpois(cam.counts,d.star,log=TRUE),na.rm=TRUE) +
        sum(dnorm(gamma.star,0,gamma.prior.var^0.5,log=TRUE))
      mh2 <- sum(dpois(cam.counts,d,log=TRUE),na.rm=TRUE) +
        sum(dnorm(gamma[i],0,gamma.prior.var^0.5,log=TRUE))
      mh <- exp(mh1-mh2)

      if(mh>runif(1)){gamma[i+1] <- gamma.star;
      accept[i+1] <- 1;
      d <- d.star} else{
        gamma[i+1] <- gamma[i];
        accept[i+1] <- 0}
    } else{
      gamma[i+1] <- gamma[i];
      accept[i+1] <- 0}

    d <- exp(gamma[i+1])
    tot.u[i+1] <- d*tot.A/(cam.A*t.steps)

    # Update tuning parms
    if(i%%tune.check == 0){
      batch_n <- batch_n+1
      delta_n <- batch_n^-1
      # delta_n <- min(0.01,batch_n^-1)
      accept.gamma.check <- mean(accept[(i-tune.check+1):i],na.rm=T)
      gamma.tune[accept.gamma.check>0.44] <- gamma.tune[accept.gamma.check>0.44] + delta_n
      gamma.tune[accept.gamma.check<=0.44] <- gamma.tune[accept.gamma.check<=0.44] - delta_n
    }

  }
  # print("MCMC complete")

  list(accept = accept,gamma = gamma,tot.u = tot.u)
}


########################################
# MCMC for MCT w/ covariates
########################################
fit.model.mcmc.MCT.cov <- function(n.iter,
                                   gamma.start,
                                   gamma.prior.var,
                                   gamma.tune,
                                   cam.counts,
                                   t.staying.dat,
                                   covariate_labels,
                                   covariates.index,
                                   cam.A,
                                   cell.A){

  # Variables that will be saved
  gamma <- matrix(,n.iter+1,sum(covariates.index==1))
  tot.u <- matrix(,n.iter+1,1)
  accept <- matrix(,n.iter+1,sum(covariates.index==1))
  gamma[1,] <- gamma.start
  colnames(gamma) <- paste0("gamma.",covariate_labels)
  colnames(tot.u) <- "Total estimate"
  colnames(accept) <- c(paste0("accept.rate.gamma.",
                              covariate_labels))
  # Initialize with landscape-scale covariates
  d <- exp(Z%*%gamma[1,])   # expected densities
  u <- d/t.steps
  tot.u[1] <- sum(u)

  tune.check <- 100
  batch_n <- 0

  # Begin MCMC loop
  # prog_bar <- txtProgressBar(min = 0,max = n.iter,style = 3,width = 50,char = "=")
  for(i in 1:n.iter){
    # setTxtProgressBar(prog_bar, i)
    #Sample gamma
    for(gg in 1:sum(covariates.index==1)){
      gamma.star <- gamma[i,]
      gamma.star[gg] <- rnorm(1,gamma[i,gg],exp(2*gamma.tune[gg]))
      d.star <- exp(Z%*%gamma.star)
      u.star <- d.star/t.steps

      # repeat estimated parms for fitting
      d.star.cams <- d.star[cam.samps] * cam.A / cell.A
      d.all.cams <- d[cam.samps] * cam.A / cell.A

      if(all(d.star.cams>0)){
        mh1 <- sum(dpois(cam.counts,d.star.cams,log=TRUE),na.rm=TRUE) +
          sum(dnorm(gamma.star,0,gamma.prior.var^0.5,log=TRUE))
        mh2 <- sum(dpois(cam.counts,d.all.cams,log=TRUE),na.rm=TRUE) +
          sum(dnorm(gamma[i,],0,gamma.prior.var^0.5,log=TRUE))
        mh <- exp(mh1-mh2)

        if(mh>runif(1)){gamma[i,] <- gamma.star;
        accept[i+1,gg] <- 1;
        u <- u.star;
        d <- d.star
        } else{
          gamma[i,] <- gamma[i,];
          accept[i+1,gg] <- 0
        }
      } else{
        gamma[i,] <- gamma[i,];
        accept[i+1,gg] <- 0}
    }
    gamma[i+1,] <- gamma[i,];
    d <- exp(Z%*%gamma[i+1,])
    u <- d/t.steps
    tot.u[i+1] <- sum(u)

    # Update tuning parms
    if(i%%tune.check == 0){
      batch_n <- batch_n+1
      delta_n <- batch_n^-1
      # delta_n <- min(0.01,batch_n^-1)
      accept.gamma.check <- colMeans(accept[(i-tune.check+1):i,1:sum(covariates.index==1)],na.rm=T)
      gamma.tune[which(accept.gamma.check>0.44)] <- gamma.tune[which(accept.gamma.check>0.44)] + delta_n
      gamma.tune[which(accept.gamma.check<=0.44)] <- gamma.tune[which(accept.gamma.check<=0.44)] - delta_n
    }
  }
  # print("MCMC complete")

  list(accept = accept,gamma = gamma,tot.u = tot.u, u = u)
}

########################################
# MCMC for STE no covariates
########################################
fit.model.mcmc.STE <- function(n.iter,gamma.start,gamma.prior.var,gamma.tune,
                               STE.dat,cam.A,cam.samps.in,censor){

  # Variables that will be saved
  gamma <- matrix(,n.iter+1,1)
  tot.u <- matrix(,n.iter+1,1)
  accept <- matrix(,n.iter+1,1)

  gamma[1] <- gamma.start
  colnames(gamma) <- "gamma"
  colnames(tot.u) <- "Total estimate"
  colnames(accept) <- "accept.rate.gamma"

  # Define censored space to event
  STE.dat[STE.dat == censor] <- NA
  STE.censor <- STE.dat
  STE.censor[] <- 0
  STE.censor[which(is.na(STE.dat))] <- censor

  u <- exp(gamma[1])
  tot.u <- u*tot.A

  tune.check <- 100
  batch_n <- 0

  # Begin MCMC loop
  # prog_bar <- txtProgressBar(min = 0,max = n.iter,style = 3,width = 50,char = "=")
  for(i in 1:n.iter){
    # setTxtProgressBar(prog_bar, i)
    #Sample gamma
    gamma.star <- rnorm(1,gamma[i],exp(2*gamma.tune))
    u.star <- exp(gamma.star)

    if(all(u.star>0)){
      mh1 <- sum(dexp(STE.dat,u.star,log=TRUE),na.rm=TRUE) +
        sum(pexp(STE.censor, u.star, lower.tail = F, log = TRUE)) +
        sum(dnorm(gamma.star,0,gamma.prior.var^0.5,log=TRUE))
      mh2 <- sum(dexp(STE.dat,u,log=TRUE),na.rm=TRUE) +
        sum(pexp(STE.censor, u, lower.tail = F, log = TRUE)) +
        sum(dnorm(gamma[i],0,gamma.prior.var^0.5,log=TRUE))
      mh <- exp(mh1-mh2)

      if(mh>runif(1)){gamma[i+1] <- gamma.star;
      accept[i+1] <- 1;
      u <- u.star} else{
        gamma[i+1] <- gamma[i];
        accept[i+1] <- 0}
    } else{
      gamma[i+1] <- gamma[i];
      accept[i+1] <- 0}

    u <- exp(gamma[i+1])
    tot.u[i+1] <- u*tot.A

    # Update tuning parms
    if(i%%tune.check == 0){
      batch_n <- batch_n+1
      delta_n <- batch_n^-1
      # delta_n <- min(0.01,batch_n^-1)
      accept.gamma.check <- mean(accept[(i-tune.check+1):i],na.rm=T)
      gamma.tune[accept.gamma.check>0.44] <- gamma.tune[accept.gamma.check>0.44] + delta_n
      gamma.tune[accept.gamma.check<=0.44] <- gamma.tune[accept.gamma.check<=0.44] - delta_n
    }

  }
  # print("MCMC complete")

  list(accept = accept,gamma = gamma,tot.u = tot.u)
}


########################################
# MCMC for SECR no covariates
########################################
# fit.model.mcmc.SECR <- function(n.iter,
#                                 psi.start,psi.prior.var,psi.tune,
#                                 sigma.start,sigma.prior.var,sigma.tune,
#                                 alpha.start,alpha.prior.var,alpha.tune,
#                                 z.dat,cam.A,cam.samps.in,censor){
#
#   # Variables that will be saved
#   psi <- matrix(,n.iter+1,M)
#   sigma <- matrix(,n.iter+1,1)
#   alpha <- matrix(,n.iter+1,1)
#   tot.u <- matrix(,n.iter+1,1)
#   accept <- matrix(,n.iter+1,1)
#
#   psi[1,] <- psi.start
#   sigma[1,] <- sigma.start
#   alpha[1] <- alpha.start
#   colnames(psi) <- "psi"
#   colnames(sigma) <- "sigma"
#   colnames(alpha) <- "alpha"
#   colnames(tot.u) <- "Total estimate"
#   colnames(accept) <- c("accept.rate.psi","accept.rate.sigma", "accept.rate.alpha")
#
#   # z[1,] <- rbern(M,psi[1])
#
#   lambda[1,] <- lambda_0 * exp(-d^2/(2*sigma[1,]^2)) * z[i,]
#
#   # u <- exp(psi[1])
#   # tot.u <- u*tot.A
#
#   tune.check <- 100
#   batch_n <- 0
#
#   # Begin MCMC loop
#   # prog_bar <- txtProgressBar(min = 0,max = n.iter,style = 3,width = 50,char = "=")
#   for(i in 1:n.iter){
#     # setTxtProgressBar(prog_bar, i)
#
#     # Maybe do a for loop through these
#     #Sample psi
#     psi.star <- rnorm(M,psi[i,],exp(2*psi.tune)) # psi.tune = 0.1
#     # z.star <- rbern(M,psi.star) # Not sure if this is correct
#
#     # if(all(u.star>0)){
#       mh1 <- sum(dbern(z.dat,psi.star,log=TRUE),na.rm=TRUE) +
#         sum(dnorm(psi.star,0,psi.prior.var^0.5,log=TRUE))
#       mh2 <- sum(dbern(z.dat,psi[i,],log=TRUE),na.rm=TRUE) +
#         sum(dnorm(psi[i,],0,psi.prior.var^0.5,log=TRUE))
#       mh <- exp(mh1-mh2)
#
#       if(mh>runif(1)){
#         psi[i+1,] <- psi.star;
#         accept[i+1] <- 1
#         } else{
#         psi[i+1,] <- psi[i,];
#         accept[i+1] <- 0
#         }
#     # } else{
#     #   psi[i+1] <- psi[i];
#     #   accept[i+1] <- 0
#     #   }
#
#     # Sample sigma
#     sigma.star <- rnorm(M, sigma[i,], exp(2*sigma.tune)) # sigma.tune = 0.1
#     p[i, j] <- psi[i]*plogis(sigma.star)*exp(-alpha1*d[i,j]*d[i,j])
#
#     lambda.star <- lambda_0 * exp(-d^2/(2*sigma.star^2)) * z[i,]
#
#     if (sigma.star > 0) {
#     mh1 <- sum(dpois(y.dat,lambda.star,log=TRUE),na.rm=TRUE) +
#       sum(dnorm(sigma.star,0,psi.prior.var^0.5,log=TRUE))
#     mh2 <- sum(dexp(y.dat,lambda[i,],log=TRUE),na.rm=TRUE) +
#       sum(dnorm(sigma[i,],0,psi.prior.var^0.5,log=TRUE))
#     mh <- exp(mh1-mh2)
#
#     if(mh>runif(1)){
#       sigma[i+1,] <- sigma.star;
#       accept[i+1] <- 1
#     } else{
#       sigma[i+1,] <- sigma[i,];
#       accept[i+1] <- 0
#     }
#     } else{
#       sigma[i+1] <- sigma[i,];
#       accept[i+1] <- 0
#     }
#
#     # Update tuning parms
#     if(i%%tune.check == 0){
#       batch_n <- batch_n+1
#       delta_n <- batch_n^-1
#       # delta_n <- min(0.01,batch_n^-1)
#       accept.psi.check <- mean(accept[(i-tune.check+1):i],na.rm=T)
#       psi.tune[accept.psi.check>0.44] <- psi.tune[accept.psi.check>0.44] + delta_n
#       psi.tune[accept.psi.check<=0.44] <- psi.tune[accept.psi.check<=0.44] - delta_n
#     }
#
#   }
#   # print("MCMC complete")
#
#   list(accept = accept,psi = psi,tot.u = tot.u)
# }

fit.model.mcmc.SECR <- function(SECR.dat,
                                psi.start,
                                sigma.start,
                                lambda_0.start,
                                S.tune,
                                sigma.tune,
                                lambda_0.tune,
                                X,
                                M,
                                x,
                                y,
                                n.iter) {

  # Variables that will be saved
  psi <- matrix(,n.iter+1,1)
  sigma <- matrix(,n.iter+1,1)
  lambda_0 <- matrix(,n.iter+1,1)
  tot.u <- matrix(,n.iter+1,1)
  accept <- matrix(,n.iter+1,2)

  psi[1] <- psi.start
  sigma[1] <- sigma.start
  lambda_0[1] <- lambda_0.start
  accept[1,] <- 0
  colnames(psi) <- "psi"
  colnames(sigma) <- "sigma"
  colnames(lambda_0) <- "lambda_0"
  colnames(tot.u) <- "Total estimate"
  colnames(accept) <- c("accept.rate.sigma", "accept.rate.lambda_0")

  tune.check <- 50
  batch_n <- 10

  #create initial values
  # home ranges
  S <- cbind(runif(M, x[1], x[2]), runif(M,y[1],y[2]))
  z <- rbinom(M,1,psi[1])
  seen <- apply(SECR.dat>0, 1, any)
  z[seen] <- 1 # set seen individuals' z=1
  tot.u[1] <- sum(z)

  #initiate distance matrix and lambda_ij matrix
  d <- e2dist(S, X)
  lambda<-lambda_0[1]*exp(-(d^2)/(2*sigma[1]^2))

  #start iterations of the chain
  for (iter in 1:n.iter) {

    #update sigma
    sigma.star <- rnorm(1, sigma[iter], exp(2*sigma.tune))
    if(sigma.star>0 & sigma.star < 10){   # automatically reject sigma.star that are < 10
    # if(sigma.star>0){   # automatically reject sigma.star that are <0
      lambda.star <- lambda_0[iter]*exp(-(d^2)/(2*sigma.star^2))
      mh1 <- sum(dpois(SECR.dat, lambda.star * z, log=TRUE))
      mh2 <- sum(dpois(SECR.dat, lambda * z, log=TRUE))
      mh <- exp(mh1-mh2)
      if(runif(1) < mh & !is.na(mh)){
        lambda <- lambda.star;
        sigma[iter+1] <- sigma.star;
        accept[iter+1,1] <- 1
      } else {
        sigma[iter+1] <- sigma[iter];
        accept[iter+1,1] <- 0
      }
    }else {
      sigma[iter+1] <- sigma[iter];
      accept[iter+1,1] <- 0
    }

    #update lambda_0
    lambda_0.star <- rnorm(1, lambda_0[iter], exp(2*lambda_0.tune))
    if(lambda_0.star>0){   #automatically reject lambda_0.star that are < 10
      lambda.star <- lambda_0.star*exp(-(d^2)/(2*sigma[iter+1]^2))
      mh1 <- sum(dpois(SECR.dat, lambda.star * z, log=TRUE))
      mh2 <- sum(dpois(SECR.dat, lambda * z, log=TRUE))
      mh <- exp(mh1-mh2)
      if(runif(1) < mh  & !is.na(mh)){
        lambda <- lambda.star
        lambda_0[iter+1] <- lambda_0.star
        accept[iter+1,2] <- 1
      } else {
        lambda_0[iter+1] <- lambda_0[iter]
        accept[iter+1,2] <- 0
      }
    } else {
      lambda_0[iter+1] <- lambda_0[iter]
      accept[iter+1,2] <- 0
    }

    #update z
    Zups <- 0
    for(i in 1:M) {
      if(!seen[i]) {#no need to update seen individuals, since their z =1
        z.star <- ifelse(z[i]==0, 1, 0)
        mh1 <- sum(dpois(SECR.dat[i,], lambda[i,] * z.star, log=TRUE)) +
          dbinom(z.star, 1, psi[iter], log=TRUE)
        mh2 <- sum(dpois(SECR.dat[i,], lambda[i,] * z[i], log=TRUE)) +
          dbinom(z[i], 1, psi[iter], log=TRUE)
        mh <- exp(mh1-mh2)

        if(runif(1) < mh & !is.na(mh)) {
          z[i] <- z.star
          Zups <- Zups + 1
        }
      }
    }#end M loop
    # print(mean(Zups))

    #update psi
    psi[iter+1] <- rbeta(1, 1+sum(z), 1 + M-sum(z))

    # update S
    Ssups <- 0
    for(i in 1:M) {
      S.star <- c(rnorm(1, S[i,1], S.tune), rnorm(1, S[i,2], S.tune))

      # Check to make sure S.star is within bounds
      inbox <-
        S.star[1] < x[2] &
        S.star[1] > x[1] &
        S.star[2] < y[2] &
        S.star[2] > y[1]

      if(inbox){
        dtmp <- sqrt((S.star[1] - X[,1])^2 + (S.star[2] - X[,2])^2)
        lambda.star<- lambda_0[iter+1]*exp(-(dtmp^2)/(2*sigma[iter+1]^2) )

        mh1 <- sum(dpois(SECR.dat[i,], lambda.star*z[i], log=TRUE))
        mh2 <- sum(dpois(SECR.dat[i,], lambda[i,]*z[i], log=TRUE))
        mh <- exp(mh1-mh2)
        if(runif(1) < mh & !is.na(mh)) {
          S[i,] <- S.star
          lambda[i,] <- lambda.star
          d[i,] <- dtmp
          Ssups <- Ssups + 1
        }
      }
    }#end M loop
    # print(Ssups)

    tot.u[iter+1] <- sum(z)

    # Update tuning parms
    if(iter%%tune.check == 0){
      batch_n <- batch_n+1
      delta_n <- batch_n^-1
      # delta_n <- min(0.01,batch_n^-1)
      accept.sigma.check <- mean(accept[(iter-tune.check+1):iter,1],na.rm=T)
      sigma.tune[accept.sigma.check>0.44] <- sigma.tune[accept.sigma.check>0.44] + delta_n
      sigma.tune[accept.sigma.check<=0.44] <- sigma.tune[accept.sigma.check<=0.44] - delta_n
      accept.lambda_0.check <- mean(accept[(iter-tune.check+1):iter,2],na.rm=T)
      lambda_0.tune[accept.lambda_0.check>0.44] <- lambda_0.tune[accept.lambda_0.check>0.44] + delta_n
      lambda_0.tune[accept.lambda_0.check<=0.44] <- lambda_0.tune[accept.lambda_0.check<=0.44] - delta_n
    }

  }
  # print("MCMC complete")

  list(accept = accept,sigma = sigma,psi = psi, lambda_0 = lambda_0, tot.u = tot.u)

}  #end of function call

SCR0pois <-
  function(y,X,M, xl,xu,yl,yu,delta, niter) {

    #create initial values
    S<-cbind(runif(M, xl, xu), runif(M,yl,yu))
    sigma<-runif(1,0.5,5)
    lam0<-runif(1,0.1,1)
    psi<-runif(1,0.2,0.8)
    z<-rbinom(M,1,psi)
    seen <- apply(y>0, 1, any)
    z[seen]<-1#set seen individuals' z=1

    #initiate distance matrix and lamij matrix
    d <- e2dist(S, X)
    lam<-lam0*exp(-(d*d)/(2*sigma*sigma))

    #set up matrix to hold results
    out<-matrix(nrow=niter, ncol=4)
    colnames(out)<-c('sigma', 'lam0', 'psi', 'N')

    #have R print starting values
    cat("\nstarting values =", c(sigma, lam0, psi, sum(z),"\n\n"))

    #start iterations of the chain

    for (iter in 1:niter) {

      #have R output the time and parameter estimates at every 100th iteration
      if(iter %% 100 == 0) {
        cat("iter", iter, format(Sys.time(), "%H:%M:%S"), "\n")
        cat("   current =", out[iter-1,], "\n")
      }

      #update sigma
      sig.cand <- rnorm(1, sigma, delta[1])
      if(sig.cand>0){   #automatically reject sig.cand that are <0
        lam.cand <- lam0*exp(-(d*d)/(2*sig.cand*sig.cand))
        ll<- sum(dpois(y, lam*z, log=TRUE))
        llcand <- sum(dpois(y, lam.cand*z, log=TRUE))
        if(runif(1) < exp( llcand  - ll) ){
          ll<-llcand
          lam<-lam.cand
          sigma<-sig.cand
        }
      }

      #update lam0
      lam0.cand <- rnorm(1, lam0, delta[2])
      if(lam0.cand>0){   #automatically reject lam0.cand that are <0
        lam.cand <- lam0.cand*exp(-(d*d)/(2*sigma*sigma))
        ll<- sum(dpois(y, lam*z, log=TRUE))
        llcand <- sum(dpois(y, lam.cand*z, log=TRUE))
        if(runif(1) < exp( llcand  - ll) ){
          ll<-llcand
          lam<-lam.cand
          lam0<-lam0.cand
        }
      }

      #update z
      zUps <- 0#set counter to monitor acceptance rate
      for(i in 1:M) {
        if(seen[i])#no need to update seen individuals, since their z =1
          next
        zcand <- ifelse(z[i]==0, 1, 0)
        llz <- sum(dpois(y[i,],lam[i,]*z[i], log=TRUE))
        llcand <- sum(dpois(y[i,], lam[i,]*zcand, log=TRUE))

        prior <- dbinom(z[i], 1, psi, log=TRUE)
        prior.cand <- dbinom(zcand, 1, psi, log=TRUE)
        if(runif(1) < exp( (llcand+prior.cand) - (llz+prior) )) {
          z[i] <- zcand
          zUps <- zUps+1
        }
      }#end M loop

      #update psi
      psi<-rbeta(1, 1+sum(z), 1 + M-sum(z))

      #update s
      Sups <- 0
      for(i in 1:M) {
        Scand <- c(rnorm(1, S[i,1], delta[3]), rnorm(1, S[i,2], delta[3]))
        inbox<-Scand[1]<xu & Scand[1]>xl & Scand[2]<yu & Scand[2]>yl
        if(inbox){
          dtmp <- sqrt((Scand[1] - X[,1])^2 + (Scand[2] - X[,2])^2)
          lam.cand<- lam0*exp(-(dtmp*dtmp)/(2*sigma*sigma) )

          llS <- sum(dpois(y[i,], lam[i,]*z[i], log=TRUE))
          llcand <- sum(dpois(y[i,], lam.cand*z[i], log=TRUE))
          if(runif(1) < exp(llcand - llS)) {
            S[i,] <- Scand
            lam[i,] <- lam.cand
            d[i,] <- dtmp
            Sups <- Sups+1
          }
        }
      }#end M loop

      #prompt R to output acceptance rates of z and S
      if(iter %% 100 == 0) {
        cat("   Acceptance rates\n")
        cat("     z =", zUps/M, "\n")
        cat("     S =", Sups/M, "\n")
      }

      out[iter,] <- c(sigma,lam0,psi,sum(z))

    }  #end of iteration loop

    return(out)

  }  #end of function call


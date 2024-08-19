#########################################
# Random walk simulations
#########################################
# Correlated random walk with grouping behavior
ABM_sim <- function(bounds,
                    t.steps,
                    speed_vals,
                    kappa,
                    clump_sizes,
                    clump.rad,
                    init_placement = 1) {
  
  # Initialize movement matrices
  animalxy.all <- data.frame(matrix(NA, 0, 4))
  colnames(animalxy.all) <- c("ID", "X", "Y", "t")
  
  # Grid lengths
  dx <- (bounds[2]-bounds[1])/q^0.5
  dy <- (bounds[2]-bounds[1])/q^0.5
  
  
  clump_IDs <- c(0, cumsum(clump_sizes))
  
  for (nc in 1:length(clump_sizes)) {
    
    clump_size <- clump_sizes[nc]
    # Initialize individual locations for all time steps
    X <- matrix(NA, clump_size, t.steps)
    Y <- matrix(NA, clump_size, t.steps)
    
    # Initialize herd centers
    if (init_placement == 1) {
      # Randomly place across landscape
      clump.X0 <- runif(1, min(bounds) + dx, max(bounds) - dx)
      clump.Y0 <- runif(1, min(bounds) + dy, max(bounds) - dy)
    } else {
      # Place according to inverse movement speeds
      t.stay.cdf <- cumsum(speed_vals^-1/sum(speed_vals^-1))
      t.stay.ind <- max(which(runif(1,0,1) > t.stay.cdf))
      clump.X0 <- ((t.stay.ind-1) %% q^0.5 + 1)/q^0.5
      clump.Y0 <- (ceiling(t.stay.ind/q^0.5))/q^0.5
    }
    
    # Initialize individuals near herd center
    X[,1] <- rtruncnorm(clump_size, min(bounds), max(bounds), clump.X0, clump.rad)
    Y[,1] <- rtruncnorm(clump_size, min(bounds), max(bounds), clump.Y0, clump.rad)
    
    # Random starting angle
    theta <- runif(1, 0, 2*pi)
    
    # corner locations of each grid cell
    grid.ints <- seq(min(bounds), max(bounds),by=max(bounds)/q^0.5)
    
    # Initialize location matrix for continuous time
    XY.all <- as.data.frame(matrix(NA, nrow = 0, ncol = 4))
    colnames(XY.all) <- c("ID","X","Y","t")
    XY.all[1:clump_size,] <- t(rbind((clump_IDs[nc]+1):clump_IDs[nc+1], X[,1], Y[,1], dt))
    
    # Loop through all time steps
    for(i in 2:(t.steps)) {
      # Turning angles for herd and all individuals
      theta <- rvonmises(1, theta, kappa)
      # A fixed concentration of 20 works for individuals
      theta.all <- rvonmises(clump_size, theta, 20)
      
      # Loop through all individuals in herd
      for(ci in 1:clump_size){
        # Time step
        t.step <- dt
        
        # Individual moves across landscape until time is depleted
        while(t.step > 0){
          # x,y indices of individual
          X.ind <- ceiling(X[ci,i-1]*(q^0.5/max(bounds)))
          Y.ind <- ceiling(Y[ci,i-1]*(q^0.5/max(bounds)))
          
          # Determine movement speed in current space
          v <- speed_vals[X.ind,Y.ind]
          step.size <- v*t.step
          
          # Step length components
          dX <- step.size*cos(as.numeric(theta.all[ci]))
          dY <- step.size*sin(as.numeric(theta.all[ci]))
          
          # Determine coordinate `corners' of grid for line segment calculations
          r.corners <- rbind(cbind(grid.ints[X.ind],grid.ints[Y.ind]),
                             cbind(grid.ints[X.ind],grid.ints[Y.ind]),
                             cbind(grid.ints[X.ind],grid.ints[Y.ind+1]),
                             cbind(grid.ints[X.ind+1],grid.ints[Y.ind]))
          
          # step lengths to add to corners
          r.length <- rbind(cbind(0,dy),cbind(dx,0),cbind(dx,0),cbind(0,dy))
          
          # check for intersections with grid
          # Don't use intersections when start point is on the line (a,b = 0),
          # or if a full step lands on a line (a,b = 1)
          a <- c()
          b <- c()
          int.check <- c()
          for(rr in 1:4){
            a[rr] <- det(rbind((c(X[ci,i-1],Y[ci,i-1])-r.corners[rr,]),
                               r.length[rr,]))/det(rbind(r.length[rr,],
                                                         c(dX,dY)))
            b[rr] <- det(rbind((c(X[ci,i-1],Y[ci,i-1])-r.corners[rr,]),
                               c(dX,dY)))/det(rbind(r.length[rr,],
                                                    c(dX,dY)))
          }
          
          int.check <- a<1 & a>0 & b<1 & b>0 |
            (a==0 & b<1 & b>0) | (a==1 & b<1 & b>0) | (b==0 & a<1 & a>0) | (b==1 & a<1 & a>0)
          
          # If individual crosses into new landscape cell, change movement speed
          if(any(int.check==T)){
            # only use first intersection crossing
            if(sum(int.check==T)>1){
              # coordinates for each intersection
              int.coords.all <- t(matrix(c(X[ci,i-1],Y[ci,i-1]),2,sum(int.check==T))) +
                a[int.check == TRUE]*t(matrix(c(dX,dY),2,sum(int.check==T)))
              
              vec.dist.all <- sqrt((t(matrix(c(X[ci,i-1]),1,sum(int.check==T)))-int.coords.all[,1])^2+
                                     (t(matrix(c(Y[ci,i-1]),1,sum(int.check==T)))-int.coords.all[,2])^2)
              int.T <- which(int.check==T)
              int.ind <- int.T[which(vec.dist.all == min(vec.dist.all))]
            }else{
              int.ind <- which(int.check == T)
            }
            
            # X,Y coordinates at intersection
            temp.X <- r.corners[int.ind,1] + b[int.ind]*r.length[int.ind,1]
            temp.Y <- r.corners[int.ind,2] + b[int.ind]*r.length[int.ind,2]
            
            # Set individual near grid line, adjust for new speed
            if (temp.X < (max(bounds)) & temp.X > (min(bounds)) &
                temp.Y < (max(bounds)) & temp.Y > (min(bounds))){
              t.step <- t.step - sqrt((X[ci,i-1]-(temp.X))^2+
                                        (Y[ci,i-1]-(temp.Y))^2)/v
              
              # perturb individual so that they are not on intersection
              X[ci,i-1] <- temp.X + 0.00001*sign(cos(theta.all[ci]))*(temp.X==r.corners[int.ind,1])
              Y[ci,i-1] <- temp.Y + 0.00001*sign(sin(theta.all[ci]))*(temp.Y==r.corners[int.ind,2])
              
              # Track all intermediate steps
              XY.all[nrow(XY.all)+1,] <- c(clump_IDs[nc] + ci, X[ci,i-1], Y[ci,i-1], i*dt - t.step)
              
            }
            # If encounter border, change turning angle
            else{
              if (temp.X >= (max(bounds)) || temp.X <= (min(bounds))){
                theta.x <- cos(theta.all[ci])
                theta.y <- sin(theta.all[ci])
                
                # Reflective angle
                theta.all[ci] <- atan(-theta.y/theta.x) + pi * (theta.x>0)
                
                # perturb individual so that they are not on intersection
                X[ci,i-1] <- temp.X + 0.00001*sign(cos(theta.all[ci]))*(temp.X==r.corners[int.ind,1])
                Y[ci,i-1] <- temp.Y + 0.00001*sign(sin(theta.all[ci]))*(temp.Y==r.corners[int.ind,2])
                
                # Designate one individual to determine herd turning angle
                if(ci == 1){
                  theta <- theta.all[ci]
                }
              }
              if(temp.Y >= (max(bounds)) || temp.Y <= (min(bounds))) {
                theta.x <- cos(theta.all[ci])
                theta.y <- sin(theta.all[ci])
                # Reflective angle
                theta.all[ci] <- atan(-theta.y/theta.x) + pi * (theta.x<0)
                
                # perturb individual so that they are not on intersection
                X[ci,i-1] <- temp.X + 0.00001*sign(cos(theta.all[ci]))*(temp.X==r.corners[int.ind,1])
                Y[ci,i-1] <- temp.Y + 0.00001*sign(sin(theta.all[ci]))*(temp.Y==r.corners[int.ind,2])
                
                # designate one individual to determine herd movement
                if(ci == 1){
                  theta <- theta.all[ci]
                }
              }
            }
          }
          # Individuals take step within cell
          else{
            X[ci,i] <- X[ci,i-1] + dX
            Y[ci,i] <- Y[ci,i-1] + dY
            t.step <- 0
            
            # Track all intermediate steps
            XY.all[nrow(XY.all)+1,] <- c(clump_IDs[nc] + ci, X[ci,i], Y[ci,i], i*dt)
          }
        }
      }
    }
    
    animalxy.all$lscape_type <- lscape_speeds$Speed[animalxy.all$XY_inds]
    
    animalxy.all <- rbind(animalxy.all, XY.all)
  }
  
  # Convert animalxy 2D coordinates to single digit coordinates
  animalxy.all <- mutate(animalxy.all,
                         XY_inds = ceiling(animalxy.all$X/dx) +
                           floor(animalxy.all$Y/dx)*q^0.5,
                         ii = 1:nrow(animalxy.all)) |>
    group_by(ID) |>
    mutate(dist = c(sqrt((X[2:length(ID)] - X[1:(length(ID)-1)])^2 + (Y[2:length(ID)] - Y[1:(length(ID)-1)])^2), NA),
           speeds = dist/c(t[2:length(ID)] - t[1:(length(ID)-1)], NA))
  
  # Escapee
  if(max(animalxy.all$X)>max(bounds) || min(animalxy.all$X)<min(bounds)){
    stop("Animal out of bounds")
  }
  
  return(animalxy.all)
}

create_lscape <- function(Design = "Random", speed_bounds = c(0.9,1.1), custom_inds) {
  
  # # Initialize speed dataframe
  speed_vals <- data.frame(Index = 1:q,
                              Speed_index = sample(nrow(speed_bounds),
                                                   q,
                                                   replace = T)) |>
    mutate(Speed = speed_bounds$Speed[Speed_index],
           Value = runif(q,
                         speed_bounds$Min[Speed_index],
                         speed_bounds$Max[Speed_index]))
  
  if(Design == "Random") {
    # Shuffle indices
    shuffle_inds <- sample(1:q, q, replace = F)
    
    if(nrow(speed_bounds) == 1){
      speed_vals[speed_inds] <- dx*runif(q, speed_bounds[1], speed_bounds[2])
    } else {
      subset_inds <- round(seq(1, q, length.out = nrow(speed_bounds) + 1))
      
      speed_inds <- c()
      # Loop through all speeds
      for (speeds in 1:nrow(speed_bounds)) {
        speed_inds[[speeds]] <- list(shuffle_inds[subset_inds[speeds]:subset_inds[speeds+1]])
        speed_vals[unlist(speed_inds[speeds])] <- dx*runif(length(unlist(speed_inds[speeds])),
                                                              speed_bounds[speeds, 1],
                                                              speed_bounds[speeds, 2])
      }
    }
  }  else if (Design == "Custom") {
    for (speeds in 1:nrow(speed_bounds)) {
      speed_inds <- unlist(custom_inds[speeds])
      speed_vals[speed_inds] <- dx*runif(length(speed_inds),
                                            speed_bounds[speeds, 1],
                                            speed_bounds[speeds, 2])
    }
  }
  # return(speed_vals)
  list(speed_vals = speed_vals, speed_inds = speed_inds)
}



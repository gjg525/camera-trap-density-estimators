#########################################
# Random walk simulations
#########################################
# Correlated random walk with grouping behavior
ABM_sim <- function(bounds,
                    t.steps,
                    speeds,
                    direction,
                    kappa,
                    road,
                    clump_sizes,
                    clump.rad,
                    init_placement = 1) {

  # Initialize movement matrices
  animalxy.all <- data.frame(Animal_ID = double(),
                             x = double(),
                             y = double(),
                             t = double())

  clump_IDs <- c(0, cumsum(clump_sizes))

  # corner locations of each grid cell
  grid.ints <- seq(min(bounds), max(bounds), by = max(bounds) / q ^ 0.5)

  # Loops through all clumps
  for (nc in 1:length(clump_sizes)) {
    # Get clulmp size for current clump
    clump_size <- clump_sizes[nc]

    # Initialize individual locations for all time steps
    X <- matrix(NA, clump_size, t.steps + 1)
    Y <- matrix(NA, clump_size, t.steps + 1)

    # Initialize herd centers
    if (init_placement == 1) {
      # Randomly place across landscape
      clump_IC <- data.frame(IC_index = sample(which(!is.na(road)), 1)) |>
        dplyr::summarise(x = ((IC_index-1) %% q^0.5 + 1)*dx,
                  y = (ceiling(IC_index/q^0.5))*dx)
    } else {
      # Place according to inverse movement speeds
      clump_IC <- data.frame(IC_speeds = cumsum(speeds^-1/sum(speeds^-1)),
                             IC_index = max(which(runif(1,0,1) > t.stay.cdf))) |>
        dplyr::summarise(x = ((IC_index-1) %% q^0.5 + 1)*dx,
                  y = (ceiling(IC_index/q^0.5))*dx)
    }

    # # Initialize individuals in grid cell
    X[,1] <- runif(clump_size)*dx + clump_IC$x - dx
    Y[,1] <- runif(clump_size)*dy + clump_IC$y - dy

    h_range_center <- c(X[, 1], Y[, 1])

    # Random starting angle
    theta <- runif(1, 0, 2*pi)

    # Temporary location matrix for continuous time
    clumpxy.all <- data.frame(Animal_ID = double(),
                               x = double(),
                               y = double(),
                               t = double())
    clumpxy.all[1:clump_size,] <- t(rbind((clump_IDs[nc]+1):clump_IDs[nc+1],
                                          X[,1],
                                          Y[,1],
                                          0))

    # Loop through all time steps
    for(i in 1:t.steps) {
      # Turning angles for herd and all individuals
      theta <- Rfast::rvonmises(1,
                         theta,
                         kappa[ceiling(X[1,i]*(q^0.5/max(bounds))),
                               ceiling(Y[1,i]*(q^0.5/max(bounds)))])
      # A fixed concentration of 20 works for individuals
      theta.all <- Rfast::rvonmises(clump_size, theta, 20)

      # Loop through all individuals in herd
      for(ci in 1:clump_size){
        # Time step
        t.step <- dt

        # Individual moves across landscape until time is depleted
        while(t.step > 0){
          # x,y indices of individual
          X.ind <- ceiling(X[ci,i]*(q^0.5/max(bounds)))
          Y.ind <- ceiling(Y[ci,i]*(q^0.5/max(bounds)))

          # Determine movement speed in current space
          v <- speeds[X.ind,Y.ind]
          step.size <- v*t.step

          # Individual follows herd unless in directional cell
          theta.bias <- direction[X.ind,Y.ind]
          theta.all[ci] <- dplyr::if_else(is.na(theta.bias),
                                  theta.all[ci],
                                  Rfast::rvonmises(1, theta.bias, kappa[X.ind,Y.ind]))

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
            a[rr] <- det(rbind((c(X[ci,i],Y[ci,i])-r.corners[rr,]), r.length[rr,]))/
              det(rbind(r.length[rr,], c(dX,dY)))
            b[rr] <- det(rbind((c(X[ci,i],Y[ci,i])-r.corners[rr,]), c(dX,dY)))/
              det(rbind(r.length[rr,], c(dX,dY)))
          }

          int.check <- a<1 & a>0 & b<1 & b>0 |
            (a==0 & b<1 & b>0) |
            (a==1 & b<1 & b>0) |
            (b==0 & a<1 & a>0) |
            (b==1 & a<1 & a>0)

          # If individual crosses into new landscape cell, change movement speed
          if(any(int.check==T)){
            # only use first intersection crossing
            if(sum(int.check==T)>1){
              # coordinates for each intersection
              int.coords.all <- t(matrix(c(X[ci,i],Y[ci,i]),2,sum(int.check==T))) +
                a[int.check == TRUE]*t(matrix(c(dX,dY),2,sum(int.check==T)))

              vec.dist.all <- sqrt((t(matrix(c(X[ci,i]),1,sum(int.check==T)))-int.coords.all[,1])^2+
                                     (t(matrix(c(Y[ci,i]),1,sum(int.check==T)))-int.coords.all[,2])^2)
              int.T <- which(int.check==T)
              int.ind <- int.T[which(vec.dist.all == min(vec.dist.all))]
            }else{
              int.ind <- which(int.check == T)
            }

            # X,Y coordinates at intersection
            temp.X <- r.corners[int.ind,1] + b[int.ind]*r.length[int.ind,1]
            temp.Y <- r.corners[int.ind,2] + b[int.ind]*r.length[int.ind,2]

            # Set individual near grid line, adjust for new speed
            # Determine whether the lscape type is "blocked" or not
            # Do this better
            tx <- ceiling((temp.X + 0.00001*sign(cos(theta.all[ci]))*
                          (temp.X==r.corners[int.ind,1]))*(q^0.5/max(bounds)))
            ty <- ceiling((temp.Y + 0.00001*sign(sin(theta.all[ci]))*
                          (temp.Y==r.corners[int.ind,2]))*(q^0.5/max(bounds)))
            tx <- min(c(tx, dim(road)[1]))
            tx <- max(c(tx, 1))
            ty <- min(c(ty, dim(road)[1]))
            ty <- max(c(ty, 1))
            if (temp.X < (max(bounds)) & temp.X > (min(bounds)) &
                temp.Y < (max(bounds)) & temp.Y > (min(bounds))
                ){
              if (!is.na(road[tx, ty])) {
                t.step <- t.step - sqrt((X[ci,i]-(temp.X))^2 +
                                        (Y[ci,i]-(temp.Y))^2) / v

                # perturb individual so that they are not on intersection
                X[ci,i] <- temp.X + 0.00001*sign(cos(theta.all[ci]))*(temp.X==r.corners[int.ind,1])
                Y[ci,i] <- temp.Y + 0.00001*sign(sin(theta.all[ci]))*(temp.Y==r.corners[int.ind,2])

                # Track all intermediate steps
                clumpxy.all[nrow(clumpxy.all)+1,] <- c(clump_IDs[nc] + ci,
                                                       X[ci,i],
                                                       Y[ci,i],
                                                       i*dt - t.step)
              }
              else {
                t.step <- t.step - sqrt((X[ci,i]-(temp.X))^2 +
                                          (Y[ci,i]-(temp.Y))^2) / v

                # Reflect off of NA cells
                theta.x <- cos(theta.all[ci])
                theta.y <- sin(theta.all[ci])

                theta.all[ci] <- dplyr::if_else(abs(tx - ceiling(X[ci,i]*(q^0.5/max(bounds)))) > 0,
                                        (atan(-theta.y/theta.x) + pi * (theta.x>0)) %% (2*pi),
                                        (atan(-theta.y/theta.x) + pi * (theta.x<0)) %% (2*pi))
                # perturb individual so that they are not on intersection
                X[ci,i] <- temp.X + 0.00001*sign(cos(theta.all[ci]))*(temp.X==r.corners[int.ind,1])
                Y[ci,i] <- temp.Y + 0.00001*sign(sin(theta.all[ci]))*(temp.Y==r.corners[int.ind,2])

                # Track all intermediate steps
                clumpxy.all[nrow(clumpxy.all)+1,] <- c(clump_IDs[nc] + ci,
                                                       X[ci,i],
                                                       Y[ci,i],
                                                       i*dt - t.step)
                # designate one individual to determine herd movement
                if(ci == 1){
                  theta <- theta.all[ci]
                }
              }
            }
            # If encounter border, change turning angle
            else{
              if (temp.X >= (max(bounds)) || temp.X <= (min(bounds))){
                t.step <- t.step - sqrt((X[ci,i]-(temp.X))^2 +
                                          (Y[ci,i]-(temp.Y))^2) / v

                theta.x <- cos(theta.all[ci])
                theta.y <- sin(theta.all[ci])

                # Reflective angle
                theta.all[ci] <- (atan(-theta.y/theta.x) + pi * (theta.x>0)) %% (2*pi)

                # perturb individual so that they are not on intersection
                X[ci,i] <- temp.X + 0.00001*sign(cos(theta.all[ci]))*(temp.X==r.corners[int.ind,1])
                Y[ci,i] <- temp.Y + 0.00001*sign(sin(theta.all[ci]))*(temp.Y==r.corners[int.ind,2])

                # Track all intermediate steps
                clumpxy.all[nrow(clumpxy.all)+1,] <- c(clump_IDs[nc] + ci,
                                                       X[ci,i],
                                                       Y[ci,i],
                                                       i*dt - t.step)
                # Designate one individual to determine herd turning angle
                if(ci == 1){
                  theta <- theta.all[ci]
                }
              }
              if(temp.Y >= (max(bounds)) || temp.Y <= (min(bounds))) {
                t.step <- t.step - sqrt((X[ci,i]-(temp.X))^2 +
                                          (Y[ci,i]-(temp.Y))^2) / v

                theta.x <- cos(theta.all[ci])
                theta.y <- sin(theta.all[ci])
                # Reflective angle
                theta.all[ci] <- (atan(-theta.y/theta.x) + pi * (theta.x<0)) %% (2*pi)

                # perturb individual so that they are not on intersection
                X[ci,i] <- temp.X + 0.00001*sign(cos(theta.all[ci]))*(temp.X==r.corners[int.ind,1])
                Y[ci,i] <- temp.Y + 0.00001*sign(sin(theta.all[ci]))*(temp.Y==r.corners[int.ind,2])

                # Track all intermediate steps
                clumpxy.all[nrow(clumpxy.all)+1,] <- c(clump_IDs[nc] + ci,
                                                       X[ci,i],
                                                       Y[ci,i],
                                                       i*dt - t.step)
                # designate one individual to determine herd movement
                if(ci == 1){
                  theta <- theta.all[ci]
                }
              }
            }
          }
          # Individuals take step within cell
          else{
            X[ci,i + 1] <- X[ci,i] + dX
            Y[ci,i + 1] <- Y[ci,i] + dY
            t.step <- 0

            # Track all intermediate steps
            clumpxy.all[nrow(clumpxy.all)+1,] <- c(clump_IDs[nc] + ci,
                                                   X[ci,i + 1],
                                                   Y[ci,i + 1],
                                                   i*dt)
          }
        }
      }
      # clump direction is determined from 1 representative animal
      theta <- theta.all[1]
    }

    animalxy.all <- rbind(animalxy.all, clumpxy.all)
  }

  # Convert animalxy 2D coordinates to single digit coordinates
  animalxy.all <- dplyr::mutate(animalxy.all,
                         lscape_index = ceiling(animalxy.all$x/dx) +
                           floor(animalxy.all$y/dx)*q^0.5,
                         ii = 1:nrow(animalxy.all)) |>
    dplyr::group_by(Animal_ID) |>
    dplyr::mutate(dist = c(sqrt((x[2:length(Animal_ID)] - x[1:(length(Animal_ID)-1)])^2 +
                         (y[2:length(Animal_ID)] - y[1:(length(Animal_ID)-1)])^2), NA),
           speeds = dist/c(t[2:length(Animal_ID)] - t[1:(length(Animal_ID)-1)], NA))

  # Escapee
  if(max(animalxy.all$x)>max(bounds) || min(animalxy.all$x)<min(bounds)){
    stop("Animal out of bounds")
  }

  # Add speed and road info
  animalxy.all$lscape_type <- lscape_speeds$Speed[animalxy.all$lscape_index]
  animalxy.all$road <- lscape_speeds$Road[animalxy.all$lscape_index]

  return(animalxy.all)
}

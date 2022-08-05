#########################################
# Random walks for nind animals
#########################################
# Correlated random walk with grouping behavior
clump.corr.walk <- function(bounds, 
                            t.steps, 
                            steplen, 
                            dx, 
                            dy, 
                            v.abm,
                            t.stay.cdf, 
                            corr.walk.kappa,
                            clump.size, 
                            clump.rad) {
  
  X <- matrix(NA, clump.size, t.steps)
  Y <- matrix(NA, clump.size, t.steps)
  
  # Randomly place herds across landscape
  clump.X0 <- runif(1,min(bounds)+dx,max(bounds)-dx)
  clump.Y0 <- runif(1,min(bounds)+dy,max(bounds)-dy)
  
  # # Initialize based on spatial covariates
  # t.stay.samp <- runif(1,0,1)
  # t.stay.ind <- max(which(t.stay.samp>t.stay.cdf))
  # clump.X0 <- ((t.stay.ind-1) %% q^0.5 +1)/q^0.5
  # clump.Y0 <- (ceiling(t.stay.ind/q^0.5))/q^0.5
  
  # Initialize clump such that individuals are near each other
  X[,1] <- rtruncnorm(clump.size,min(bounds),max(bounds),clump.X0,clump.rad)
  Y[,1] <- rtruncnorm(clump.size,min(bounds),max(bounds),clump.Y0,clump.rad)
  
  # Initial starting angle
  theta <- runif(1, 0, 2*pi)
  
  # corner locations of each grid cell
  grid.ints <- seq(min(bounds), max(bounds),by=max(bounds)/q^0.5)
  
  XY.all <- as.data.frame(matrix(NA, nrow = 0, ncol = 4))
  colnames(XY.all) <- c("ID","X","Y","t")
  XY.all[1:clump.size,] <- t(rbind(1:clump.size, X[,1], Y[,1], 1))
  
  foo <- 5
  for(i in 2:(t.steps)) {
    # Correlated random walk for group
    theta <- rvonmises(1, theta, corr.walk.kappa)
    theta.all <- rvonmises(clump.size, theta, 20)
    
    # Define movement for each individual in clump
    for(ci in 1:clump.size){
      # time step
      t.step <- dt
      
      while(t.step > 0){
        # x,y indices of individual
        X.ind <- ceiling(X[ci,i-1]*(q^0.5/max(bounds)))
        Y.ind <- ceiling(Y[ci,i-1]*(q^0.5/max(bounds)))
        
        # Determine motility in current space
        mu.s <- v.abm[X.ind,Y.ind]
        step.size <- mu.s*t.step
        
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
        # Account for points that start on grid
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
        
        # If individual's trajectory crosses grid lines, change movement speed
        if(any(int.check==T)){
          # only use shortest distance intersection
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
          
          # set individual near grid line, adjust for new speed
          if (temp.X < (max(bounds)) & temp.X > (min(bounds)) &
              temp.Y < (max(bounds)) & temp.Y > (min(bounds))){
            t.step <- t.step - sqrt((X[ci,i-1]-(temp.X))^2+
                                      (Y[ci,i-1]-(temp.Y))^2)/mu.s
            
            # perturb individual so that they are not on intersection
            X[ci,i-1] <- temp.X + 0.00001*sign(cos(theta.all[ci]))*(temp.X==r.corners[int.ind,1])
            Y[ci,i-1] <- temp.Y + 0.00001*sign(sin(theta.all[ci]))*(temp.Y==r.corners[int.ind,2])
            
            # Track all intermediate steps
            XY.all[nrow(XY.all)+1,] <- c(ci, X[ci,i-1], Y[ci,i-1], i - t.step)

          }
          # If encounter border, change turning angle
          else{
            if (temp.X >= (max(bounds)) || temp.X <= (min(bounds))){
              theta.x <- cos(theta.all[ci])
              theta.y <- sin(theta.all[ci])
              # Reflective angle
              theta.all[ci] <- atan(-theta.y/theta.x) + pi * (theta.x>0)
              
              # designate one individual to determine herd movement
              if(ci == 1){
                theta <- theta.all[ci]
              }
            }
            if(temp.Y >= (max(bounds)) || temp.Y <= (min(bounds))) {
              theta.x <- cos(theta.all[ci])
              theta.y <- sin(theta.all[ci])
              # Reflective angle
              theta.all[ci] <- atan(-theta.y/theta.x) + pi * (theta.x<0)
              
              # designate one individual to determine herd movement
              if(ci == 1){
                theta <- theta.all[ci]
              }
            }
          }
        }
        # Individuals take step as usual
        else{
          X[ci,i] <- X[ci,i-1] + dX
          Y[ci,i] <- Y[ci,i-1] + dY
          t.step <- 0
          
          # Track all intermediate steps
          XY.all[nrow(XY.all)+1,] <- c(ci, X[ci,i], Y[ci,i], i)
        }
      }
      
      
    }
  }
  # out <- data.frame(t(rbind(as.vector(t(X)),as.vector(t(Y)))))
  # return(out)
  
  list(animalxy = data.frame(t(rbind(as.vector(t(X)),as.vector(t(Y))))),
       XY.all = XY.all)
}


# Correlated random walk with clumping behavior
animalxy <- data.frame(matrix(NA, t.steps*sum(clump.size),4))
animalxy.all <- data.frame(matrix(NA, 0, 4))
colnames(animalxy) <- c("X", "Y", "t", "ID")
colnames(animalxy.all) <- c("ID", "X", "Y", "t")
t.inds <- 1
ind.index <- 1
# print("ABM progress")
# abm_pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
#                          max = num.clumps, # Maximum value of the progress bar
#                          style = 3,    # Progress bar style (also available style = 1 and style = 2)
#                          width = 50,   # Progress bar width. Defaults to getOption("width")
#                          char = "=")   # Character used to create the bar
for(nc in 1:num.clumps){
  # setTxtProgressBar(abm_pb, nc)
  abm.out <- clump.corr.walk(bounds, t.steps, steplen, dx, dy, v.abm,t.stay.cdf, corr.walk.kappa,clump.size[nc],clump.rad = clump.rad)
  XY <- abm.out$animalxy
  XY.all <- abm.out$XY.all
  XY.all <- XY.all[order(XY.all$ID),]
  XY.all$ID <- XY.all$ID + ind.index - 1
  tot.inds <- t.steps*clump.size[nc]
  animalxy[t.inds:(t.inds+(t.steps*clump.size[nc])-1),] <-
    cbind(XY,rep(1:t.steps, times = clump.size[nc]),rep(ind.index:(ind.index+clump.size[nc]-1), each = t.steps))
  animalxy.all <- rbind(animalxy.all, XY.all)
  t.inds <- t.inds+(t.steps*clump.size[nc])
  ind.index <- ind.index + clump.size[nc]
}
# print("ABM sims complete")

# Escapee
if(max(animalxy[,1:2])>max(bounds) || min(animalxy[,1:2])<min(bounds)){
  stop("Animal index out of bounds")
}

# Check for NA values
if(any(is.na(animalxy))){
  stop("NA vals in simulated data")
}

# Create raster object for simulation distributions
u.abm.all <- raster(, nrows = q^0.5, ncols = q^0.5, xmn = 0,
                    xmx = 1, ymn = 0, ymx = 1, crs = NA)
u.abm.all[] <- 0
u.abm.all <- stack(replicate(t.steps,u.abm.all))

for (tt in 1:t.steps) {
  u.temp = matrix(0, nrow = q^0.5, ncol = q^0.5)
  ind.seq <- seq(tt,t.steps*nind,t.steps)
  for(inds in 1:nind) {
    x.round <- ceiling(animalxy[ind.seq[inds],1]*(q^0.5/max(bounds)))
    y.round <- ceiling(animalxy[ind.seq[inds],2]*(q^0.5/max(bounds)))
    
    # Transpose for converting matrix to raster
    u.temp[x.round,y.round] <- u.temp[x.round,y.round]+1
  }
  u.abm.all[[tt]] <- raster(t(u.temp))
}

# Convert animalxy 2D coordinates to single digit coordinates
animalxy.inds.2D <- ceiling(animalxy[,1:2])
animalxy.inds <- animalxy.inds.2D[,1]+(animalxy.inds.2D[,2]-1)*q^0.5

animalxy.all <- mutate(animalxy.all,
                       XY_inds = ceiling(animalxy.all$X) +
                         floor(animalxy.all$Y)*q^0.5,
                       ii = 1:nrow(animalxy.all))

# # Plot ABM simulations
b.df <- data.frame(c(bounds,dx))
g <- ggplot( ) +
  geom_path(data = animalxy,
            aes(animalxy[,1], animalxy[,2], col = as.factor(animalxy[,4]))) +
  geom_rect(data = b.df,
            aes(xmin = bounds[1], xmax = bounds[2], ymin = bounds[1], ymax = bounds[2]),
            color="black",
            fill= NA) +
  theme(legend.position = "none")
g


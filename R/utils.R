calc_tri_area <- function(x,y) {
  # Calculate area of triangle using vertices
  tri_area <- 0.5 * abs(x[1]*(y[2] - y[3]) + x[2]*(y[3] - y[1]) + x[3] * (y[1] - y[2]))
}

# Biased sample design (lscape_speeds)
sample_speeds <- function(cam.dist.set) {
  ps <- ncam*c(0.1,0.1,0.1)
  ps[cam.dist.set-1] <- ncam*0.8
  cam.samps <- c(sample(lscape_speeds |>
                          filter(Speed == "Slow") |>
                          pull(Index),
                        ps[1],
                        replace=F),
                 sample(lscape_speeds |>
                          filter(Speed == "Medium") |>
                          pull(Index),
                        ps[2],
                        replace=F),
                 sample(lscape_speeds |>
                          filter(Speed == "Fast") |>
                          pull(Index),
                        ps[3],
                        replace=F))

  # cam.samps <- c(sample(which(lscape_speeds$Speed == "Slow"), ncam*slow, replace = F),
  #                sample(which(lscape_speeds$Speed == "Medium"), ncam*medium, replace = F),
  #                sample(which(lscape_speeds$Speed == "Fast"), ncam*fast, replace = F))

}

# Biased sample design (road)
sample_roads <- function(on = 0.4, off = 0.6) {
  cam.samps <- c(sample(which(lscape_speeds$Road == "On Trail"), ncam*on, replace = F),
                 sample(which(lscape_speeds$Road == "Off Trail"), ncam*off, replace = F))
}

###################################
# Data collection functions
###################################

# Determine if a point lies within a triangle
calc_in_triangle <- function(p_viewshed, point) {

  A_1 <- calc_tri_area(t(rbind(p_viewshed[1:2,1], point[1])),
                       t(rbind(p_viewshed[1:2,2], point[2])))
  A_2 <- calc_tri_area(t(rbind(p_viewshed[c(1,3),1], point[1])),
                       t(rbind(p_viewshed[c(1,3),2], point[2])))
  A_3 <- calc_tri_area(t(rbind(p_viewshed[2:3,1], point[1])),
                       t(rbind(p_viewshed[2:3,2], point[2])))

  tri_area <- calc_tri_area(t(p_viewshed[,1]), t(p_viewshed[, 2]))

  in_triangle <- ifelse((A_1 + A_2 + A_3 - 0.000001) > tri_area,
                        FALSE,
                        TRUE)
}

# Calculate where an animal's trajectory intersects with a camera viewshed and the time spent within viewshed
calc_intersects <- function(p_viewshed, p_animal, speed) {
  X <- c()
  Y <- c()
  t_tot <- 0 # time spent in camera

  # Get lengths of all edges
  r.length <- rbind(p_viewshed[2,] - p_viewshed[1,],
                    p_viewshed[3,] - p_viewshed[2,],
                    p_viewshed[1,] - p_viewshed[3,])

  for (nn in 1:(nrow(p_animal)-1)) {
    X.temp <- c()
    Y.temp <- c()

    animal.step <- p_animal[nn+1,] - p_animal[nn,]

    # check for intersections with viewshed
    # Don't use intersections when start point is on the line (a,b = 0),
    # or if a full step lands on a line (a,b = 1)
    a <- c()
    b <- c()
    int.check <- c()
    for(rr in 1:3){
      a[rr] <- det(data.matrix(rbind(p_animal[nn,] - p_viewshed[rr,], r.length[rr,])))/
        det(data.matrix(rbind(r.length[rr,], animal.step)))
      b[rr] <- det(data.matrix(rbind((p_animal[nn,] - p_viewshed[rr,]), animal.step)))/
        det(data.matrix(rbind(r.length[rr,], animal.step)))
    }

    int.check <- a<1 & a>0 & b<1 & b>0 |
      (a==0 & b<1 & b>0) |
      (a==1 & b<1 & b>0) |
      (b==0 & a<1 & a>0) |
      (b==1 & a<1 & a>0)

    # Check if animal is already inside viewshed
    in_cam <- calc_in_triangle(p_viewshed, p_animal[nn,])

    if(in_cam) {
      # Individual steps outside of camera
      if(any(int.check==T)){

        int.ind <- which(int.check == T)

        # X,Y coordinates at intersection
        X.temp <- p_viewshed[int.ind,1] + b[int.ind]*r.length[int.ind,1]
        Y.temp <- p_viewshed[int.ind,2] + b[int.ind]*r.length[int.ind,2]

        # Find distance from start point to intersection
        dist_1 <- sqrt((X.temp - p_animal[nn, 1])^2 +
                         (Y.temp - p_animal[nn, 2])^2)

        # Calculate time spent in camera
        t_tot <- t_tot + dist_1/speed[nn]

        # plot(rbind(p_viewshed, p_viewshed[1,]), type = "b", col = "black")
        # lines(p_animal, col = "blue")
        # lines(X.temp,Y.temp, col = "red", pch = 18)
      } else {
        # Individual remains in camera

        # Find distance from start point to end point
        dist_1 <- sqrt((p_animal[nn + 1, 1] - p_animal[nn, 1])^2 +
                         (p_animal[nn + 1, 2] - p_animal[nn, 2])^2)

        # Calculate time spent in camera
        t_tot <- t_tot + dist_1/speed[nn]

        # plot(rbind(p_viewshed, p_viewshed[1,]), type = "b", col = "black")
        # lines(p_animal, col = "blue")
        # lines(X.temp,Y.temp, col = "red", pch = 18)
      }
    } else{
      # Individual starts outside of camera
      if(sum(int.check) == 1) {
        # Individual steps into and stays in camera

        int.ind <- which(int.check == T)

        # X,Y coordinates at intersection
        X.temp <- p_viewshed[int.ind,1] + b[int.ind]*r.length[int.ind,1]
        Y.temp <- p_viewshed[int.ind,2] + b[int.ind]*r.length[int.ind,2]

        # Find distance from intersection to end point
        dist_1 <- sqrt((p_animal[nn + 1, 1] - X.temp)^2 +
                         (p_animal[nn + 1, 2] - Y.temp)^2)

        # Calculate time spent in camera
        t_tot <- t_tot + dist_1/speed[nn]

      } else if(sum(int.check) == 2) {
        # Individual passes through camera

        int.ind <- which(int.check == T)

        # X,Y coordinates at intersection
        X.temp <- p_viewshed[int.ind,1] + b[int.ind]*r.length[int.ind,1]
        Y.temp <- p_viewshed[int.ind,2] + b[int.ind]*r.length[int.ind,2]

        # Find distance between points of intersection
        dist_1 <- sqrt((X.temp[1] - X.temp[2])^2 +
                         (Y.temp[1] - Y.temp[2])^2)

        # Calculate time spent in camera
        t_tot <- t_tot + dist_1/speed[nn]
      }
    }
    # Accumulate all x,y coords
    X <- c(X, X.temp)
    Y <- c(Y, Y.temp)

  }
  return(list(cbind(X, Y), t_tot))
}

###################################
# Landscape creator helper
###################################
# Define individual cells
Cell_bias <- function(lscape_data, Index, Dir, Kappa, Road_ID) {
  lscape_data$Direction[Index] <- Dir
  lscape_data$Kappa[Index] <- Kappa
  lscape_data$Road[Index] <- Road_ID

  return(lscape_data)
}


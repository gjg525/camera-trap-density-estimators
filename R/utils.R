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
calc_tri_area <- function(x,y) {
  # Calculate area of triangle using vertices
  tri_area <- 0.5 * abs(x[1]*(y[2] - y[3]) + x[2]*(y[3] - y[1]) + x[3] * (y[1] - y[2]))
}

calc_tri_signs <- function(p1, p2, p3) {
  tri_sign <- (p1[1] - p3[1]) * (p2[2] - p3[2]) - (p2[1] - p3[1]) * (p1[2] - p3[2])
}

# Determine if a point lies within a triangle using where the point lies relative to the triangle half-planes
calc_in_triangle <- function(p_viewshed, point) {

  A_1 <- calc_tri_signs(point, p_viewshed[1,], p_viewshed[2,])
  A_2 <- calc_tri_signs(point, p_viewshed[3,], p_viewshed[1,])
  A_3 <- calc_tri_signs(point, p_viewshed[2,], p_viewshed[3,])

  neg_vals <- (A_1 < 0) || (A_2 < 0) || (A_3 < 0)
  pos_vals <- (A_1 > 0) || (A_2 > 0) || (A_3 > 0)

  in_triangle <- ifelse(!(neg_vals && pos_vals),
                        TRUE,
                        FALSE)
}

# # Determine if a point lies within a triangle using areas (slower)
# calc_in_triangle <- function(p_viewshed, point) {
# 
#   A_1 <- calc_tri_area(t(rbind(p_viewshed[1:2,1], point[1])),
#                         t(rbind(p_viewshed[1:2,2], point[2])))
#   A_2 <- calc_tri_area(t(rbind(p_viewshed[c(1,3),1], point[1])),
#                         t(rbind(p_viewshed[c(1,3),2], point[2])))
#   A_3 <- calc_tri_area(t(rbind(p_viewshed[2:3,1], point[1])),
#                         t(rbind(p_viewshed[2:3,2], point[2])))
# 
#   tri_area <- calc_tri_area(t(p_viewshed[,1]), t(p_viewshed[, 2]))
# 
#   in_triangle <- ifelse((A_1 + A_2 + A_3 - 0.000001) > tri_area,
#                         FALSE,
#                         TRUE)
# }

# Calculate where an animal's trajectory intersects with a camera viewshed and the time spent within viewshed
calc_intersects <- function(p_viewshed, p_animal, speed, t) {
  X <- c()
  Y <- c()
  t_stay <- c()
  in_cam_all <- c()
  t_in_all <- c()
  t_out_all <- c()
  # Capture encounters as an individual leaves the viewshed
  encounter <- 0
  num_encounters <- 0
  
  # Get lengths of all edges
  r.length <- rbind(p_viewshed[2,] - p_viewshed[1,],
                    p_viewshed[3,] - p_viewshed[2,],
                    p_viewshed[1,] - p_viewshed[3,])

  for (nn in 1:(nrow(p_animal)-1)) {
    X.temp <- c()
    Y.temp <- c()
    t_tot <- 0 # time spent in camera
    t_in <- NA
    t_out <- NA

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
        # t_tot <- t_tot + dist_1/speed[nn]
        t_tot <- dist_1/speed[nn]
        
        # Calculate enter/exit times
        t_in <- NA
        t_out <- t[nn] + t_tot
        
        num_encounters <- num_encounters + 1
      } else {
        # Individual remains in camera

        # Find distance from start point to end point
        dist_1 <- sqrt((p_animal[nn + 1, 1] - p_animal[nn, 1])^2 +
                         (p_animal[nn + 1, 2] - p_animal[nn, 2])^2)

        # Calculate time spent in camera
        # t_tot <- t_tot + dist_1/speed[nn]
        t_tot <- dist_1/speed[nn]
        
        # Encounter if individual remains in viewshed on final step
        if (t[nn] == t[length(t) - 1]) {
          num_encounters <- num_encounters + 1
        }
        
        # # Calculate enter/exit times
        # t_in <- NA
        # t_out <- NA
        
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
        # t_tot <- t_tot + dist_1/speed[nn]
        t_tot <- dist_1/speed[nn]
        
        # Find distance from start point to intersection
        dist_2 <- sqrt((p_animal[nn, 1] - X.temp)^2 +
                         (p_animal[nn, 2] - Y.temp)^2)
        
        # Calculate enter/exit times
        t_in <- t[nn] + dist_2/speed[nn]
        # t_out <- NA
        
        # Encounter if individual remains in viewshed on final step
        if (t[nn] == t[length(t) - 1]) {
          num_encounters <- num_encounters + 1
        }
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
        # t_tot <- t_tot + dist_1/speed[nn]
        t_tot <- dist_1/speed[nn]
        
        # Find distance from start point to first intersection
        dist_2 <- sqrt((p_animal[nn, 1] - X.temp[1])^2 +
                         (p_animal[nn, 2] - Y.temp[1])^2)
        
        # Calculate enter/exit times
        t_in <- t[nn] + dist_2/speed[nn]
        t_out <- t_in + t_tot
        
        num_encounters <- num_encounters + 1
      }
    }
    
    if (length(t_in) > 1){
      print(t[nn])
    }
    # Accumulate all values
    X <- c(X, X.temp)
    Y <- c(Y, Y.temp)
    t_stay <- c(t_stay, t_tot)
    in_cam_all <- c(in_cam_all, in_cam)
    encounter <- c(encounter, num_encounters)
    t_in_all <- c(t_in_all, t_in)
    t_out_all <- c(t_out_all, t_out)
    
  }

  return(list(cbind(X, Y), t_stay, in_cam_all, encounter, t_in_all, t_out_all))
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


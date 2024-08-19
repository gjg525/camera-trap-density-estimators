calc_local_area <- function(cover_all, xs, ys, ts) {

  nbuf <- 2 * dx
  maxx <- max(xs)
  minx <- min(xs)
  lx <- maxx - minx
  maxy <- max(ys)
  miny <- min(ys)
  ly <- maxy - miny
  jcolmin <- floor((minx) / dx) + 1 - nbuf
  irowmin <- floor((miny) / dx) + 1 - nbuf
  jcolmax <- floor((maxx) / dx) + 1 + nbuf
  irowmax <- floor((maxy) / dx) + 1 + nbuf
  idxx <- seq(jcolmin, jcolmax)
  idxy <- seq(irowmin, irowmax)
  Area <- length(idxx) * length(idxy)
  cover_small <- cover_all[idxy, idxx]
  
  return(list(cover_small, Area))
}

load(file = paste0("Sim_results/Sim_Original_all_vars.RData"))

# Threshold parms
n_thresh <- 20
n_win <- 100

num_runs <- 1000

# Define number of clumps
num.clumps <- 100

# Define clump sizes for every clump
clump_sizes <- rep(1,num.clumps)
nind <- sum(clump_sizes)

# Landscape parms
q <- 30^2             # Number grid cells
bounds <- c(0, q^0.5) # Sampling area boundaries
t.steps <- 500        # Number of time steps
dt <- 1               # Time step size

# Grid cell lengths
dx <- (bounds[2]-bounds[1])/q^0.5
dy <- (bounds[2]-bounds[1])/q^0.5

xs <- animalxy.all %>% 
  dplyr::filter(Animal_ID == 1) %>% 
  dplyr::pull(x)

ys <- animalxy.all %>% 
  dplyr::filter(Animal_ID == 1) %>% 
  dplyr::pull(y)

ts <- animalxy.all %>% 
  dplyr::filter(Animal_ID == 1) %>% 
  dplyr::pull(t)

dts <- diff(ts)

dsqrd <- diff(xs)^2 + diff(ys)^2

mubar <- mean(.25 * dsqrd / (max(.25, dts)), na.rm = T)

# find cover types visited
js <- floor((xs) / dx) + 1  # row, column indices for visited pixels
is <- floor((ys) / dx) + 1  
cover_visited <- rep(0, length(is))

cover_all <- matrix(lscape_speeds$Speed, nrow = 30, ncol = 30, byrow = T)
for (j in 1:length(is)) {
  cover_visited[j] <- cover_all[is[j], js[j]]
}

nobs <- length(cover_visited)  # total number of observations

# find relevant cover types with enough visitation
visited <- unique(cover_visited)

for (j in 1:length(visited)) {
  nj[j] <- sum(cover_visited == visited[j])
  # visited[j] <- (nj > n_thresh) * visited[j]
}

visited <- visited[nj > n_thresh] 

local_area <- calc_local_area(cover_all, xs, ys, ts)
cover_small <- local_area[[1]]
Area <- local_area[[2]]

if (nobs > n_win) {
  # sometimes there are not enough observations, so can't do the windowed calculation
  
  dts2 <- 0.5 * (c(dts, 0) + c(0, dts))  # approx delta ts for each position
  DjT <- matrix(0, nrow = length(visited), ncol = nobs - n_win)
  PjT <- DjT  # Pj = Dj; Dj = DjT;
  Aj <- rep(0, length(visited))
  
  for (j in 1:length(visited)) {
    Aj[j] <- sum(sum(cover_small == visited[j]))  # modulo dx^2
    for (nt in 1:(nobs - n_win)) {
      iwin <- nt:(nt + n_win - 1) # window indices
      Twin <- sum(dts2[iwin], na.rm = TRUE) # total time of window
      PjT[j, nt] <- sum((cover_visited[iwin] == visited[j]) * dts2[iwin]) / max(1, Twin)
      DjT[j, nt] <- (PjT[j, nt] > 0) * 
        mubar * 
        Aj[j] / 
        Area / 
        max(0.02 / Twin, PjT[j, nt])
    }
  }
  
  # in julian days after collaring, the times will be:
  jd <- ts[(n_win / 2 + 1):(nobs - n_win / 2)] - 
    floor(ts[n_win / 2] / 365) * 365
} else {
  # not enough observations
  jd <- 0  # jd will be zero and DjT empty if not enough obs
}

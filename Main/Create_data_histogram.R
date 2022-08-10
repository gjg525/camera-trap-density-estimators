# # Main script (run aux_functions, MLE_functions, and Initializations prior)

## NOTE tidyr may cause errors in extract function
library(tidyverse)
# library(pbapply) 
library(msm) 
library(MASS)
# library(Directional)
library(raster)
# library(fields)
library(RColorBrewer)
# library(rasterVis)
library(truncnorm)
# library(mvtnorm)
# library(gridExtra)
library(fBasics)
library(coda)
# library(Rlab)
# library(fitdistrplus)
# library(patchwork)
library(crayon)
# library(GMCM)
library(Rfast)
library(dplyr)
library(beepr)
options(error = beep)

gc() # take out the trash
set.seed(1)

fig_dir <- "C:/Users/guen.grosklos/Google Drive/Missoula_postdoc/Code/All_models/"

source("./Main/MCMC_functions.R")
source("./Main/MLE_functions.R")


###########################
# Initialize
###########################
# Define camera sample variations
# Cam sample design (1: random, 2: 80% slow, 3: 80% medium, 4: 80% fast)
cs.all <- 1
cam.dist.labels <- c("random","slow","med","fast")

# Define landscape variations
# 1: all slow, 2: all medium, 3: all fast, 4: equal slow, medium, fast 5: 80% fast
lv.all <- 4
lv.labels <- c("slow_lscape_all","med_lscape_all","fast_lscape_all","","fast_lscape")

# Variable parms
num.clumps <- 100
clump.size <- rep(1,num.clumps)
nind<-sum(clump.size)
num.runs <- 100 # number of repeated runs

# Correlated random walk parameter
# 0 is uncorrelated random walk, inf is ideal gas model (5 is good correlation)
corr.walk.kappa <- 5

# # Number cameras
ncam <-  900

# # Legend labels
leg1<-c("EEDE", "REST", "TTE", "MCT", "STE")
leg.props<-c("EEDE (MLE)", "REST (MLE)", "TTE (MLE)", "MCT (MLE)",
             "EEDE (MCMC)", "REST (MCMC)", "TTE (MCMC)", "MCT (MCMC)")

# Agent based model parms
q <- 30^2   # number grid cells
bounds <- c(0, q^0.5)
t.steps <- 200
dt <- 1 # time step size
dx <- (bounds[2]-bounds[1])/q^0.5 # space step size in x
dy <- (bounds[2]-bounds[1])/q^0.5 # space step size in y 
clump.rad <- dx/2 # tightness of clumping behavior

tot.A <- (bounds[2]-bounds[1])^2    # total area
cam.A <- ((bounds[2]-bounds[1])/q^0.5)^2 # For now, camera area is same as cell area

# Fitting parms for MCMC
n.iter <- 5000
burn.in <- 2000

# moving speed bounds for slow, medium, fast speeds
# speed.bounds <- dx/rbind(c(5.1,4.9),c(3.1,2.9),c(1.1,0.9))
speed.bounds <- dx/rbind(c(10.1,9.9),c(2.6,2.4),c(0.6,0.4))

# indices for all covariates including intercept in spatial.covariates
# Assumes 3 covariate types (slow, medium, fast)
covariates.index <- c(0,rep(1,3))

# Data collection parms
# Occasion length for TTE
JJ <- 20  
# TTE and staying time censor
t.censor <- JJ
num.occ <- t.steps/JJ
# STE censor
STE.censor <- ncam*cam.A

#########################################
# Initialize landscape
#########################################
# Define slow, medium, and fast movement speed cells
# z.slow <- matrix(0,q^0.5,q^0.5)
# z.med <- matrix(0,q^0.5,q^0.5)
# z.fast <- matrix(runif(q,speed.bounds[3,1],speed.bounds[3,2]),q^0.5,q^0.5)

# Initialize landscape with equal slow, medium, and fast cells
# num.slow.inds <- round(q/3)
# num.med.inds <- round(q/3)

# # Calculate the number of fast cells
# num.fast.inds <- q-num.slow.inds-num.med.inds
# sample.inds <- sample(1:q,num.slow.inds+num.med.inds,replace = F)
# fast.inds <- 1:q
# fast.inds<-fast.inds[-sample.inds]
# slow.inds <- sample.inds[1:num.slow.inds]
# med.inds <- sample.inds[(num.slow.inds+1):(num.slow.inds+num.med.inds)]
# z.fast[sample.inds] <- 0
# z.slow[slow.inds] <- runif(num.slow.inds,speed.bounds[1,1],speed.bounds[1,2])
# z.med[med.inds] <- runif(num.med.inds,speed.bounds[2,1],speed.bounds[2,2])




# Define slow, medium, and fast movement speed cells
z.slow <- matrix(0,q^0.5,q^0.5)
z.med <- matrix(0,q^0.5,q^0.5)
z.fast <- matrix(0,q^0.5,q^0.5)

sample.inds <- sample(1:q,q,replace = F)

num.slow.inds <- round(q/3)
num.med.inds <- round(q/3)
num.fast.inds <- q-num.slow.inds-num.med.inds
slow.inds <- sample.inds[1:num.slow.inds]
med.inds <- sample.inds[(num.slow.inds+1):(num.slow.inds+num.med.inds)]
fast.inds <- sample.inds[(num.slow.inds+num.med.inds+1):q]
z.slow[slow.inds] <- runif(num.slow.inds,speed.bounds[1,1],speed.bounds[1,2])
z.med[med.inds] <- runif(num.med.inds,speed.bounds[2,1],speed.bounds[2,2])
z.fast[fast.inds] <- runif(num.fast.inds,speed.bounds[3,1],speed.bounds[3,2])


## Create rasterstack object for covariates
spatial.covariates <- raster(, nrows = q^0.5, ncols = q^0.5, xmn = 0, xmx = 1,
                             ymn = 0, ymx = 1, crs = NA)
spatial.covariates$cell <- c(1:q)
spatial.covariates$slow <- c(z.slow)/dx
spatial.covariates$med <- c(z.med)/dx
spatial.covariates$fast <- c(z.fast)/dx

Z <- as.matrix(spatial.covariates[[which(covariates.index==1)]])

# Define movement speeds for each cell
v.abm <- matrix(z.slow+z.med+z.fast,q^0.5,q^0.5)

# Define staying time
t.stay.abm <- 4/(dx*v.abm)

# Create cumulative distribution to distribute population according to stay time
t.stay.scale <- t.stay.abm/sum(t.stay.abm)
t.stay.cdf <- rep(0,q+1)
for (ii in 1:q){
  t.stay.cdf[ii+1] <- sum(t.stay.scale[1:ii]) 
}

#########################################
# Random walks for nind animals
#########################################
# Run agent-based model
source("./Main/ABM_sims.R")

# # 1-cell sized cameras, randomly selected
cam.samps <- sample(1:q,ncam, replace=FALSE)

# Make sure cam.samps adds up correctly
if(length(cam.samps)!=ncam){
  stop("camera samples not equal to number cameras")
}

# Camera statistics
cam.counts <- u.abm.all[cam.samps] # num counts at each camera for each time
cam.counts.sum <- rowSums(cam.counts) # num counts at each camera 

# indices for each covariate type
cam.slow <- which(cam.samps %in% slow.inds)
cam.med <- which(cam.samps %in% med.inds)
cam.fast <- which(cam.samps %in% fast.inds)
cam.props <- c(length(cam.slow),length(cam.med),length(cam.fast))/ncam
cam.props.rounds <- round(cam.props*100)/100

# number of counts from cameras for each covariate type
cam.counts.slow <- cam.counts[cam.slow,]
cam.counts.med <- cam.counts[cam.med,]
cam.counts.fast <- cam.counts[cam.fast,]
cam.counts.props <- c(sum(cam.counts.slow),sum(cam.counts.med),sum(cam.counts.fast))/sum(cam.counts)

# number of counts across whole landscape for each covariate type
slow.counts <- sum(u.abm.all[slow.inds])/t.steps
med.counts <- sum(u.abm.all[med.inds])/t.steps
fast.counts <- sum(u.abm.all[fast.inds])/t.steps

# Check that abm distributions match predicted distributions based on staying time
stay.time.distribution <- c(mean(t.stay.abm[slow.inds]),
                            mean(t.stay.abm[med.inds]),
                            mean(t.stay.abm[fast.inds]))
stay.time.distribution.scale <-stay.time.distribution/sum(stay.time.distribution)
abm.distribution <- c(slow.counts,med.counts,fast.counts)/c(num.slow.inds,num.med.inds,num.fast.inds)
abm.distribution.scale <- abm.distribution/sum(abm.distribution)

################################
# Collect data
################################
source("./Main/collect_data.R")

row_count_df <- data.frame(counts = c(c(cam.counts.sum[cam.slow]),
                                  c(cam.counts.sum[cam.med]),
                                  c(cam.counts.sum[cam.fast])),
                       speed = rep(c("Slow","Medium","Fast"),
                                   each = length(c(row_counts[cam.slow]))))


ggplot(row_count_df, aes(x=counts, fill = speed)) +
  geom_density(position = "identity", alpha = 0.4, adjust = 3) +
  labs(x = "Count", y = "Density", fill="Landscape Type") +
  xlim(c(0,100)) +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        legend.text=element_text(size=20),
        legend.title = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
# ggsave(file = paste(fig_dir,"figs/Counts_dat_histogram.eps", sep = ""), device = cairo_ps)

enc_df <- data.frame(encounter = c(c(num.encounters.dat[cam.slow]),
                                c(num.encounters.dat[cam.med]),
                                c(num.encounters.dat[cam.fast])),
                     speed = rep(c("Slow","Medium","Fast"),each = 300))

ggplot(enc_df, aes(x=encounter, fill = speed)) +
  # geom_histogram(position = "identity", alpha = 0.2, bins = 5)
  geom_density(position = "identity", alpha = 0.4, adjust = 2) +
  labs(x = "Number of Encounters", y = "Density", fill="Landscape Type") +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        legend.text=element_text(size=20),
        legend.title = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
# ggsave(file = paste(fig_dir,"figs/Encounter_dat_histogram.eps", sep = ""), device = cairo_ps)


stay_df <- data.frame(Stay_time = c(c(t.staying.dat[cam.slow,]),
                                    c(t.staying.dat[cam.med,]),
                                    c(t.staying.dat[cam.fast,])),
                      speed = rep(c("Slow","Medium","Fast"),
                                  each = c(length(c(t.staying.dat[cam.slow,])))))

ggplot(stay_df, aes(x=Stay_time, fill = speed)) +
  # geom_histogram(position = "identity", alpha = 0.4, bins = 20)
  geom_density(position = "identity", alpha = 0.4, adjust = 5) +
  xlim(c(0,20)) +
  labs(x = "Staying Time", y = "Density", fill="Landscape Type") +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        legend.text=element_text(size=20),
        legend.title = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
# ggsave(file = paste(fig_dir,"figs/Stay_time_dat_histogram.eps", sep = ""), device = cairo_ps)

TTE_df <- data.frame(TTE = c(c(TTE.dat[cam.slow,]),
                                c(TTE.dat[cam.med,]),
                                c(TTE.dat[cam.fast,])),
                      speed = rep(c("Slow","Medium","Fast"),
                                  each = length(c(TTE.dat[cam.slow,]))))

ggplot(TTE_df, aes(x=TTE, fill = speed)) +
  # geom_histogram(position = "identity", alpha = 0.2, bins = 5)
  geom_density(position = "identity", alpha = 0.4, adjust = 2) +
  labs(x = "Time to Encounter", y = "Density", fill="Landscape Type") +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        legend.text=element_text(size=20),
        legend.title = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
# ggsave(file = paste(fig_dir,"figs/TTE_dat_histogram.eps", sep = ""), device = cairo_ps)


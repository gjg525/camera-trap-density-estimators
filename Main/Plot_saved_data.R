library(RColorBrewer)
library(GMCM)

# Read saved files for plotting
fig_dir <- "C:/Users/guen.grosklos/Google Drive/Missoula_postdoc/Code/All_models/"

# Initializations
nind <- 100
leg1<-c("EEDE", "REST", "TTE", "MCT", "STE")
leg.props<-c("EEDE", "REST", "TTE", "MCT")

# Cam sample design (1: random, 2: 80% slow, 3: 80% medium, 4: 80% fast)
cs.all <- 4
cam.dist.labels <- c("random","slow","med","fast")
cam.dist.labels.caps <- c("Random","Slow","Medium","Fast")

# Define landscape variations
# 1: all slow, 2: all medium, 3: all fast, 4: equal slow, medium, fast 5: 80% fast
lv.all <- 5
lv.labels <- c("_slow_lscape_all","_med_lscape_all","_fast_lscape_all","","_fast_lscape")

for (cam.dist.set in cs.all){
  for (lscape_var in lv.all){
    
# Labels for saving figures
means_label <- paste(cam.dist.labels[cam.dist.set],"_cams_means",lv.labels[lscape_var], sep = "")
props_label <- paste(cam.dist.labels[cam.dist.set],"_cams_props",lv.labels[lscape_var], sep = "")


# # # Load .csv files
# D.all.Means.mat <- as.matrix(read.csv(paste(fig_dir,"sim_data/",means_label,".csv", sep = "")))
# D.all.Sds.mat <- as.matrix(read.csv(paste(fig_dir,"sim_data/",means_label,"_sds.csv", sep = "")))
D.all.Means.mat.all <- as.matrix(read.csv(paste(fig_dir,"sim_data/all_",means_label,".csv", sep = "")))
D.all.Sds.mat.all <- as.matrix(read.csv(paste(fig_dir,"sim_data/all_",means_label,"_SD.csv", sep = "")))
all.props.Means <- as.matrix(read.csv(paste(fig_dir,"sim_data/",props_label,".csv", sep = "")))
all.props.Sds <- as.matrix(read.csv(paste(fig_dir,"sim_data/",props_label,"_sds.csv", sep = "")))

all.props.Means <- matrix(as.numeric(all.props.Means[c(1,3,5,7,9),2:4]),nrow = 5, ncol = 3)
all.props.Sds <- matrix(as.numeric(all.props.Sds[c(1,3,5,7,9),2:4]),nrow = 5, ncol = 3)
rownames(all.props.Means) <- c("ABM (Truth)",leg.props)
colnames(all.props.Means) <- c("Slow","Medium","fast")

####################################
# Calculate Summaries
####################################
D.all.Means <- colMeans(D.all.Means.mat.all[,2:ncol(D.all.Means.mat.all)])
D.all.Sds <- apply(D.all.Means.mat.all[,2:ncol(D.all.Means.mat.all)], 2, sd)

D.all.MLE.Means <- D.all.Means[1:5]
D.all.MLE.Sds <- D.all.Sds[1:5]
D.all.MLE.cov.Means <- D.all.Means[6:10]
D.all.MLE.cov.Sds <- D.all.Sds[6:10]
D.all.MCMC.Means <- D.all.Means[11:15]
D.all.MCMC.Sds <- D.all.Sds[11:15]
D.all.MCMC.cov.Means <- D.all.Means[16:20]
D.all.MCMC.cov.Sds <- D.all.Sds[16:20]

# Standard deviation results
SD.all.Means <- colMeans(D.all.Sds.mat.all[,2:ncol(D.all.Sds.mat.all)])
SD.all.Sds <- apply(D.all.Sds.mat.all[,2:ncol(D.all.Sds.mat.all)], 2, sd)

SD.all.MCMC.Means <- SD.all.Means[11:15]
SD.all.MCMC.Sds <- SD.all.Sds[11:15]
SD.all.MCMC.cov.Means <- SD.all.Means[16:20]
SD.all.MCMC.cov.Sds <- SD.all.Sds[16:20]

if (any(lv.all == 1:3)) {
  # setEPS()
  # postscript(paste(fig_dir,"figs/",means_label,".eps", sep = ""),width=8,height=5)
  op <- par(mar=c(5, 6, 4, 2) + 0.1)
  cam.props.label <- paste("Camera Bias: ", cam.dist.labels.caps[cam.dist.set], sep = "")
  plot(seq(0.9,4.9,by=1), 
       D.all.MCMC.Means, 
       ylim=c(min(nind, D.all.Means[11:20] - D.all.Sds[11:20],na.rm=T),
              max(nind, D.all.Means[11:20] + D.all.Sds[11:20],na.rm=T)),
       xlim=c(0.5,5.5),xlab="Method", 
       ylab=paste(cam.props.label,"\n Mean Estimates"), 
       pch=16, col="black", cex=1.8, xaxt = "n", cex.lab = 1.5, cex.axis = 1.3)
  points(seq(1.1,5.1,by=1),D.all.MCMC.cov.Means, col="black", pch=1,cex=1.8)
  arrows(x0=seq(0.9,4.9,by=1), y0=D.all.MCMC.Means-D.all.MCMC.Sds, 
         x1=seq(0.9,4.9,by=1), y1=D.all.MCMC.Means+D.all.MCMC.Sds, 
         code=3, angle=90, length=0.1, col="black", lwd=2)
  arrows(x0=seq(1.1,5.1,by=1), y0=D.all.MCMC.cov.Means-D.all.MCMC.cov.Sds, 
         x1=seq(1.1,5.1,by=1), y1=D.all.MCMC.cov.Means+D.all.MCMC.cov.Sds, 
         code=3, angle=90, length=0.1, col="black", lwd=2)
  lines(c(0.3,5.7), c(nind,nind), type = "l", lty = 2, lwd = 3.5)
  axis(1, at = c(1:5), labels = leg1, cex.axis = 1.3)
  legend("topright", c("Non-Covariate","Covariate"),
         pch=c(16,1), col=c("black","black"))
  par(op)
  # dev.off()
  
} else{
  # setEPS()
  # postscript(paste(fig_dir,"figs/",means_label,".eps", sep = ""),width=8,height=5)
  op <- par(mar=c(5, 6, 4, 2) + 0.1)
  cam.props.label <- paste("Camera Bias: ", cam.dist.labels.caps[cam.dist.set], sep = "")
  plot(seq(0.9,4.9,by=1), 
       D.all.MCMC.Means, 
       ylim=c(min(nind, D.all.Means[11:20] - D.all.Sds[11:20],na.rm=T),
              max(nind, D.all.Means[11:20] + D.all.Sds[11:20],na.rm=T)),
       xlim=c(0.5,5.5),xlab="Method", 
       ylab=paste(cam.props.label,"\n Mean Estimates"), 
       pch=16, col="black", cex=1.8, xaxt = "n", cex.lab = 1.5, cex.axis = 1.3)
  points(seq(1.1,5.1,by=1),D.all.MCMC.cov.Means, col="black", pch=1,cex=1.8)
  arrows(x0=seq(0.9,4.9,by=1), y0=D.all.MCMC.Means-D.all.MCMC.Sds, 
         x1=seq(0.9,4.9,by=1), y1=D.all.MCMC.Means+D.all.MCMC.Sds, 
         code=3, angle=90, length=0.1, col="black", lwd=2)
  arrows(x0=seq(1.1,5.1,by=1), y0=D.all.MCMC.cov.Means-D.all.MCMC.cov.Sds, 
         x1=seq(1.1,5.1,by=1), y1=D.all.MCMC.cov.Means+D.all.MCMC.cov.Sds, 
         code=3, angle=90, length=0.1, col="black", lwd=2)
  lines(c(0.3,5.7), c(nind,nind), type = "l", lty = 2, lwd = 3.5)
  axis(1, at = c(1:5), labels = leg1, cex.axis = 1.3)
  legend("topright", c("Non-Covariate","Covariate"),
         pch=c(16,1), col=c("black","black"))
  par(op)
  # dev.off()
  
  # setEPS()
  # postscript(paste(fig_dir,"figs/",props_label,".eps", sep = ""),width=8,height=5)
  bar.p <- barplot(all.props.Means, 
                   names = c("Slow", "Medium", "Fast"), 
                   beside = T, 
                   legend.text=T,
                   # col = c("black",brewer.pal(4, "Set2"),brewer.pal(4, "Set2")),
                   col = c("black",brewer.pal(4, "Set2")),
                   # border = c("black",brewer.pal(8, "Paired")),
                   # alpha = c(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5),
                   # col=c("black","blue","red","green","violet","blue","red","green","violet"),
                   density=c(300,300,300,300,300),
                   angle=c(0,0,0,0,0),
                   ylim=c(0,1) , ylab="Proportion Occupied",xlab = "Movement Speeds",
                   args.legend = list(x = "topright",
                                      inset = c( 0.035, 0)),
                   cex=1.3, 
                   cex.lab = 1.5, 
                   cex.axis = 1.3)
  arrows(bar.p,all.props.Means+all.props.Sds, 
         bar.p, 
         all.props.Means-all.props.Sds, 
         angle=90, 
         code=3, 
         length=0.05)
  # dev.off()
  
  # setEPS()
  # postscript(paste(fig_dir,"figs/SD_",means_label,".eps", sep = ""),width=8,height=5)
  op <- par(mar=c(5, 6, 4, 2) + 0.1)
  cam.props.label <- paste("Camera Bias: ", cam.dist.labels.caps[cam.dist.set], sep = "")
  plot(seq(0.9,4.9,by=1), 
       SD.all.MCMC.Means, 
       ylim=c(min(SD.all.Means[11:20] - SD.all.Sds[11:20],na.rm=T),
              max(SD.all.Means[11:20] + SD.all.Sds[11:20],na.rm=T)),
       xlim=c(0.5,5.5),xlab="Method", 
       ylab=paste(cam.props.label,"\n SD Estimates"), 
       pch=16, col="black", cex=1.8, xaxt = "n", cex.lab = 1.5, cex.axis = 1.3)
  points(seq(1.1,5.1,by=1),SD.all.MCMC.cov.Means, col="black", pch=1,cex=1.8)
  arrows(x0=seq(0.9,4.9,by=1), y0=SD.all.MCMC.Means-SD.all.MCMC.Sds, 
         x1=seq(0.9,4.9,by=1), y1=SD.all.MCMC.Means+SD.all.MCMC.Sds, 
         code=3, angle=90, length=0.1, col="black", lwd=2)
  arrows(x0=seq(1.1,5.1,by=1), y0=SD.all.MCMC.cov.Means-SD.all.MCMC.cov.Sds, 
         x1=seq(1.1,5.1,by=1), y1=SD.all.MCMC.cov.Means+SD.all.MCMC.cov.Sds, 
         code=3, angle=90, length=0.1, col="black", lwd=2)
  axis(1, at = c(1:5), labels = leg1, cex.axis = 1.3)
  legend("topright", c("Non-Covariate","Covariate"),
         pch=c(16,1), col=c("black","black"))
  par(op)
  # dev.off()

}
  }
}

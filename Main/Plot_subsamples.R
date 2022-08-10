
# file directory for saving
fig_dir <- "C:/Users/guen.grosklos/Google Drive/Missoula_postdoc/Code/All_models/"

means_label <- "random_cams_means"

results_out <- data.frame(matrix(NA, nrow = 0, ncol = 7))
IDs <- c("mean", "mean_cov", "sd", "sd_cov", "CI", "CI_cov")
ncam_all <- c(5, 10, 20, 50, 100)

# mcmc.means <- matrix(NA, nrow = 5, ncol = 5)
# colnames(mcmc.means) <- c("EEDE", "REST", "TTE", "MCT", "STE")
# mcmc.means.cov <- matrix(NA, nrow = 5, ncol = 5)
# colnames(mcmc.means.cov) <- c("EEDE", "REST", "TTE", "MCT", "STE")
# mcmc.sds <- matrix(NA, nrow = 5, ncol = 5)
# colnames(mcmc.sds) <- c("EEDE", "REST", "TTE", "MCT", "STE")
# mcmc.sds.cov <- matrix(NA, nrow = 5, ncol = 5)
# colnames(mcmc.sds.cov) <- c("EEDE", "REST", "TTE", "MCT", "STE")
# mcmc.CI <- matrix(NA, nrow = 5, ncol = 5)
# colnames(mcmc.CI) <- c("EEDE", "REST", "TTE", "MCT", "STE")
# mcmc.CI.cov <- matrix(NA, nrow = 5, ncol = 5)
# colnames(mcmc.CI.cov) <- c("EEDE", "REST", "TTE", "MCT", "STE")

for (nca in 1:length(ncam_all)) {
  # # Load .csv files
  D.all <- as.matrix(read.csv(paste(fig_dir,"sim_data/",means_label,nca, "_cams.csv", sep = "")))
  SD.all <- as.matrix(read.csv(paste(fig_dir,"sim_data/",means_label,nca, "_cams_SD.csv", sep = "")))
  
  D.all[is.infinite(D.all)] <- NA
  names(D.all)<−NULL
  
  D.all.mcmc <- D.all[,12:16]
  D.all.mcmc.cov <- D.all[,17:21]
  SD.all.mcmc <- SD.all[,12:16]
  SD.all.mcmc.cov <- SD.all[,17:21]
  
  results_out <- rbind(results_out, 
                       c(ncam_all[nca], IDs[1], colMeans(D.all.mcmc, na.rm = T)),
                       c(ncam_all[nca], IDs[2], colMeans(D.all.mcmc.cov, na.rm = T)),
                       c(ncam_all[nca], IDs[3], apply(D.all.mcmc, 2, sd, na.rm = T)),
                       c(ncam_all[nca], IDs[4], apply(D.all.mcmc.cov, 2, sd, na.rm = T)),
                       c(ncam_all[nca], IDs[5], colMeans(SD.all.mcmc, na.rm = T)),
                       c(ncam_all[nca], IDs[6], colMeans(SD.all.mcmc.cov, na.rm = T))
                       )
  
  mcmc.means[nca,] <- colMeans(D.all.mcmc, na.rm = T)
  mcmc.means.cov[nca,] <- colMeans(D.all.mcmc.cov, na.rm = T)
  mcmc.sds[nca,] <- apply(D.all.mcmc, 2, sd, na.rm = T)
  mcmc.sds.cov[nca,] <- apply(D.all.mcmc.cov, 2, sd, na.rm = T)
  
  mcmc.CI[nca,] <- colMeans(SD.all.mcmc, na.rm = T)
  mcmc.CI.cov[nca,] <- colMeans(SD.all.mcmc.cov, na.rm = T)
  
}
colnames(results_out) <- c("num_cams", "ID", "EEDE", "REST", "TTE", "MCT", "STE")

ggplot(data = results_out, aes(x = num_cams, y = EEDE)) +
  geom_line()

plot(mcmc.means[,2:5])


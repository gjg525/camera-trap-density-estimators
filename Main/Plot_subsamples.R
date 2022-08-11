
# file directory for saving
fig_dir <- "C:/Users/guen.grosklos/Google Drive/Missoula_postdoc/Code/All_models/"

means_label <- "random_cams_means"

IDs <- c("mean", "mean_cov", "sd", "sd_cov", "CI", "CI_cov")
Model <- c("EEDE", "REST", "TTE", "MCT", "STE")
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

results_out <- data.frame(matrix(NA, nrow = 0, ncol = 4))
for (nca in 1:length(ncam_all)) {
  # # Load .csv files
  D.all <- as.matrix(read.csv(paste(fig_dir,"sim_data/",means_label,nca, "_cams.csv", sep = "")))
  SD.all <- as.matrix(read.csv(paste(fig_dir,"sim_data/",means_label,nca, "_cams_SD.csv", sep = "")))
  
  D.all[is.infinite(D.all)] <- NA
  names(D.all)<- NULL
  
  D.all.mcmc <- D.all[,12:16]
  D.all.mcmc.cov <- D.all[,17:21]
  SD.all.mcmc <- SD.all[,12:16]
  SD.all.mcmc.cov <- SD.all[,17:21]
  
  results_out <- rbind(results_out, 
                       cbind(ncam_all[nca], IDs[1], matrix(Model), matrix(colMeans(D.all.mcmc, na.rm = T))),
                       cbind(ncam_all[nca], IDs[2], matrix(Model), matrix(colMeans(D.all.mcmc.cov, na.rm = T))),
                       cbind(ncam_all[nca], IDs[3], matrix(Model), matrix(apply(D.all.mcmc, 2, sd, na.rm = T))),
                       cbind(ncam_all[nca], IDs[4], matrix(Model), matrix(apply(D.all.mcmc.cov, 2, sd, na.rm = T))),
                       cbind(ncam_all[nca], IDs[5], matrix(Model), matrix(colMeans(SD.all.mcmc, na.rm = T))),
                       cbind(ncam_all[nca], IDs[6], matrix(Model), matrix(colMeans(SD.all.mcmc.cov, na.rm = T)))
                       )
  
  # mcmc.means[nca,] <- colMeans(D.all.mcmc, na.rm = T)
  # mcmc.means.cov[nca,] <- colMeans(D.all.mcmc.cov, na.rm = T)
  # mcmc.sds[nca,] <- apply(D.all.mcmc, 2, sd, na.rm = T)
  # mcmc.sds.cov[nca,] <- apply(D.all.mcmc.cov, 2, sd, na.rm = T)
  # 
  # mcmc.CI[nca,] <- colMeans(SD.all.mcmc, na.rm = T)
  # mcmc.CI.cov[nca,] <- colMeans(SD.all.mcmc.cov, na.rm = T)
  
}
# colnames(results_out) <- c("num_cams", "ID", "EEDE", "REST", "TTE", "MCT", "STE")
colnames(results_out) <- c("num_cams", "ID", "Model", "Est")

results_out %>% 
  dplyr::filter(ID == "mean") %>% 
ggplot(aes(x = num_cams, y = Est, color = Model)) +
  geom_point()



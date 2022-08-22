library(tidyverse)

# file directory for saving
fig_dir <- "C:/Users/guen.grosklos/Google Drive/Missoula_postdoc/Code/All_models/"

means_label <- "random_cams_means"

IDs <- c("mean", "mean_cov", "sd", "sd_cov", "SD", "SD_cov", "SD_sd", "SD_cov_sd")
Model <- c("EEDE", "REST", "TTE", "MCT", "STE")
ncam_all <- c(5, 10, 20, 50, 100)
fig_colors <- c("#2ca25f", "#fc8d59", "#67a9cf", "#f768a1", "#bae4b3", "#fed98e")

results_out <- data.frame(matrix(NA, nrow = 0, ncol = 4))
for (nca in 1:length(ncam_all)) {
  # # Load .csv files
  D.all <- as.matrix(read.csv(paste(fig_dir,"sim_data/",means_label,nca, "_cams.csv", sep = "")))
  SD.all <- as.matrix(read.csv(paste(fig_dir,"sim_data/",means_label,nca, "_cams_SD.csv", sep = "")))
  
  D.all[is.infinite(D.all)] <- NA
  names(D.all)<- NULL
  # Omit estimates that didn't converge
  SD.all[SD.all > 500] <- NA
  
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
                       cbind(ncam_all[nca], IDs[6], matrix(Model), matrix(colMeans(SD.all.mcmc.cov, na.rm = T))),
                       cbind(ncam_all[nca], IDs[7], matrix(Model), matrix(apply(SD.all.mcmc, 2, sd, na.rm = T))),
                       cbind(ncam_all[nca], IDs[8], matrix(Model), matrix(apply(SD.all.mcmc.cov, 2, sd, na.rm = T)))
  )
  
  # mcmc.means[nca,] <- colMeans(D.all.mcmc, na.rm = T)
  # mcmc.means.cov[nca,] <- colMeans(D.all.mcmc.cov, na.rm = T)
  # mcmc.sds[nca,] <- apply(D.all.mcmc, 2, sd, na.rm = T)
  # mcmc.sds.cov[nca,] <- apply(D.all.mcmc.cov, 2, sd, na.rm = T)
  # 
  # mcmc.SD[nca,] <- colMeans(SD.all.mcmc, na.rm = T)
  # mcmc.SD.cov[nca,] <- colMeans(SD.all.mcmc.cov, na.rm = T)
  
}
# colnames(results_out) <- c("num_cams", "ID", "EEDE", "REST", "TTE", "MCT", "STE")
colnames(results_out) <- c("num_cams", "ID", "Model", "Est")
results_out$Est <- as.numeric(results_out$Est)
results_out$num_cams <- as.numeric(results_out$num_cams)
results_out <- results_out[!is.na(results_out$Est),]
results_out$Model <- factor(results_out$Model, levels = c("EEDE", "REST", "TTE", "MCT", "STE"))


ggplot(data = results_out[results_out$ID == "mean",], 
       aes(x = num_cams[ID == "mean"], color = Model[ID == "mean"])) +
  geom_ribbon(aes(ymin=results_out$Est[results_out$ID == "mean"] - results_out$Est[results_out$ID == "sd"],
                  ymax=results_out$Est[results_out$ID == "mean"] + results_out$Est[results_out$ID == "sd"], 
                  fill = results_out$Model[results_out$ID == "mean"]),
              linetype = 2,
              alpha=0.6) +
  geom_line(aes(y = Est[ID == "mean"]), size = 1.5) +
  geom_point(aes(y = Est[ID == "mean"]), shape = 4, stroke = 2, color = "black") +
  geom_hline(yintercept=100, linetype="dashed", size=1) +
  labs(x = "Number of Cameras",
       y = paste("Non-Covariate \n Posterior Mean"),
       fill = "Model") +
  guides(color = "none") +
  theme(text = element_text(size = 20),
        legend.title=element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = c(0.87, 0.77)) +
  scale_color_manual(values = fig_colors[2:5]) +
  scale_fill_manual(values = fig_colors[2:5])
# ggsave(paste(fig_dir,"figs/Noncov_subcam_mean.eps", sep = ""), device = cairo_ps)

ggplot(data = results_out[results_out$ID == "mean_cov",], 
       aes(x = num_cams[ID == "mean_cov"], color = Model[ID == "mean_cov"])) +
  geom_ribbon(aes(ymin=results_out$Est[results_out$ID == "mean_cov"] - results_out$Est[results_out$ID == "sd_cov"],
                  ymax=results_out$Est[results_out$ID == "mean_cov"] + results_out$Est[results_out$ID == "sd_cov"], 
                  fill = results_out$Model[results_out$ID == "mean_cov"]),
              linetype = 2,
              alpha=0.6) +
  geom_line(aes(y = Est[ID == "mean_cov"]), size = 1.5) +
  geom_point(aes(y = Est[ID == "mean_cov"]), shape = 4, stroke = 2, color = "black") +
  geom_hline(yintercept=100, linetype="dashed", size=1) +
  labs(x = "Number of Cameras",
       y = paste("Covariate \n Posterior Mean"),
       fill = "Model") +
  guides(color = "none") +
  theme(text = element_text(size = 20),
        legend.title=element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = c(0.87, 0.77)) +
  scale_color_manual(values = fig_colors[1:4]) +
  scale_fill_manual(values = fig_colors[1:4]) 
# ggsave(paste(fig_dir,"figs/Cov_subcam_mean.eps", sep = ""), device = cairo_ps)


ggplot(data = results_out[results_out$ID == "SD",], 
       aes(x = num_cams[ID == "SD"], color = Model[ID == "SD"])) +
  geom_ribbon(aes(ymin=results_out$Est[results_out$ID == "SD"] - results_out$Est[results_out$ID == "SD_sd"],
                  ymax=results_out$Est[results_out$ID == "SD"] + results_out$Est[results_out$ID == "SD_sd"], 
                  fill = results_out$Model[results_out$ID == "SD"]),
              linetype = 2,
              alpha=0.6) +
  geom_line(aes(y = Est[ID == "SD"]), size = 1.5) +
  geom_point(aes(y = Est[ID == "SD"]), shape = 4, stroke = 2, color = "black") +
  labs(x = "Number of Cameras",
       y = paste("Non-Covariate \n Posterior SD"),
       fill = "Model") +
  guides(color = "none") +
  theme(text = element_text(size = 20),
        legend.title=element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = c(0.87, 0.77)) +
  scale_color_manual(values = fig_colors[2:5]) +
  scale_fill_manual(values = fig_colors[2:5])
# ggsave(paste(fig_dir,"figs/Noncov_subcam_sd.eps", sep = ""), device = cairo_ps)


ggplot(data = results_out[results_out$ID == "SD_cov",], 
       aes(x = num_cams[ID == "SD_cov"], color = Model[ID == "SD_cov"])) +
  geom_ribbon(aes(ymin=results_out$Est[results_out$ID == "SD_cov"] - results_out$Est[results_out$ID == "SD_cov_sd"],
                  ymax=results_out$Est[results_out$ID == "SD_cov"] + results_out$Est[results_out$ID == "SD_cov_sd"], 
                  fill = results_out$Model[results_out$ID == "SD_cov"]),
              linetype = 2,
              alpha=0.6) +
  geom_line(aes(y = Est[ID == "SD_cov"]), size = 1.5) +
  geom_point(aes(y = Est[ID == "SD_cov"]), shape = 4, stroke = 2, color = "black") +
  labs(x = "Number of Cameras",
       y = paste("Covariate \n Posterior SD"),
       fill = "Model") +
  guides(color = "none") +
  # scale_x_continuous(expand = c(0, 0)) +
  theme(text = element_text(size = 20),
        legend.title=element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = c(0.87, 0.77)) +
  scale_color_manual(values = fig_colors[1:4]) +
  scale_fill_manual(values = fig_colors[1:4])
# ggsave(paste(fig_dir,"figs/Cov_subcam_sd.eps", sep = ""), device = cairo_ps)




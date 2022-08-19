library(RColorBrewer)
library(GMCM)
library(tidyverse)

# Read saved files for plotting
fig_dir <- "C:/Users/guen.grosklos/Google Drive/Missoula_postdoc/Code/All_models/"
fig_colors <- c("#2ca25f", "#fc8d59", "#67a9cf", "#f768a1", "#bae4b3", "#fed98e")

# Initializations
nind <- 100
leg1<-c("EEDE", "REST", "TTE", "MCT", "STE")
leg.props<-c("EEDE", "REST", "TTE", "MCT")

# Cam sample design (1: random, 2: 80% slow, 3: 80% medium, 4: 80% fast)
cs.all <- 1
cam.dist.labels <- c("random","slow","med","fast")
cam.dist.labels.caps <- c("Random","Slow","Medium","Fast")

# Define landscape variations
# 1: all slow, 2: all medium, 3: all fast, 4: equal slow, medium, fast 5: 80% fast
lv.all <- 4
lv.labels <- c("_slow_lscape_all","_med_lscape_all","_fast_lscape_all","","_fast_lscape")

for (cam.dist.set in cs.all){
  for (lscape_var in lv.all){
    
# Labels for saving figures
means_label <- paste(cam.dist.labels[cam.dist.set],"_cams_means",lv.labels[lscape_var], sep = "")
props_label <- paste(cam.dist.labels[cam.dist.set],"_cams_props",lv.labels[lscape_var], sep = "")
cam.props.label <- paste("Camera Bias: ", cam.dist.labels.caps[cam.dist.set], sep = "")


# # # Load .csv files
D.all.Means.mat.all <- as.matrix(read.csv(paste(fig_dir,"sim_data/all_",means_label,".csv", sep = "")))
D.all.Sds.mat.all <- as.matrix(read.csv(paste(fig_dir,"sim_data/all_",means_label,"_SD.csv", sep = "")))
all.props.Means <- as.matrix(read.csv(paste(fig_dir,"sim_data/",props_label,".csv", sep = "")))
all.props.Sds <- as.matrix(read.csv(paste(fig_dir,"sim_data/",props_label,"_sds.csv", sep = "")))

# Omit MLE results
all.props.Means <- matrix(as.numeric(all.props.Means[c(1,3,5,7,9),2:4]),nrow = 5, ncol = 3)
all.props.Sds <- matrix(as.numeric(all.props.Sds[c(1,3,5,7,9),2:4]),nrow = 5, ncol = 3)
colnames(all.props.Means) <- c("Slow","Medium","fast")
all.props.df <- data.frame(all.props.Means) %>% 
  mutate(Model = c("ABM (Truth)",leg.props))


colnames(all.props.Sds) <- c("Slow","Medium","fast")
all.props.sd.df <- data.frame(all.props.Sds) %>% 
  mutate(Model = c("ABM (Truth)",leg.props))

all.props.sd.df <- all.props.sd.df %>% 
  pivot_longer(cols = !Model,
               names_to = 'Speed',
               values_to = 'sd')

all.props.df <- all.props.df %>% 
  pivot_longer(cols = !Model,
               names_to = 'Speed',
               values_to = 'Proportions') %>% 
  mutate(sd = all.props.sd.df$sd)

all.props.df$Speed <- factor(all.props.df$Speed, levels = c("Slow", "Medium", "Fast"))
all.props.df$Model <- factor(all.props.df$Model, levels = c("ABM (Truth)", "EEDE", "REST", "TTE", "MCT", "STE"))
####################################
# Calculate Summaries
####################################
# Convert all data from wide to long format
D.all.df <- data.frame(D.all.Means.mat.all[,12:21])
D.all.df <- D.all.df %>% 
  pivot_longer(cols = colnames(D.all.df),
               names_to='Model',
               values_to='Est') %>% 
  mutate(Covariate = NA)
D.all.df$Covariate[grepl(".2", D.all.df$Model, fixed = TRUE)] <- "Non-Covariate"
D.all.df$Covariate[grepl(".3", D.all.df$Model, fixed = TRUE)] <- "Covariate"
D.all.df$Model <- gsub(x = D.all.df$Model, pattern = "(.2)", replacement = "")
D.all.df$Model <- gsub(x = D.all.df$Model, pattern = "(.3)", replacement = "")
D.all.df$Model[grepl(".3", D.all.df$Model, fixed = TRUE)] <- "Yes"
D.all.df$Est[is.na(D.all.df$Est)] <- 100
D.all.df$Model <- factor(D.all.df$Model, levels = c("EEDE", "REST", "TTE", "MCT", "STE"))
D.all.df$Covariate <- factor(D.all.df$Covariate, levels = c("Non-Covariate", "Covariate"))
  
# Convert all data from wide to long format
SD.all.df <- data.frame(D.all.Sds.mat.all[,12:21])
SD.all.df <- SD.all.df %>% 
  pivot_longer(cols = colnames(SD.all.df),
               names_to='Model',
               values_to='Est') %>% 
  mutate(Covariate = NA)
SD.all.df$Covariate[grepl(".2", SD.all.df$Model, fixed = TRUE)] <- "Non-Covariate"
SD.all.df$Covariate[grepl(".3", SD.all.df$Model, fixed = TRUE)] <- "Covariate"
SD.all.df$Model <- gsub(x = SD.all.df$Model, pattern = "(.2)", replacement = "")
SD.all.df$Model <- gsub(x = SD.all.df$Model, pattern = "(.3)", replacement = "")
SD.all.df$Model[grepl(".3", SD.all.df$Model, fixed = TRUE)] <- "Yes"
SD.all.df$Est[is.na(SD.all.df$Est)] <- 0
SD.all.df$Model <- factor(SD.all.df$Model, levels = c("EEDE", "REST", "TTE", "MCT", "STE"))
SD.all.df$Covariate <- factor(SD.all.df$Covariate, levels = c("Non-Covariate", "Covariate"))

if (any(lv.all == 1:3)) {
  ggplot(D.all.df, aes(x=Model, y=Est, fill=Covariate)) +
    geom_boxplot(lwd = .1, fatten = .1) +
    geom_boxplot(data=subset(D.all.df, D.all.df$Model %in% c("REST", "TTE", "MCT")), colour = "black") +
    geom_boxplot(data=subset(D.all.df, D.all.df$Model == "EEDE"), colour=c("white","black")) +
    geom_boxplot(data=subset(D.all.df, D.all.df$Model == "STE"), colour=c("black","white")) +
    geom_hline(yintercept=100, linetype="dashed", size=1) +
    labs(x = "Model",
         y = paste(cam.props.label,"\n Posterior Mean")) +
    theme(text = element_text(size = 20),
          legend.title=element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.position = c(0.17, 0.84),
          legend.background = element_blank(),
          legend.spacing.y = unit(0, "mm"), 
          legend.box.background = element_rect(colour = "black")) 
  # ggsave(paste(fig_dir,"figs/",means_label,"_box.eps", sep = ""), device = cairo_ps)
  
  
} else{
  ggplot(D.all.df, aes(x=Model, y=Est, fill=Covariate)) +
    geom_boxplot(lwd = .1, fatten = .1) +
    geom_boxplot(data=subset(D.all.df, D.all.df$Model %in% c("REST", "TTE", "MCT")), colour = "black") +
    geom_boxplot(data=subset(D.all.df, D.all.df$Model == "EEDE"), colour=c("white","black")) +
    geom_boxplot(data=subset(D.all.df, D.all.df$Model == "STE"), colour=c("black","white")) +
    geom_hline(yintercept=100, linetype="dashed", size=1) +
    labs(x = "Model",
         y = paste(cam.props.label,"\n Posterior Mean")) +
    theme(text = element_text(size = 20),
          legend.title=element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.position = c(0.17, 0.84),
          legend.background = element_blank(),
          legend.spacing.y = unit(0, "mm"), 
          legend.box.background = element_rect(colour = "black"))
    # ggsave(paste(fig_dir,"figs/",means_label,"_box.eps", sep = ""), device = cairo_ps)
  
  ggplot(all.props.df, aes(x=Speed, y = Proportions, fill=Model)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=Proportions-sd, ymax=Proportions+sd), width=.2,
                  position=position_dodge(.9), size = .7) +
    scale_y_continuous(limits=c(0, max(all.props.df$Proportions+all.props.df$sd) +.01), expand = c(0, 0)) +
    labs(x = "Speed",
         y = "Proportion Occupied") +
    scale_fill_manual(values = c("black",fig_colors[1:4])) +
    theme(text = element_text(size = 20),
          legend.title=element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.position = c(0.87, 0.7),
          legend.background = element_blank(),
          legend.spacing.y = unit(0, "mm"), 
          legend.box.background = element_rect(colour = "black")) 
  # ggsave(paste(fig_dir,"figs/",props_label,".eps", sep = ""), device = cairo_ps)
  
  ggplot(SD.all.df, aes(x=Model, y=Est, fill=Covariate)) +
    geom_boxplot() +
    labs(x = "Model",
         y = paste(cam.props.label,"\n Posterior SD")) +
    scale_y_continuous(limits=c(0, max(SD.all.df$Est) +2), expand = c(0, 0)) +
    theme(text = element_text(size = 20),
          legend.title=element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.position = c(0.17, 0.84),
          legend.background = element_blank(),
          legend.spacing.y = unit(0, "mm"), 
          legend.box.background = element_rect(colour = "black")) 
  # ggsave(paste(fig_dir,"figs/SD_",means_label,"_box.eps", sep = ""), device = cairo_ps)
  

}
  }
}

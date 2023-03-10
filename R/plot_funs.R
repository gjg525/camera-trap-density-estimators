

################################
## Plot ABM simulations
################################
plot_ABM <- function() {
  # Sample area border
  b.df <- data.frame(c(bounds,dx))

  # camera cell locations
  cam.df <- tibble(cam = 1:ncam,
                   xmin = dx*(cam.samps - 1) %% q^0.5,
                   ymin = dx*ceiling(cam.samps / q^0.5) - dx,
                   xmax = xmin + dx,
                   ymax = ymin + dy,
                   a = (ymax - ymin) * (xmax - xmin)
  )
  # Triangle viewshed cameras
  dt.triangle <- data.frame(group = rep(1:ncam, each = 3),
                            polygon.x = tri_cam_samps$x,
                            polygon.y = tri_cam_samps$y)


  g <- ggplot(lscape_speeds) +
    # geom_tile(aes(X-.5*dx, Y-.5*dy, fill = Speed)) +
    scale_fill_manual(values= c("grey20", "grey50", "grey80")) +
    geom_path(data = animalxy.all,
              aes(x, y, col = as.factor(Animal_ID))) +
    geom_rect(data = b.df,
              aes(xmin = bounds[1], xmax = bounds[2], ymin = bounds[1], ymax = bounds[2]),
              color="black",
              fill= NA) +
    # geom_rect(data = cam.df,
    #           aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax),
    #           color="black",
    #           fill= NA) +
    geom_polygon(data = dt.triangle,
                 aes(x=polygon.x, y=polygon.y, group=group),
                 color = "black",
                 fill = NA) +
    labs(x = "X", y = "Y", fill="Landscape Type") +
    guides(colour = "none",
           alpha = "none") +
    theme(panel.background = element_blank())
  g + coord_fixed()
}

plot_space_use <- function() {
  # u.abm.all <- raster(, nrows = q^0.5, ncols = q^0.5, xmn = 0,
  #                     xmx = q^0.5, ymn = 0, ymx = q^0.5, crs = NA)
  # u.abm.all[] <- 0
  # u.abm.all <- stack(replicate(t.steps,u.abm.all))

  u.abm.all = matrix(0, nrow = q^0.5, ncol = q^0.5)

  # animalxy.discrete <- animalxy.all |>
  #   filter(t %in% 1:t.steps)
  for(xx in 1:nrow(animalxy.all)) {
    x.round <- ceiling(animalxy.all$x[xx]*(q^0.5/max(bounds)))
    y.round <- ceiling(animalxy.all$y[xx]*(q^0.5/max(bounds)))

    # Transpose for converting matrix to raster
    u.abm.all[x.round,y.round] <- u.abm.all[x.round,y.round]+1
  }

  u_plot <- levelplot(u.abm.all,
                      layers = 1,
                      main = list(expression("Animal space-use"), cex = 1),
                      cuts = 254,
                      margin = FALSE,
                      scales = list(draw = FALSE),
                      xlab = "x",
                      ylab = "y",
                      col.regions = colorRampPalette(rev(brewer.pal(11, "Spectral")), bias = 1))

  grid.arrange(u_plot, ncol = 1, nrow = 1)
}


################################
## Plot data density plots
################################
# # Density plot for count data
plot_count_data <- function(fill = "speed") {
  ggplot(count_data, aes(x=count, fill = .data[[fill]])) +
    geom_density(position = "identity", alpha = 0.4, adjust = 3) +
    labs(x = "Count", y = "Frequency", fill="Landscape Type") +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=22),
          legend.text=element_text(size=20),
          legend.title = element_text(size=20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1))

  # ggplot(count_data, aes(x=count, fill = .data[[fill]])) +
  #   geom_histogram(position = "identity", alpha = 0.4) +
  #   labs(x = "Count", y = "Frequency", fill="Landscape Type") +
  #   theme(axis.text=element_text(size=20),
  #         axis.title=element_text(size=22),
  #         legend.text=element_text(size=20),
  #         legend.title = element_text(size=20),
  #         panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(),
  #         panel.background = element_blank(),
  #         panel.border = element_rect(colour = "black", fill=NA, size=1))

}

# # Density plot for encounter data
plot_encounter_data <- function(fill = "speed") {
  ggplot(encounter_data, aes(x=encounter, fill = .data[[fill]])) +
    geom_density(position = "identity", alpha = 0.4, adjust = 2) +
    labs(x = "Number of Encounters", y = "Frequency", fill="Landscape Type") +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=22),
          legend.text=element_text(size=20),
          legend.title = element_text(size=20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1))

  # ggplot(encounter_data, aes(x=encounter, fill = .data[[fill]])) +
  #   geom_histogram(position = "identity", alpha = 0.4) +
  #   labs(x = "Count", y = "Frequency", fill="Landscape Type") +
  #   theme(axis.text=element_text(size=20),
  #         axis.title=element_text(size=22),
  #         legend.text=element_text(size=20),
  #         legend.title = element_text(size=20),
  #         panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(),
  #         panel.background = element_blank(),
  #         panel.border = element_rect(colour = "black", fill=NA, size=1))

}

# # Density plot for staying time data
plot_staytime_data <- function(fill = "speed") {
  ggplot(stay_time_raw, aes(x=t_stay, fill = .data[[fill]])) +
    geom_density(position = "identity", alpha = 0.4, adjust = 5) +
    labs(x = "Staying Time", y = "Frequency", fill="Landscape Type") +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=22),
          legend.text=element_text(size=20),
          legend.title = element_text(size=20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1))
}

# # Density plot for TTE data
plot_TTE_data <- function(fill = "speed") {
  ggplot(TTE_data_raw, aes(x=TTE, fill = .data[[fill]])) +
    geom_density(position = "identity", alpha = 0.4, adjust = 2) +
    labs(x = "Time to Encounter", y = "Frequency", fill = "Landscape Type") +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=22),
          legend.text=element_text(size=20),
          legend.title = element_text(size=20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1))
}

################################
## Plot single-run results
################################

plot_onerun_results <- function() {
  ggplot(D.all, aes(x = Model)) +
    geom_point(aes(y = Est)) +
    geom_hline(yintercept=nind, linetype="dashed") +
    geom_errorbar(aes(ymin = Est - SD, ymax = Est + SD), width = 0.2)
}


################################
## Plot multi-run results
################################

plot_multirun_means <- function() {
  ggplot(D.all, aes(x = Model, y = Est, fill = Covariate)) +
    geom_boxplot(lwd = .1, fatten = .1) +
    geom_boxplot(data=subset(D.all, D.all$Model %in% c("REST", "TTE", "MCT")), colour = "black") +
    geom_boxplot(data=subset(D.all, D.all$Model == "TDST"), colour=c("black","white")) +
    geom_boxplot(data=subset(D.all, D.all$Model == "STE"), colour=c("white","black")) +
    geom_hline(yintercept=100, linetype="dashed", size=1) +
    labs(x = "Model",
         # y = paste(cam.props.label[lscape_var+1], " Landscape \n Abundance Estimates")) +
         y = "Abundance Estimates") +
    scale_fill_manual(values=c('grey40','Grey')) +
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
}

plot_grouped_multirun_means <- function(Unused_cov = "None", 
                                        Filter_model = "None", 
                                        Cov = "Covariate") {

  # Filter with covariate labels not used
  D.all %>% 
    dplyr::filter(Model != Filter_model) %>% 
    dplyr::filter(Covariate != Unused_cov) %>%
    ggplot(aes(x = Model, y = Est, fill = Run)) +
    geom_boxplot(lwd = 0.5, fatten = .5, outlier.size = 1) +
    geom_hline(yintercept=nind, linetype="dashed", size=1) +
    labs(x = "Model",
         y = paste0(Cov,"\n", "Abundance Estimates")) +
    scale_fill_manual(values=if(sim_num == 4){ 
      c('grey40', 'grey60','grey80')} 
      else {c('grey20', 'grey40', 'grey60','grey80')}
    ) +
    theme(text = element_text(size = 20),
          legend.title=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.position = "none",
          legend.background = element_blank(),
          legend.spacing.y = unit(0, "mm"),
          legend.box.background = element_rect(colour = "black"))
}

plot_multirun_sds <- function() {
  ggplot(D.all, aes(x = Model, y = SD, fill = Covariate)) +
    geom_boxplot(lwd = .1, fatten = .1) +
    geom_boxplot(data=subset(D.all, D.all$Model %in% c("REST", "TTE", "MCT")), colour = "black") +
    geom_boxplot(data=subset(D.all, D.all$Model == "TDST"), colour=c("black","white")) +
    geom_boxplot(data=subset(D.all, D.all$Model == "STE"), colour=c("white","black")) +
    labs(x = "Model",
         # y = paste(cam.props.label[lscape_var+1], " Landscape \n Abundance Estimates")) +
         y = "SD") +
    scale_fill_manual(values=c('grey40','Grey')) +
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
}

plot_multirun_CV <- function() {
  ggplot(D.all, aes(x = Model, y = SD/Est, fill = Covariate)) +
    geom_boxplot(lwd = .1, fatten = .1) +
    geom_boxplot(data=subset(D.all, D.all$Model %in% c("REST", "TTE", "MCT")), colour = "black") +
    geom_boxplot(data=subset(D.all, D.all$Model == "TDST"), colour=c("black","white")) +
    geom_boxplot(data=subset(D.all, D.all$Model == "STE"), colour=c("white","black")) +
    labs(x = "Model",
         y = "CV") +
    scale_fill_manual(values=c('grey40','Grey')) +
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
}

plot_grouped_multirun_CV <- function(Unused_cov = "None", 
                                     Filter_model = "None",
                                     Title) {
  D.all %>% 
    dplyr::filter(Model != Filter_model) %>% 
    dplyr::filter(Covariate != Unused_cov) %>%
    ggplot(aes(x = Model, y = SD/Est, fill = Run)) +
    geom_boxplot(lwd = 0.5, fatten = .5, outlier.size = 1) +
    labs(x = "Model",
         y = "CV") +
    guides(fill = guide_legend(title = Title)) +
    scale_fill_manual(values=if(sim_num == 4){ 
      c('grey40', 'grey60','grey80')} 
      else {c('grey20', 'grey40', 'grey60','grey80')}
    ) +
    theme(text = element_text(size = 20),
          legend.title=element_text(size=18),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.position = if(Unused_cov == "None"){c(0.856, 0.8)}
          else{c(0.142, 0.762)},
          legend.background = element_blank(),
          legend.spacing.y = unit(0, "mm"))
}

plot_multirun_hist <- function() {
  mean_data <- D.all |>
    group_by(Model) |>
    summarise(Means = mean(Est),
              SDs = sd(Est))

  ggplot(D.all, aes(x = Est, fill = Model)) +
    geom_density(position = "identity", alpha = 0.3, adjust = 1) +
    # geom_vline(data = mean_data, aes(xintercept = Means, color = Model)) +
    labs(x = "Estimate", y = "Density", fill="Model") +
    theme(text = element_text(size = 20),
          legend.title=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.background = element_blank(),
          legend.spacing.y = unit(0, "mm"),
          legend.box.background = element_rect(colour = "black"))
}

plot_model_proportions <- function() {
  # number of counts across whole landscape for each covariate type
  slow.counts <- animal_data %>% filter(lscape_type == "Slow") %>% pull(t_spent)
  med.counts <- animal_data %>% filter(lscape_type == "Medium") %>% pull(t_spent)
  fast.counts <- animal_data %>% filter(lscape_type == "Fast") %>% pull(t_spent)
  
  Prop_all <- D.all %>% 
    dplyr::filter(Covariate == "Covariate") %>% 
    unnest_wider(Prop_speeds, names_sep="_") %>% 
    rename(Slow = Prop_speeds_1,
           Medium = Prop_speeds_2,
           Fast = Prop_speeds_3) %>% 
    select(Model, Slow, Medium, Fast) %>% 
    pivot_longer(!Model, names_to = "Speed", values_to = "Proportions") %>% 
    group_by(Model, Speed) %>% 
    summarise(Means = mean(Proportions),
              SDs = sd(Proportions),
              .groups = 'drop') %>% 
    add_row(Model = "ABM",
            Speed = c("Slow", "Medium", "Fast"),
            Means = c(slow.counts, med.counts, fast.counts))
  
  
  ggplot(Prop_all, aes(x=Speed, y = Means, fill = Model)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=Means-SDs, ymax=Means+SDs), width=.2,
                  position=position_dodge(.9), size = .7) +
    # scale_y_continuous(limits=c(0, max(Prop_all$Means) +.01), expand = c(0, 0)) +
    labs(x = "Landscape Type",
         y = "Relative Distributions") +
    scale_fill_manual(values = c("grey40",fig_colors[1:4])) +
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
}

plot_ABM_stay_proportions <- function() {
  # number of counts across whole landscape for each covariate type
  slow.counts <- animal_data %>% filter(lscape_type == "Slow") %>% pull(t_spent)
  med.counts <- animal_data %>% filter(lscape_type == "Medium") %>% pull(t_spent)
  fast.counts <- animal_data %>% filter(lscape_type == "Fast") %>% pull(t_spent)
  
  Prop_all <- stay_time_raw %>% 
    dplyr::group_by(speed) |> 
    dplyr::summarise(Model = "Stay Time", 
                     Means = mean(t_stay),
                     SDs = sd(t_stay),
                     .groups = 'drop') %>% 
    dplyr::add_row(Model = "ABM",
                   speed = c("Slow", "Medium", "Fast"),
                   Means = c(slow.counts, med.counts, fast.counts)) |> 
    dplyr::group_by(Model) |> 
    dplyr::mutate(Means = Means / sum(Means))
    
  
  
  ggplot(Prop_all, aes(x=speed, y = Means/sum(Means), fill = Model)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    # geom_errorbar(aes(ymin=Means-SDs, ymax=Means+SDs), width=.2,
    #               position=position_dodge(.9), size = .7) +
    labs(x = "Landscape Type",
         y = "Relative Distributions") +
    scale_fill_manual(values = c("grey40",fig_colors[1:4])) +
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
}

# Old code. Use for color schemes
######################

# library(RColorBrewer)
# library(GMCM)
# library(tidyverse)
# 
# # Read saved files for plotting
# fig_dir <- "C:/Users/guen.grosklos/Google Drive/Missoula_postdoc/Code/All_models/"
# fig_colors <- c("#2ca25f", "#fc8d59", "#67a9cf", "#f768a1", "#bae4b3", "#fed98e")
# 
# # Initializations
# nind <- 100
# leg1<-c("EEDE", "REST", "TTE", "MCT", "STE")
# leg.props<-c("EEDE", "REST", "TTE", "MCT")
# 
# # Cam sample design (1: random, 2: 80% slow, 3: 80% medium, 4: 80% fast)
# cs.all <- 1
# cam.dist.labels <- c("random","slow","med","fast")
# cam.dist.labels.caps <- c("Random","Slow","Medium","Fast")
# 
# # Define landscape variations
# # 1: all slow, 2: all medium, 3: all fast, 4: equal slow, medium, fast 5: 80% fast
# lv.all <- 1:3
# lv.labels <- c("_slow_lscape_all","_med_lscape_all","_fast_lscape_all","","_fast_lscape")
# 
# for (cam.dist.set in cs.all){
#   for (lscape_var in lv.all){
#     
#     # Labels for saving figures
#     means_label <- paste(cam.dist.labels[cam.dist.set],"_cams_means",lv.labels[lscape_var], sep = "")
#     props_label <- paste(cam.dist.labels[cam.dist.set],"_cams_props",lv.labels[lscape_var], sep = "")
#     cam.props.label <- paste("Camera Bias: ", cam.dist.labels.caps[cam.dist.set], sep = "")
#     
#     
#     # # # Load .csv files
#     D.all.Means.mat.all <- as.matrix(read.csv(paste(fig_dir,"sim_data/all_",means_label,".csv", sep = "")))
#     D.all.Sds.mat.all <- as.matrix(read.csv(paste(fig_dir,"sim_data/all_",means_label,"_SD.csv", sep = "")))
#     all.props.Means <- as.matrix(read.csv(paste(fig_dir,"sim_data/",props_label,".csv", sep = "")))
#     all.props.Sds <- as.matrix(read.csv(paste(fig_dir,"sim_data/",props_label,"_sds.csv", sep = "")))
#     
#     # Omit MLE results
#     all.props.Means <- matrix(as.numeric(all.props.Means[c(1,3,5,7,9),2:4]),nrow = 5, ncol = 3)
#     all.props.Sds <- matrix(as.numeric(all.props.Sds[c(1,3,5,7,9),2:4]),nrow = 5, ncol = 3)
#     colnames(all.props.Means) <- c("Slow","Medium","Fast")
#     all.props.df <- data.frame(all.props.Means) %>% 
#       mutate(Model = c("ABM (Truth)",leg.props))
#     
#     
#     colnames(all.props.Sds) <- c("Slow","Medium","Fast")
#     all.props.sd.df <- data.frame(all.props.Sds) %>% 
#       mutate(Model = c("ABM (Truth)",leg.props))
#     
#     all.props.sd.df <- all.props.sd.df %>% 
#       pivot_longer(cols = !Model,
#                    names_to = 'Speed',
#                    values_to = 'sd')
#     
#     all.props.df <- all.props.df %>% 
#       pivot_longer(cols = !Model,
#                    names_to = 'Speed',
#                    values_to = 'Proportions') %>% 
#       mutate(sd = all.props.sd.df$sd)
#     
#     all.props.df$Speed <- factor(all.props.df$Speed, levels = c("Slow", "Medium", "Fast"))
#     all.props.df$Model <- factor(all.props.df$Model, levels = c("ABM (Truth)", "EEDE", "REST", "TTE", "MCT", "STE"))
#     ####################################
#     # Calculate Summaries
#     ####################################
#     # Convert all data from wide to long format
#     D.all.df <- data.frame(D.all.Means.mat.all[,12:21])
#     D.all.df <- D.all.df %>% 
#       pivot_longer(cols = colnames(D.all.df),
#                    names_to='Model',
#                    values_to='Est') %>% 
#       mutate(Covariate = NA)
#     D.all.df$Covariate[grepl(".2", D.all.df$Model, fixed = TRUE)] <- "Non-Covariate"
#     D.all.df$Covariate[grepl(".3", D.all.df$Model, fixed = TRUE)] <- "Covariate"
#     D.all.df$Model <- gsub(x = D.all.df$Model, pattern = "(.2)", replacement = "")
#     D.all.df$Model <- gsub(x = D.all.df$Model, pattern = "(.3)", replacement = "")
#     D.all.df$Model[grepl(".3", D.all.df$Model, fixed = TRUE)] <- "Yes"
#     D.all.df$Est[is.na(D.all.df$Est)] <- 100
#     D.all.df$Model <- factor(D.all.df$Model, levels = c("EEDE", "REST", "TTE", "MCT", "STE"))
#     D.all.df$Covariate <- factor(D.all.df$Covariate, levels = c("Non-Covariate", "Covariate"))
#     
#     # Convert all data from wide to long format
#     SD.all.df <- data.frame(D.all.Sds.mat.all[,12:21])
#     SD.all.df <- SD.all.df %>% 
#       pivot_longer(cols = colnames(SD.all.df),
#                    names_to='Model',
#                    values_to='Est') %>% 
#       mutate(Covariate = NA)
#     SD.all.df$Covariate[grepl(".2", SD.all.df$Model, fixed = TRUE)] <- "Non-Covariate"
#     SD.all.df$Covariate[grepl(".3", SD.all.df$Model, fixed = TRUE)] <- "Covariate"
#     SD.all.df$Model <- gsub(x = SD.all.df$Model, pattern = "(.2)", replacement = "")
#     SD.all.df$Model <- gsub(x = SD.all.df$Model, pattern = "(.3)", replacement = "")
#     SD.all.df$Model[grepl(".3", SD.all.df$Model, fixed = TRUE)] <- "Yes"
#     SD.all.df$Est[is.na(SD.all.df$Est)] <- 0
#     SD.all.df$Model <- factor(SD.all.df$Model, levels = c("EEDE", "REST", "TTE", "MCT", "STE"))
#     SD.all.df$Covariate <- factor(SD.all.df$Covariate, levels = c("Non-Covariate", "Covariate"))
#     
#     if (any(lv.all %in% 1:3)) {
#       ggplot(D.all.df, aes(x=Model, y=Est, fill=Covariate)) +
#         geom_boxplot(lwd = .1, fatten = .1) +
#         geom_boxplot(data=subset(D.all.df, D.all.df$Model %in% c("REST", "TTE", "MCT")), colour = "black") +
#         geom_boxplot(data=subset(D.all.df, D.all.df$Model == "EEDE"), colour=c("white","black")) +
#         geom_boxplot(data=subset(D.all.df, D.all.df$Model == "STE"), colour=c("black","white")) +
#         geom_hline(yintercept=100, linetype="dashed", size=1) +
#         labs(x = "Model",
#              y = paste(cam.dist.labels.caps[lscape_var+1], " Landscape \n Abundance Estimates")) +
#         scale_fill_manual(values=c('grey40','Grey')) +
#         theme(text = element_text(size = 20),
#               legend.title=element_blank(), 
#               panel.grid.major = element_blank(), 
#               panel.grid.minor = element_blank(),
#               panel.background = element_blank(), 
#               axis.line = element_line(colour = "black"),
#               panel.border = element_rect(colour = "black", fill=NA, size=1),
#               legend.position = c(0.17, 0.84),
#               legend.background = element_blank(),
#               legend.spacing.y = unit(0, "mm"), 
#               legend.box.background = element_rect(colour = "black")) 
#       # ggsave(paste(fig_dir,"figs/",means_label,"_box.eps", sep = ""), device = cairo_ps)
#       
#       
#     } else{
#       ggplot(D.all.df, aes(x=Model, y=Est, fill=Covariate)) +
#         geom_boxplot(lwd = .1, fatten = .1) +
#         geom_boxplot(data=subset(D.all.df, D.all.df$Model %in% c("REST", "TTE", "MCT")), colour = "black") +
#         geom_boxplot(data=subset(D.all.df, D.all.df$Model == "EEDE"), colour=c("white","black")) +
#         geom_boxplot(data=subset(D.all.df, D.all.df$Model == "STE"), colour=c("black","white")) +
#         geom_hline(yintercept=100, linetype="dashed", size=1) +
#         labs(x = "Model",
#              y = paste(cam.props.label,"\n Abundance Estimates")) +
#         scale_fill_manual(values=c('grey40','Grey')) +
#         theme(text = element_text(size = 20),
#               legend.title=element_blank(), 
#               panel.grid.major = element_blank(), 
#               panel.grid.minor = element_blank(),
#               panel.background = element_blank(), 
#               axis.line = element_line(colour = "black"),
#               panel.border = element_rect(colour = "black", fill=NA, size=1),
#               legend.position = c(0.17, 0.84),
#               legend.background = element_blank(),
#               legend.spacing.y = unit(0, "mm"), 
#               legend.box.background = element_rect(colour = "black"))
#       # ggsave(paste(fig_dir,"figs/",means_label,"_box.eps", sep = ""), device = cairo_ps)
#       
#       ggplot(all.props.df, aes(x=Speed, y = Proportions, fill = Model)) +
#         geom_bar(stat="identity", color="black", position=position_dodge()) +
#         geom_errorbar(aes(ymin=Proportions-sd, ymax=Proportions+sd), width=.2,
#                       position=position_dodge(.9), size = .7) +
#         scale_y_continuous(limits=c(0, max(all.props.df$Proportions+all.props.df$sd) +.01), expand = c(0, 0)) +
#         labs(x = "Landscape Type",
#              y = "Relative Distributions") +
#         scale_fill_manual(values = c("grey40",fig_colors[1:4])) +
#         theme(text = element_text(size = 20),
#               legend.title=element_blank(), 
#               panel.grid.major = element_blank(), 
#               panel.grid.minor = element_blank(),
#               panel.background = element_blank(), 
#               axis.line = element_line(colour = "black"),
#               panel.border = element_rect(colour = "black", fill=NA, size=1),
#               legend.position = c(0.87, 0.7),
#               legend.background = element_blank(),
#               legend.spacing.y = unit(0, "mm"), 
#               legend.box.background = element_rect(colour = "black")) 
#       # ggsave(paste(fig_dir,"figs/",props_label,".eps", sep = ""), device = cairo_ps)
#       
#       ggplot(SD.all.df, aes(x=Model, y=Est, fill=Covariate)) +
#         geom_boxplot() +
#         labs(x = "Model",
#              y = paste(cam.props.label,"\n SD Results")) +
#         scale_y_continuous(limits=c(0, max(SD.all.df$Est) +2), expand = c(0, 0)) +
#         scale_fill_manual(values=c('grey40','Grey')) +
#         theme(text = element_text(size = 20),
#               legend.title=element_blank(), 
#               panel.grid.major = element_blank(), 
#               panel.grid.minor = element_blank(),
#               panel.background = element_blank(), 
#               axis.line = element_line(colour = "black"),
#               panel.border = element_rect(colour = "black", fill=NA, size=1),
#               legend.position = c(0.17, 0.84),
#               legend.background = element_blank(),
#               legend.spacing.y = unit(0, "mm"), 
#               legend.box.background = element_rect(colour = "black")) 
#       # ggsave(paste(fig_dir,"figs/SD_",means_label,"_box.eps", sep = ""), device = cairo_ps)
#       
#       
#     }
#   }
# }


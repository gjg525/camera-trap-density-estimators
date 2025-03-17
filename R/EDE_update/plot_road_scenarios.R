library(dplyr)
library(ggplot2)
tot_animals <- 25

fig_colors <- c("#2ca25f", "#fc8d59", "#67a9cf", "#f768a1", "#bae4b3", "#fed98e")

file_names <- c(
  # "random",
  "slow",
  "fast"
)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

source("./R/EDE_update/utils.R")
source("./R/EDE_update/plot_funs.R")


D_all_road <- tibble::tibble()
for (ii in 1:length(file_names)) {
  print(ii)
  all_results <- loadRData(paste0("Sim_results/results_road_",
                                  file_names[ii],
                                  "_cam.RData"))

  D_all_road <- dplyr::bind_rows(
    D_all_road,
    all_results[[5]] %>%
      dplyr::mutate(
        SampDesign = paste0(file_names[ii], "_cam")
      )
  )
}

# save(D_all_road, file = "G:/My Drive/Missoula_postdoc/TDST_Model/Data/D_all_road.RData")

D_all_road %>% 
  dplyr::filter(Est < 250) %>%
  ggplot2::ggplot(ggplot2::aes(x = SampDesign, y = Est, fill = Model)) +
  ggplot2::geom_boxplot(lwd = 0.5, fatten = .5, outlier.shape = NA) +
  ggplot2::geom_hline(yintercept=tot_animals, linetype="dashed", size = 0.7) +
  ggplot2::labs(x = "Camera Sample Design",
                y = "Posterior Means") +
  ggplot2::scale_fill_manual(values= fig_colors[1:5]) +
  ggplot2::scale_color_manual(values = c('grey0', 'grey40', 'grey60')) +
  ggplot2::theme(text = ggplot2::element_text(size = 16),
                 panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(),
                 axis.line = ggplot2::element_line(colour = "black"),
                 panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
                 legend.position = "none",
                 legend.background = ggplot2::element_blank(),
                 legend.spacing.y = ggplot2::unit(0, "mm"),
                 legend.box.background = ggplot2::element_rect(colour = "black"))

# prerun
# source("02_0_FFQ_preProcessed.R")
# source("02_diversity_food_frequency.R")
# source("03_bmi_and_adiposity.R")

rm(list=ls()) 

library(ggplot2)
library(patchwork)

# Figure 3
species_scores <- read.delim("04.biometrics/species_scores.tsv") 
rownames(species_scores) <- species_scores$food

file_list = list.files("03.dietary_assesment/", 
                       pattern = "^Fig3.*\\.rds$", 
                       full.names = TRUE)
panel_list <- panel_names <- gsub("^03\\.dietary_assesment/Fig3|\\.rds$", "", file_list)

plot_ <- lapply(file_list, readRDS)
names(plot_) <- panel_list

col_A <- (plot_[["A"]] / plot_spacer()) + plot_layout(heights = c(0.6, 0.4))
col_B <- (plot_[["B"]] / plot_spacer()) +   plot_layout(heights = c(0.85, 0.15))
col_C <- plot_[["C"]]

# Combine columns horizontally
row1 <- (col_A | col_B | col_C) + 
  plot_layout(widths = c(0.3, 0.45, 0.25), guides = "collect") & 
  theme(legend.position = "bottom", legend.title = element_blank()) 
row2 <- plot_[["D"]]

figure3 <- (row1 / row2) + 
  plot_layout(heights = c(0.7, 0.3)) + 
  plot_annotation(tag_levels = list(c("A", "B", "C", "D"))) & 
  theme(
    legend.position = "bottom",
    legend.text = element_text(size=9),
    legend.key.size = unit(4, "mm"),
    axis.text = element_text(size=10),
    plot.tag = element_text(size = 20, face = "bold"),
    plot.margin = margin(0, 1, 0, 0, "mm"),
    panel.spacing = unit(0, "mm"),
  )
figure3

figout="03.dietary_assesment/Figure3_collated"
for(dev in c("svg","png","pdf")){
  outpath=paste0(figout,".",dev)
  ggsave(plot = figure3, filename = outpath, device = dev, 
         width = 210, height = 210, units = "mm", dpi = 300,
         scale = 1.2)
}

# Figure 4
file_list = list.files("04.biometrics/", 
                       pattern = "^Fig4.*\\.rds$", 
                       full.names = TRUE)
panel_list <- panel_names <- gsub("^04\\.biometrics/Fig4|\\.rds$", "", file_list)

plot_ <- lapply(file_list, readRDS)
names(plot_) <- panel_list


row1 <- (plot_[["A1"]] + plot_[["A2"]]) + 
  plot_layout(widths = c(1.2, 2), guides = "collect") + 
  plot_annotation(
    title = "Underweight and Overweight Status by Cohort",
    theme = theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
                  axis.text = element_text(size=9),
                  title=element_text(size=10),
                  legend.position = "top",
                  legend.text = element_text(size=9),
                  legend.key.size = unit(3, "mm"),
                  panel.spacing = unit(1, "mm"),
                  plot.margin = margin(-15,-7,5,-35, "pt"))
  )

col_B <- (plot_[["B"]] / plot_spacer()) + plot_layout(heights = c(1, 0))
col_C <- (plot_[["C"]] / plot_spacer()) + plot_layout(heights = c(0.85,0.15))
row2 <- (col_B | col_C) + 
  plot_layout(widths = c(2.2, 0.8), guides = "collect") & 
  theme(legend.position = "bottom")

col_D <- (plot_[["D"]] / plot_spacer()) + plot_layout(heights = c(1,0))
col_E <- (plot_[["E"]] / plot_spacer()) + plot_layout(heights = c(0.85,0.15))
row3 <- (col_D | col_E) + 
  plot_layout(widths = c(2.2, 0.8), guides = "collect") & 
  theme(legend.position = "bottom")

# Combine rows and set vertical height ratios
figure4 <- (wrap_elements(row1) / row2 / row3) + 
  plot_layout(heights = c(0.9,1,1)) + 
  plot_annotation(tag_levels = list(c("A", "B", "C", "D", "E"))) & 
  theme(
    legend.position = "bottom",
    legend.text = element_text(size=9),
    legend.key.size = unit(4, "mm"),
    axis.text = element_text(size=10),
    plot.tag = element_text(size = 20, face = "bold"),
    plot.margin = margin(0, 1, 0, 0, "mm"),
    panel.spacing = unit(0, "mm")
  )

figure4

figout="04.biometrics/Figure4_collated"
for(dev in c("svg","png","pdf")){
  outpath=paste0(figout,".",dev)
  ggsave(plot = figure4, filename = outpath, device = dev, 
         width = 210, height = 270, units = "mm", dpi = 300,
         scale = 0.9)
}

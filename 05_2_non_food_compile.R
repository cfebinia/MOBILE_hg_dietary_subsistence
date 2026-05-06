rm(list=ls()) 

library(ggplot2)
library(patchwork)

# prerun
# source("sugar.R")
# source("smoking.R")
# source("medicine.R")

# format

d_width=0.7

# 1. Load panels and set individual guide configurations
panelA <- readRDS("fig_out/05.non_food/Panel_A.rds") + 
  scale_y_continuous(limits = c(-5,450))+
  guides(fill = guide_legend(ncol = 1), col = guide_legend(ncol = 1), shape = guide_legend(ncol = 1))
panelB <- readRDS("fig_out/05.non_food/Panel_B.rds") + 
  guides(fill = guide_legend(ncol = 1), col = guide_legend(ncol = 1), shape = guide_legend(ncol = 1))
panelC <- readRDS("fig_out/05.non_food/Panel_C.rds") + 
  scale_y_continuous(limits = c(-5,75))+
  guides(fill = guide_legend(ncol = 1), col = guide_legend(ncol = 1), shape = guide_legend(ncol = 1))
panelD <- readRDS("fig_out/05.non_food/Panel_D.rds") + 
  guides(fill = guide_legend(ncol = 1), col = guide_legend(ncol = 1), shape = guide_legend(ncol = 1))
panelE <- readRDS("fig_out/05.non_food/Panel_E.rds") + 
  guides(fill = guide_legend(ncol = 1), col = guide_legend(ncol = 1), shape = guide_legend(ncol = 1))

# 2. Define layout rows
row1 <- (panelA + panelB) + 
  plot_layout(widths = c(1, 1), guides = "collect") & 
  plot_annotation(tag_levels = list(c("A", "B"))) & 
  theme(plot.margin = margin(-5,-7,5,-5, "pt"),
        plot.tag = element_text(size = 20, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 9),
        legend.key.size = unit(4, "mm"))
row2 <- (panelC + panelD) + 
  plot_layout(widths = c(2, 3), guides = "collect") & 
  plot_annotation(tag_levels = list(c("C", "D"))) & 
  theme(plot.margin = margin(-15,-7,5,-5, "pt"),
        plot.tag = element_text(size = 20, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 9),
        legend.key.size = unit(4, "mm"))
row3 <- panelE + 
  plot_annotation(tag_levels = list(c("E"))) & 
  theme(plot.margin = margin(-20,15,5,-5, "pt"),
        plot.tag = element_text(size = 20, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 9),
        legend.key.size = unit(4, "mm"))

# 3. Assemble final figure
figure5 <- (wrap_elements(row1) / wrap_elements(row2) / wrap_elements(row3)) + 
  plot_layout(heights = c(5,5,4)) +
  theme(
    plot.title = element_text(size=9),
    plot.subtitle = element_text(size=8),
    plot.caption = element_text(size=8),
    axis.text = element_text(size = 8),
    panel.spacing = unit(2, "mm")
  )

figure5


figout="fig_out/05.non_food/Figure5_collated"
for(dev in c("svg","png","pdf")){
  outpath=paste0(figout,".",dev)
  ggsave(plot = figure5, filename = outpath, device = dev, 
         width = 210, height = 260, units = "mm", dpi = 300, scale=1.1)
}


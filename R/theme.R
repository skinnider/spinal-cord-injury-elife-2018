library(tidyverse)

# color palettes
ryb8 = c('#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1',
         '#4575b4')
flat.ui = c("#2C3E50", "#E74C3C", "#ECF0F1", "#3498DB", "#2980B9")
jordan = c("#F58D62", "#5EC5E2", "#8F6CA9", "#1482B2", "#F37F81")

# function to programatically darken colors
darken = function(color, factor = 1.4) {
  col = col2rgb(color)
  col = col / factor
  col = rgb(t(col), maxColorValue = 255)
  col
}

# ggplot theme
sci_theme = theme_bw() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        panel.grid = element_blank(),
        panel.border = element_blank(), 
        axis.line.y = element_line(colour = "grey50"),
        axis.line.x = element_line(colour = "grey50"), 
        axis.ticks = element_line(colour = "grey50"),
        legend.position = "top",
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.4, "cm"),
        plot.title = element_text(size = 10, hjust = 0.5))

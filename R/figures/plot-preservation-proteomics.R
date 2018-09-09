# Plot module preservation at the proteomic level. 
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(RColorBrewer)
source("R/theme.R")

# read module conservation
pres = read.delim("data/modules/preservation.tsv")
# remove modules '0' and '0.1' and unused columns
pres = pres[pres$module != 0 & pres$module != 0.1, 1:3]
# categorize Z
pres$Z = cut(pres$Zsummary.pres, breaks=c(-5, 1.9, 4.9, 9.9, 50), 
             labels = c("Z < 2", "Z > 2", "Z > 5", "Z > 10"))

# read preservation at the proteomic level
prot = read.delim("data/proteomics/preservation.txt")
prot %<>% dplyr::select(module, Zsummary.pres)

# plot composite
prot %<>%
  dplyr::filter(!module %in% c(0, 0.1)) %>%
  dplyr::mutate(Z = cut(Zsummary.pres, breaks=c(-5, 1.9, 4.9, 9.9, 50), 
                        labels=c("Z < 2", "Z > 2", "Z > 5", "Z > 10")))
p3 = ggplot(prot, aes(x = module, y = Zsummary.pres, fill = Z)) + 
  geom_hline(aes(yintercept = 2), linetype = 'dotted', color = 'grey50') +
  geom_hline(aes(yintercept = 5), linetype = 'dotted', color = 'grey50') +
  geom_hline(aes(yintercept = 10), linetype = 'dotted', color = 'grey50') +
  geom_col() + 
  geom_hline(aes(yintercept = 0), color = 'grey50') +
  scale_fill_manual(name = "Zsummary", values=ryb8[5:8]) +
  scale_x_continuous(breaks = 1:15, expand = c(0, 0)) +
  scale_y_continuous(expression(Z[summary] * ", proteome")) +
  sci_theme + 
  theme(axis.line.x = element_blank(),
        axis.line.y = element_line(color = 'grey50'), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'top')
p3
# stack on top of heatmap
p4 = ggplot(pres, aes(x=module, y=array, fill=Z)) + 
  geom_tile() + 
  scale_x_continuous(breaks = c(1:15), labels = c(1:15), expand = c(0, 0)) + 
  scale_y_discrete(labels = c("Rat", "Human", "Mouse"), expand = c(0, 0)) +
  scale_fill_manual(name = "Zsummary", values = ryb8[5:8], guide = F) +
  labs(x = "Module", y = "Array") + 
  sci_theme + 
  theme(axis.line = element_blank(), 
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "top",
        axis.title.x = element_blank())
p4

# plot together 
library(gtable)
library(grid)
g1 = ggplotGrob(p3)
g2 = ggplotGrob(p4)
g = rbind(g1, g2, size = "first")
g$widths = unit.pmax(g1$widths, g2$widths)
panels = g$layout$t[grep("panel", g$layout$name)]
g$heights[panels] = unit(c(2.5, 1), "null")
grid.newpage()
grid.draw(g)
ggsave("figures/figure-5D.pdf", g, width = 8.63, height = 7.2, units = 'cm')

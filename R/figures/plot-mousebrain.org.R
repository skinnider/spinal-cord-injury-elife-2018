# Plot enrichment of human spinal cord modules for cell types of the mouse
# spinal cord in scRNA-seq data.
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(tidyverse)
library(RColorBrewer)
source("R/theme.R")

# read cell type enrichments 
datasets = list.files("data/expression/mousebrain.org", pattern = "pSI_*",
                      full.names = T)
names(datasets) = paste("Level", 1:3)

# load, combine datasets, set names
enrs = datasets %>%
  map(read.delim) %>%
  map2_df(names(datasets), ~ mutate(.x, ID = .y))

# MHT across modules
enrs$p = enrs$p * n_distinct(enrs$module)
enrs$p = ifelse(enrs$p > 1, 1, enrs$p)

# plot -log(p)
enrs$logP = -log(enrs$p)

# fix cell names
enrs %<>% mutate(cell = gsub("\\.\\.", ", ", cell),
                 cell = gsub("\\.", " ", cell)) %>%
  # cap -log(P)
  mutate(logP = ifelse(logP > 50, 50, logP))

# plot
p = ggplot(enrs, aes(y = factor(module), x = factor(cell), fill = logP,
                     color = ID)) +
  geom_tile(color = "white") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_distiller(palette = "OrRd", direction = 1,
                       labels = c(seq(0, 40, 10), ">50")) +
  labs(fill = "-log(P)", y = "Module") +
  facet_wrap(~ID, ncol = 3, scales = "free") +
  # scale_fill_distiller(palette = "RdBu", direction = 1) +
  sci_theme +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        legend.title = element_text(face = "bold"),
        strip.background = element_rect(fill = 'grey40', color = 'grey40'),
        strip.text = element_text(color = 'white', size = 8))
p
ggsave("figures/figure-S7.pdf", p, width = 17.4, height = 10, 
       units = "cm", useDingbats = F)

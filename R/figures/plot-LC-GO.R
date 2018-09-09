# Plot enrichment of LC genes for shared GO terms. 
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
source("R/theme.R")

# plot
results = read.delim("data/GO/LC-GO-enrichment.txt") %>%
  mutate(breadth = factor(breadth))
p = ggplot(results, aes(x = breadth, y = enrichment, fill = ontology,
                        color = ontology)) + 
  geom_col(position = 'dodge', width = 0.8) + 
  labs(x = "Breadth", y = "Odds ratio") + 
  scale_fill_manual(values = ryb8[c(1, 5, 8)], name = "Ontology") +
  scale_color_manual(values = darken(ryb8[c(1, 5, 8)], 1.3), 
                     name = "Ontology") +
  scale_y_log10(expand = c(0.01, 0)) +
  sci_theme
p
ggsave("figures/figure-2C.pdf", p, width = 7.8, height = 5,
       units = "cm", useDingbats = F)

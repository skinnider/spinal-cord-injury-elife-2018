# Plot the number of co-complexed LC genes vs. random sets of complexes.
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
source("R/theme.R")

# read data
results = read.delim("data/PPI/co-complex.txt")
rnd = filter(results, class == 'shuffled')
obs = filter(results, class == 'observed')
p = ggplot(rnd, aes(x = co_complex)) + 
  geom_vline(data = obs, aes(xintercept = co_complex), linetype = 'dotted') +
  geom_density(fill = jordan[4], color = jordan[4], alpha = 0.6) +
  scale_y_continuous(expand = c(0, 0)) +
  # scale_x_continuous(limits = c(230, 705)) +
  xlab("Co-complex Memberships") +
  sci_theme + 
  theme(axis.line.y = element_blank(),
        axis.line.x = element_line(color = 'grey50'),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())
p
ggsave("figures/figure-S1H.pdf", p, height = 3.5, width = 4, 
       units = "cm", useDingbats = F)

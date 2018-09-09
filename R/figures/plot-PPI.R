# Plot the number of protein-protein interactions (PPIs) and the size of the
# largest connected component (LCC) formed by LC genes in human interactomes,
# compared to rewired networks.
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
source("R/theme.R")

# read data
results = read.delim("data/PPI/PPI-LCC.txt")

# plot each interactome separately 
interactomes = unique(results$interactome)
for (interactome_name in interactomes) {
  subset = filter(results, interactome == interactome_name)
  rnd = filter(subset, class == "rewired")
  obs = filter(subset, class == "observed")
  
  # plot LCC
  p1 = ggplot(rnd, aes(x = LCC)) + 
    geom_vline(data = obs, aes(xintercept = LCC), linetype = "dotted") +
    geom_density(aes(fill = '1', color = '1')) +
    scale_fill_manual(values = ryb8[7], guide = F) +
    scale_color_manual(values = c(darken(ryb8[7], 1.3), "grey70"), guide = F) +
    scale_y_continuous(expand = c(0.01, 0)) +
    scale_linetype_manual(values = c("dotted"), guide = F) +
    xlab("Largest connected component") +
    sci_theme + 
    theme(axis.line.x = element_line(), 
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank())
  p1
  output = ifelse(interactome_name == "Menche-2015",
                  "figures/figure-2E.pdf", 
                  paste0("figures/figure-S1-LCC-", interactome_name, ".pdf"))
  ggsave(output, p1, height = 3.5, width = 4, units = 'cm', useDingbats = F)
  
  # plot PPIs
  p2 = ggplot(rnd, aes(x = PPI)) + 
    geom_vline(data = obs, aes(xintercept = PPI), linetype = "dotted") +
    geom_density(aes(fill = '1', color = '1')) +
    scale_fill_manual(values = ryb8[5], guide = F) +
    scale_color_manual(values = c(darken(ryb8[5], 1.3), "grey70"), guide = F) +
    scale_y_continuous(expand = c(0.01, 0)) +
    scale_linetype_manual(values = c("dotted"), guide = F) +
    xlab("Protein-protein interactions") +
    sci_theme + 
    theme(axis.line.x = element_line(), 
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank())
  p2
  output = ifelse(interactome_name == "Menche-2015",
                  "figures/figure-2D.pdf", 
                  paste0("figures/figure-S1-PPI-", interactome_name, ".pdf"))
  ggsave(output, p2, height = 3.5, width = 4, units = 'cm', useDingbats = F)
}

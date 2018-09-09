# Plot the overlap between literature-curated genes identified by different
# experimental techniques, injury models, and model organisms, and the M3/M7
# genes.
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
source("R/theme.R")
library(RColorBrewer)
library(magrittr)
library(drlib)
library(cowplot)

# read LC genes
proteins = read.delim("data/literature_review/Table_S1.txt") %>%
  # fix experimental techniques and models
  mutate(method = fct_recode(
    technique.subtype,
    "western blot" = "WB",
    "western blot" = "immunoblot",
    "ELISA" = "immunosorbent assay",
    "MS/MS" = "LC-MS/MS",
    "MS/MS" = "NanoLC-HR- MS/MS",
    "MS/MS" = "NanoLC-HR-MS/MS",
    "MS/MS" = "CAX-PAGELC-MS/MS",
    "expression profiling" = "microarray"),
    model = fct_recode(
      type,
      "balloon-compression" = "balloon compression",
      "ischemia" = "ischemia/repurfusion",
      "clip-compression" = "clip compression",
      "balloon-compression" = "aorta occlusion"))

# read modules
modules = read.delim("data/modules/GTEx-modules.tsv")

# count genes in M3 and M7 by method
methods = proteins %>%
  mutate(M3 = ifelse(ensembl %in% modules$gene[modules$module == 3], 1, 0)) %>%
  mutate(M7 = ifelse(ensembl %in% modules$gene[modules$module == 7], 1, 0)) %>%
  select(ensembl, method, M3, M7) %>%
  group_by(method) %>%
  distinct(ensembl, .keep_all = T) %>%
  select(-ensembl) %>%
  summarise(M3 = sum(M3) / n(), M7 = sum(M7) / n()) %>%
  gather(Module, value, -method)

# adjust order
colors = colorRampPalette(rev(brewer.pal(9, "OrRd")))(5)

# plot
p1 = ggplot(methods, aes(x = reorder_within(method, -value, Module, fun = mean),
                         y = value, fill = Module, order = Module, 
                         color = Module)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.85) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = darken(colors)) +
  scale_y_continuous("% in module", expand = c(0.02, 0),
                     labels = function(breaks) breaks * 100) +
  scale_x_reordered() +
  facet_wrap(~ Module, scales = "free") +
  guides(color = F, fill = F) +
  sci_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
      axis.line = element_line(colour = "grey50"),
        axis.title.x = element_blank(),
        strip.background = element_rect(fill = 'grey50', color = 'grey50'),
        strip.text = element_text(color = 'white', size = 8))
p1

# count genes in M3 and M7 by model
methods = proteins %>%
  mutate(M3 = ifelse(ensembl %in% modules$gene[modules$module == 3], 1, 0)) %>%
  mutate(M7 = ifelse(ensembl %in% modules$gene[modules$module == 7], 1, 0)) %>%
  select(ensembl, model, M3, M7) %>%
  group_by(model) %>%
  distinct(ensembl, .keep_all = T) %>%
  select(-ensembl) %>%
  summarise(M3 = sum(M3) / n(), M7 = sum(M7) / n()) %>%
  gather(Module, value, -model)

# plot
p2 = ggplot(methods, aes(x = reorder_within(model, -value, Module, fun = mean),
                          y = value, fill = Module, order = Module, 
                          color = Module)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.85) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = darken(colors)) +
  scale_y_continuous("% in module", expand = c(0.02, 0),
                     labels = function(breaks) breaks * 100) +
  scale_x_reordered() +
  facet_wrap(~ Module, scales = "free") +
  guides(color = F, fill = F) +
  sci_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(colour = "grey50"),
        axis.title.x = element_blank(),
        strip.background = element_rect(fill = 'grey50', color = 'grey50'),
        strip.text = element_text(color = 'white', size = 8))
p2

# count genes in M3 and M7 by species
species = proteins %>%
  mutate(M3 = ifelse(ensembl %in% modules$gene[modules$module == 3], 1, 0)) %>%
  mutate(M7 = ifelse(ensembl %in% modules$gene[modules$module == 7], 1, 0)) %>%
  select(ensembl, species, M3, M7) %>%
  group_by(species) %>%
  distinct(ensembl, .keep_all = T) %>%
  select(-ensembl) %>%
  summarise(M3 = sum(M3) / n(), M7 = sum(M7) / n()) %>%
  gather(Module, value, -species) %>%
  mutate(species = factor(species, levels = c(9606, 9986, 10090, 10116),
                          labels = c("human", "rabbit", "mouse", "rat")))

# plot 
p3 = ggplot(species, aes(x = reorder_within(species, -value, Module, fun = mean),
                         y = value, fill = Module, order = Module, 
                         color = Module)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.85) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = darken(colors)) +
  scale_y_continuous("% in module", expand = c(0.02, 0),
                     labels = function(breaks) breaks * 100) +
  scale_x_reordered() +
  facet_wrap(~ Module, scales = "free") +
  guides(color = F, fill = F) +
  sci_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(colour = "grey50"),
        axis.title.x = element_blank(),
        strip.background = element_rect(fill = 'grey50', color = 'grey50'),
        strip.text = element_text(color = 'white', size = 8))
p3

# combine plots
p = p1 + p2 + p3 + plot_layout(widths = c(1.8, 1.4, 1))
p
ggsave("figures/figure-S2.pdf", p, width = 17.4, height = 5, units = 'cm',
       useDingbats = F)

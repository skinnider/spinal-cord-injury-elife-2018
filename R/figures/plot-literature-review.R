# Plot figure panels from Figure 1 and Figure S1 relating to the literature
# review of SCI genes.
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
source("R/theme.R")

# read data
dat = read.delim("data/literature_review/Table_S1.txt")

# get summary statistics
n_genes = n_distinct(dat$ensembl)
n_upregulated = dat %>% filter(direction == 'upregulated') %>% 
  pull(ensembl) %>% n_distinct()
n_downregulated = dat %>% filter(direction == 'downregulated') %>% 
  pull(ensembl) %>% n_distinct()
n_phosphorylated = dat %>% filter(grepl("phosphorylated", direction)) %>% 
  pull(ensembl) %>% n_distinct()
n_multistudy = dat %>% group_by(ensembl) %>% filter(n_distinct(pmid) > 1) %>% 
  pull(ensembl) %>% n_distinct()

# plot Figure 2A: # of studies
xvals = seq_len(5)
pmids = dat %>% 
  group_by(ensembl) %>%
  summarise(pmids = n_distinct(pmid)) %>%
  ungroup()
n_studies = data.frame(studies = xvals, genes = map_int(
  xvals, ~ sum(pmids$pmids >= .))) %>%
  mutate(studies = factor(studies))
p = ggplot(n_studies, aes(x = studies, y = genes, fill = studies,
                          color = studies)) + 
  geom_col(width = 0.7) + 
  scale_fill_manual(values = ryb8[c(1:4, 4)], guide = F) +
  scale_color_manual(values = darken(ryb8[c(1:4, 4)], 1.3), guide = F) +
  labs(x = "Number of studies", y = "Genes") +
  scale_y_continuous(expand = c(0.01, 0)) +
  scale_x_discrete(expand = c(0.15, 0), labels = c(1, 2, 3, 4, "5+")) +
  sci_theme
p
ggsave("figures/figure-2A.pdf", p, width = 3.7, height = 4.5,
       units = "cm", useDingbats = F)

# plot Figure 2B: experimental methods
dat %<>% mutate(method = fct_recode(
  technique.subtype, 
  "western blot" = "WB",
  "western blot" = "immunoblot",
  "ELISA" = "immunosorbent assay",
  "MS/MS" = "LC-MS/MS",
  "MS/MS" = "NanoLC-HR- MS/MS",
  "MS/MS" = "NanoLC-HR-MS/MS",
  "MS/MS" = "CAX-PAGELC-MS/MS",
  "expression profiling" = "microarray"))
methods = dat %>% 
  group_by(method) %>%
  summarise(genes = n_distinct(ensembl)) %>%
  ungroup() %>%
  arrange(-genes)
p = ggplot(methods, aes(x = reorder(method, -genes), y = genes, fill = method,
                        color = method, order = genes)) +
  geom_col(width = 0.85) +
  scale_fill_manual(values = setNames(ryb8[c(1:8, 8, 8)], methods$method), 
                    guide = F) +
  scale_color_manual(values = setNames(darken(ryb8[c(1:8, 8, 8)], 1.3), 
                                       methods$method), guide = F) +
  ylab("Genes") +
  scale_y_continuous(expand = c(0.01, 0)) +
  scale_x_discrete(expand = c(0.08, 0)) +
  sci_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(colour = "grey50"),
        axis.title.x = element_blank())
p
ggsave("figures/figure-2B.pdf", p, width = 6.5, height = 5.5,
       units = "cm", useDingbats = F)

# plot Figure S1A: experimental models
dat %<>% mutate(model = fct_recode(
  type, 
  "balloon-compression" = "balloon compression",
  "ischemia" = "ischemia/repurfusion",
  "clip-compression" = "clip compression",
  "balloon-compression" = "aorta occlusion"))
models = dat %>%
  group_by(model) %>%
  summarise(genes = n_distinct(ensembl)) %>%
  arrange(-genes)
p = ggplot(models, aes(x = reorder(model, -genes), y = genes, fill = model, 
                      order = model, color = model)) +
  geom_col(width = 0.85) +
  scale_fill_manual(values = setNames(
    ryb8[c(1, 2, 2, 3, 3, 4, 4)], models$model), guide = F) +
  scale_color_manual(values = setNames(
    darken(ryb8[c(1, 2, 2, 3, 3, 4, 4)], 1.3), models$model), guide = F) +
  ylab("Genes") +
  scale_y_continuous(expand=c(0.01, 0)) +
  sci_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(colour = "grey50"),
        axis.title.x = element_blank())
p
ggsave("figures/figure-S1A.pdf", p, width = 4.5, height = 5.5,
       units = "cm", useDingbats = F)

# plot Figure S1B: species
species = dat %>% 
  group_by(species) %>% 
  summarise(genes = n_distinct(ensembl)) %>%
  arrange(-genes) %>%
  mutate(species = factor(species))
p = ggplot(species, aes(x = reorder(species, -genes), y = genes, fill = species,
                        color = species)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = setNames(
    ryb8[c(1:4)], species$species), guide = F) +
  scale_color_manual(values = setNames(
    darken(ryb8[c(1:4)], 1.3), species$species), guide = F) +
  ylab("Genes") +
  scale_y_continuous(expand = c(0.01, 0)) +
  scale_x_discrete(expand = c(0.02 ,0),
                   labels = c("rat", "human", "mouse", "rabbit")) +
  sci_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(colour = "grey50"),
        axis.title.x = element_blank())
p
ggsave("figures/figure-S1B.pdf", p, width = 3.5, height = 4.5,
       units = "cm", useDingbats = F)

# plot Figure S1C: time points
dat %<>% mutate(time = fct_recode(
  time,
  '1d' = '24hr',
  '7d' = '1-8d',
  '1hr' = '1h',
  '8hr' = '8h',
  '21d' = '3w',
  '56d' = '8w'))
times = dat %>%
  mutate(time = factor(time, levels = c(
    "1hr", "2hr", "4hr", "6hr", "8hr", "12hr", "1d", "3d", "4d", "5d", "7d",
    "8d", "10d", "12d", "14d", "21d", "28d", "35d", "42d", "46d", "56d"))) %>%
  drop_na() %>%
  group_by(time) %>%
  summarise(genes = n_distinct(ensembl)) %>%
  ungroup() %>%
  arrange(time)
# tag phase: acute, subacute, chronic
acute = c("1hr", "2hr", "4hr", "6hr", "8hr", "12hr", "1d")
subacute = c("3d", "4d", "5d", "7d", "8d", "10d", "12d", "14d")
chronic = c("21d", "28d", "35d", "42d", "46d","56d")
times %<>% 
  mutate(phase = ifelse(time %in% acute, "acute", ifelse(
    time %in% subacute, "subacute", "chronic")))
p = ggplot(times, aes(x = time, y = genes, fill = phase, color = phase)) + 
  geom_col() + 
  scale_fill_manual(values = ryb8[c(1, 7, 5)], name = "Phase") +
  scale_color_manual(values = darken(ryb8[c(1, 7, 5)], 1.3), name = "Phase") +
  labs(x = "Time point", y = "Genes") + 
  scale_y_continuous(expand = c(0.01, 0)) +
  scale_x_discrete(expand = c(0.01, 0)) + 
  sci_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.line = element_line(colour = "grey50"))
p
ggsave("figures/figure-S1C.pdf", p, width = 10, height = 5,
       units = "cm", useDingbats = F)

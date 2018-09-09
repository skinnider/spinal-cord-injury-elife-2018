# Plot the most accurate hub gene biomarkers in our RNA-seq datset.
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(org.Hs.eg.db)
source("R/theme.R")

# read data
dat = read.delim("data/biomarker/hubs-rnaseq.txt")

# code moderate/severe
dat %<>% dplyr::mutate(severity = ifelse(
  severity == 100, "mild/moderate", "severe"))

# map to gene names
map = AnnotationDbi::select(
  org.Hs.eg.db, keys = unique(dat$hub), keytype = 'ENSEMBL',
  columns = c("SYMBOL", "GENENAME")) %>%
  dplyr::rename(hub = ENSEMBL) %>% 
  group_by(hub) %>%
  dplyr::slice(1)
dat %<>% left_join(map, by = 'hub') 

# what biomarkers have accuracy > 0.9 in RNAseq?
acc90 = dat %>% 
  dplyr::filter(accuracy >= 0.9) %>% 
  group_by(hub) %>%
  dplyr::filter(n() == 2) %>%
  ungroup()

# plot biomarkers with accuracy > 0.9 in RNAseq
acc90$severity %<>% replace(., . == 'mild/moderate', 'moderate')
p = ggplot(acc90, aes(x = reorder(SYMBOL, -accuracy), y = accuracy, 
                      fill = severity)) + 
  geom_col() + 
  facet_grid(~ severity) + 
  scale_fill_brewer('RdBu', guide = F) + 
  scale_y_continuous("Accuracy", expand = c(0, 0), limits = c(0, 1.05)) +
  sci_theme +
  theme(axis.text.x = element_text(hjust = 1, angle = 45),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_line(color = 'grey50'),
        strip.text = element_text(size = 8))
p
ggsave("figures/figure-5I.pdf", p, width = 6, height = 4.5, units = 'cm')

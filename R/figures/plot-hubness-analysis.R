# Plot the relationship between correlation to the M3 eigengene and 
# predictive accuracy as a SCI biomarker. 
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(cowplot)
source("R/theme.R")

# read sleuth 
sleuth = read.delim("data/expression/rnaseq/sleuth/sleuth-severity.txt") %>%
  dplyr::rename(gene = ens_human)

# read module eigengenes
kME = read.delim("data/modules/kME.txt") %>%
  rownames_to_column('gene') %>% 
  dplyr::filter(gene %in% sleuth$gene)

# read module memberships
modules = read.delim("data/modules/GTEx-modules.tsv")

# is the 'hubness' of M3 genes a good indicator of severity?
all = sleuth %>%
  left_join(kME, by = 'gene')
m3_genes = modules$gene[modules$module == 3]
m3 = all %>%
  dplyr::filter(gene %in% m3_genes)
cor = cor.test(m3$M3, m3$b, method = 'spearman')

# plot
dat = all %>%
  dplyr::mutate(in_M3 = gene %in% m3_genes)
pal = ryb8[c(6, 1)]
base = ggplot(dat, aes(x = M3, y = b, color = in_M3, group = in_M3)) + 
  geom_point(size = 0.05) + 
  geom_smooth(data = m3, alpha = 0, method = 'lm',
              aes(color = gene %in% m3_genes, group = gene %in% m3_genes)) + 
  scale_x_continuous("M3 hubness") + 
  scale_y_continuous(expression(paste(beta, ", SCI severity"))) + 
  scale_color_manual("Module", values = pal, labels = c("Other", "M3")) +
  sci_theme
base

# add densities
y_density <- axis_canvas(base, axis = "y", coord_flip = TRUE) +
  geom_density(data = dat, aes(x = b, fill = in_M3, color = in_M3),
               alpha = 0.5) +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  coord_flip()
x_density <- axis_canvas(base, axis = "x") +
  geom_density(data = dat, aes(x = M3, fill = in_M3, color = in_M3), 
               alpha = 0.5) +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal)

# create the combined plot
combined = base %>%
  insert_yaxis_grob(., y_density, position = "right") %>%
  insert_xaxis_grob(., x_density, position = "top")

# show the result
ggdraw(combined)

# save
ggsave("figures/figure-5C.pdf", width = 5.5, height = 6, 
       units = 'cm', useDingbats = F)

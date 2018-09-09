# Plot expression of ANXA1 in our RNA-seq and proteomics data.
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(org.Hs.eg.db)
source("R/theme.R")

# read RNA-seq
rnaseq = read.delim("data/expression/rnaseq/sleuth/sleuth-norm-filt.txt")
targets = read.delim("data/expression/rnaseq/targets.txt")
rnaseq = as.data.frame(cbind(targets, t(rnaseq)))
# plot annexin A1
anxa1 = "ENSG00000135046"
rnaseq$label %<>% Hmisc::capitalize() %>% 
  factor(levels = c('Sham', 'Moderate', 'Severe'))
rnaseq$title = "ANXA1"
p = ggplot(rnaseq, aes(x = label, y = ENSG00000135046, color = label, 
                       fill = label)) + 
  facet_grid(~ title) +
  geom_boxplot(alpha = 0.4, outlier.size = 0, outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.2, size = 0.4) +
  scale_x_discrete("Injury") + # , labels = c("Sham", "Mild", "Severe")) +
  scale_y_continuous("TPM") +
  scale_color_manual(values = ryb8[c(1, 3, 8)], guide = F) +
  scale_fill_manual(values = ryb8[c(1, 3, 8)], guide = F) +
  sci_theme + 
  theme(axis.line = element_line(color = 'grey50'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(), 
        strip.background = element_rect(fill = 'grey40', color = 'grey40'),
        strip.text = element_text(color = 'white', size = 8, face = 'italic'))
p
ggsave("figures/figure-5J.pdf", p, width = 4, height = 4.5, 
       units = 'cm', useDingbats = F)

# read proteomics
prot = read.delim("data/proteomics/proteomics-processed.txt") %>%
  column_to_rownames('Ortholog')
names = list(rownames(prot), colnames(prot))
prot %<>% as.matrix() %>%
  preprocessCore::normalize.quantiles()
dimnames(prot) = names
targets = read.delim("data/proteomics/targets.txt")
prot = as.data.frame(cbind(t(prot), targets))
# plot annexin A1
anxa1 = "ENSG00000135046"
prot$title = "ANXA1"
prot$label %<>% Hmisc::capitalize() %>% factor(
  levels = c('Sham', 'Moderate', 'Severe'))
scientific_10 <- function(x) {
  parse(text = gsub("\\+", "", gsub(
    "e", " %*% 10^", scales::scientific_format()(x))))
}
p = ggplot(prot, aes(x = label, y = ENSG00000135046, color = label, 
                       fill = label)) + 
  facet_grid(~ title) + 
  geom_boxplot(alpha = 0.4, outlier.size = 0, outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.2, size = 0.4) +
  scale_x_discrete("Injury") +
  scale_y_continuous("Abundance", label = scientific_10) +
  scale_color_manual(values = ryb8[c(1, 3, 8)], guide = F) +
  scale_fill_manual(values = ryb8[c(1, 3, 8)], guide = F) +
  sci_theme + 
  theme(axis.line = element_line(color = 'grey50'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(), 
        strip.background = element_rect(fill = 'grey40', color = 'grey40'),
        strip.text = element_text(color = 'white', size = 8))
p
ggsave("figures/figure-5K.pdf", p, width = 4.3, height = 4.5, 
       units = 'cm', useDingbats = F)

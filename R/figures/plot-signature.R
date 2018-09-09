# Plot the consensus network signature of SCI, and contrast it with 
# RNA-seq, proteomics, and NT-3. 
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
source("R/theme.R")

# read consensus signature
results = read.delim("data/expression/DE/consensus.txt") %>%
  dplyr::filter(experiment == 'consensus')

# read RNA-seq data from sleuth
rnaseq = read.delim("data/expression/rnaseq/sleuth/sleuth-severity-module-regulation.txt")
rnaseq = rnaseq %>%
  dplyr::select(module, p, sig, direction) %>%
  dplyr::mutate(experiment = "RNA-seq") %>%
  dplyr::filter(sig) 

# read proteomics data 
prot = read.delim("data/proteomics/regulation.txt")
prot = prot %>%
  dplyr::select(module, p, sig, direction) %>%
  dplyr::mutate(experiment = "proteomics") %>%
  filter(sig)

# read NT-3 data
nt3 = read.delim("data/expression/DE/E-GEOD-69334/regulation.tsv") %>%
  dplyr::filter(condition == 'treatmentTmiddle' & sig) %>%
  dplyr::mutate(experiment = 'NT-3') %>%
  dplyr::select(module, p, sig, direction, experiment)

# read Anderson data
anderson = read.delim("data/expression/geo/GSE76097/regulation.tsv") %>%
  filter(sig) %>%
  mutate(experiment = "reduced axonal\ndieback") %>%
  dplyr::select(module, p, sig, direction, experiment)

# combine all
dat = bind_rows(results, rnaseq, prot, nt3, anderson)

# add n.s. for other experiments
experiments = setdiff(unique(dat$experiment), 'consensus')
for (experiment in experiments) {
  for (i in unique(dat$module)) {
    if (sum(dat$module == i & dat$experiment == experiment) == 0) {
      dat %<>% rbind(list(i, 0, F, "n.s.", experiment))
    }
  }
}

# set experiment as factor
dat$experiment = factor(dat$experiment, levels = c(
  'reduced axonal\ndieback', 'NT-3', 'proteomics', 'RNA-seq', 'consensus'))
# plot
p = ggplot(dat, aes(x = module, y = experiment, fill = direction)) +
  geom_tile() +
  scale_x_continuous("Module", breaks=c(1:15), labels=c(1:15),
                     expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_manual(name = "Direction", values=ryb8[c(1, 5, 8)]) +
  labs(x="Module", y="Experiment") +
  sci_theme +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "top")
p
ggsave("figures/figure-5B.pdf", p, width = 8.75, height = 4.6,
       useDingbats = F, units = "cm")

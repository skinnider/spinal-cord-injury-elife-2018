# Plot enrichment of each module for literature-curated SCI proteins.
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
source("R/theme.R")
source("R/expression/wgcna/wgcna-functions.R")
source("R/expression/DE/DE-functions.R")

# read modules
modules = read.delim("data/modules/GTEx-modules.tsv")
module.genes = getModuleGenes(modules)

# read SCI proteins
sci = read.delim("data/literature_review/Table_S1.txt")
sci.genes = unique(sci$ensembl)

# find significant modules 
enrich = enr(sci.genes, modules, correct="bonferroni")
enrich$logp = -log10(enrich$enrich.p)
enrich$sig = enrich$enrich.p < 0.05

# plot
p = ggplot(enrich, aes(x = module, y = logp, fill = sig, colour = sig)) +
  geom_bar(stat="identity", size = 0.5) + 
  geom_hline(aes(yintercept = -log10(0.05)), linetype = "dotted") +
  scale_fill_manual("Module", values = ryb8[c(5, 8)], guide = F) +
  scale_color_manual("Module", values = darken(ryb8[c(5, 8)], 1.3), guide = F) +
  scale_x_continuous(expand=c(0,0), breaks=c(1:15)) +
  scale_y_continuous(expand=c(0.01,0)) +
  labs(x="Module", y="-log(P)") + 
  sci_theme + 
  theme(legend.position = "right", 
        axis.line = element_line(color = "grey50"))
p
ggsave("figures/figure-3B.pdf", p, width = 5.6, height = 4.8, 
       units = "cm", useDingbats = F)

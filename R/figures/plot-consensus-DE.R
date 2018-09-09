# Plot consensus upregulated and downregulated modules in spinal cord injury
# across experimentla mouse and rat datasets. 
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
source("R/theme.R")

# description of experiments:
# - E-GEOD-464: rat U34 data (Di Giovanni et al, Ann Neurol 2003), 542 samples
# - E-GEOD-5296: unpublished mouse dataset, 90 samples
# - E-GEOD-45006: rat 230 data (Chamankhah et al, BMC Genomics 2013), 24 samples
# - E-GEOD-45376: mouse RNA-seq data (Chen et al, PLoS ONE 2013), 8 samples
# - E-GEOD-69334: rat 230 data (Duan et al, PNAS 2015), 167 samples 
experiments = c("E-GEOD-464", "E-GEOD-5296", "E-GEOD-45006", "E-GEOD-45376", 
                "E-GEOD-69334")
results = data.frame()
for (experiment in experiments) {
  exptFile = paste("data/expression/DE", experiment, "regulation.tsv", 
                   sep = "/")
  reg = read.delim(exptFile)
  if ("condition" %in% colnames(reg))
    if (experiment == "E-GEOD-464") {
      reg = reg[reg$condition == "injuryT9", ]
    } else {
      reg = reg[reg$condition == "injury", ]
    }
  reg = reg[reg$sig, ]
  reg = reg[, c("module", "p", "sig", "direction")]
  reg$experiment = experiment
  results = rbind(results, reg)
}

# add missing data 
for (i in 1:15)
  for (expt in experiments)
    if (nrow(results[results$experiment == expt & results$module == i, ]) == 0)
      results = rbind(results, list(i, 1, FALSE, "n.s.", expt))

# add consensus row
for (i in 1:15) {
  directions = unique(results$direction[results$module == i])
  consensus = ifelse(length(directions) > 1, "n.s.", directions[1])
  results = rbind(results, list(i, 1, FALSE, consensus, "consensus"))
}

# write this data
write.table(results, "data/expression/DE/consensus.txt", quote = F,
            row.names = F, sep = "\t")

# plot heatmap 
expt.labs = c("Consensus", "GSE45006", "GSE45376", "GSE464", "GSE5296", 
              "GSE69334")
p = ggplot(results, aes(x = module, y = experiment, fill = direction)) + 
  geom_tile() + 
  geom_hline(aes(yintercept = 1.5), colour = "grey40", size = 0.5) +
  scale_x_continuous(breaks = 1:15, labels = 1:15, expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0), labels = expt.labs) +
  scale_fill_manual(name = "Direction", values=ryb8[c(1, 5, 8)]) +
  labs(x = "Module", y = "Experiment") + 
  sci_theme + 
  theme(axis.line = element_blank(), 
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "top")
p
ggsave("figures/figure-3E.pdf", p, width = 5.9, height = 4, 
       useDingbats = F, units = "cm")

# plot line chart
results = results[results$experiment != "consensus", ]
results = results[results$module %in% c(1, 2, 3, 6, 7, 8, 11), ]
results$logp = -log(results$p)
results$module = factor(results$module)
p = ggplot(results, aes(x = experiment, y = logp, color = module)) + 
  geom_line(aes(group = module)) + 
  geom_hline(aes(yintercept = -log10(0.05)), linetype = "dotted") +
  scale_color_manual("Module", values = ryb8[-6], labels = 
                       paste("M", c(1, 2, 3, 6, 7, 8, 11), sep="")) +
  scale_x_discrete(labels = expt.labs[-1], expand = c(0,0)) +
  ylab("-log(P)") + 
  sci_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "grey50"),
        legend.position = "right")
p
ggsave("figures/figure-3F.pdf", p, width = 6.1, height = 4.8, 
       units = "cm", useDingbats = F)

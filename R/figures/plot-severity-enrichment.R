# Plot module enrichment for genes correlated with injury severity in GSE464.
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
source("R/theme.R")

# read data 
results = read.delim("data/expression/DE/E-GEOD-464/correlation-enrichment.tsv")
results = results[results$p < 1,]
# add missing data 
for (i in 1:15)
  if (sum(results$module == i) == 0)
    results = rbind(results, list(i, 1, "n.s.", F))
# calculate significance (1% FWER) and -logP
results$sig = results$p < 0.01
results$direction[results$sig == F] = "n.s."
results$logp = -log10(results$p)

# plot
p = ggplot(results, aes(x = module, y = logp, fill = direction, 
                        colour = direction)) +
  geom_bar(stat = "identity", size = 0.5) + 
  geom_hline(aes(yintercept = -log10(0.01)), linetype = "dotted") +
  scale_fill_manual("Module", values = ryb8[c(1, 8, 5)]) +
  scale_color_manual("Module", values = darken(ryb8[c(1, 8, 5)], 1.3)) +
  scale_x_continuous(expand = c(0.01,0), breaks = 1:15) +
  scale_y_continuous(expand = c(0.01,0)) +
  labs(x = "Module", y = "-log(P)") + 
  sci_theme + 
  theme(legend.position = "top", 
        axis.line = element_line(color = "grey50"))
p
ggsave("figures/figure-5A.pdf", width = 5.6, height = 4.8, 
       units = "cm", useDingbats = F)

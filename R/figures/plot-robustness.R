# Plot the robustness of spinal cord module enrichment for LC genes to random
# addition or removal.
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(tidyverse)
source("R/theme.R")

results = read.delim("data/modules/robustness.tsv.gz")
results$module = factor(results$module)
results$logp = -log10(results$p)
results = subset(results, module == 3 | module == 7)
summary = results %>%
  group_by(idx, module, condition) %>%
  summarise(median = median(logp), iqr = IQR(logp)) %>%
  ungroup()
p = ggplot(summary, aes(x = idx, y = median, colour = module, 
                        linetype = condition)) +
  geom_linerange(color = flat.ui[3], aes(ymin = median - iqr, 
                                         ymax = median + iqr)) +
  geom_line() +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  scale_color_manual("Module", values = c(ryb8[c(1,8)]), 
                     labels = c("M3", "M7")) +
  scale_linetype_manual(
    "Condition", labels=c("Random\naddition", "Seed\nremoval"), 
    values = c("solid", "dotted")) +
  scale_x_continuous(expand=c(0,0)) +
  labs(x="Index", y="-log(P)") + 
  sci_theme + 
  theme(legend.position = "right", 
        axis.line = element_line(color = "grey50"))
p
ggsave("figures/figure-3C.pdf", p, width = 6.5, height = 4.8, 
       useDingbats = F, units = "cm")

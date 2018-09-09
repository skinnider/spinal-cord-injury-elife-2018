# Create a heatmap to illustrate module conservation across datasets (human) and
# species (mouse, rat). 
setwd("~/git/spinal-cord-injury-elife-2018")

# read preservation
pres = read.delim("data/modules/preservation.tsv")
# remove modules '0' and '0.1' and unused columns
pres = pres[pres$module != 0 & pres$module != 0.1, 1:3]
# categorize Z
pres$Z = cut(pres$Zsummary.pres, breaks=c(-5, 1.9, 4.9, 9.9, 50), 
             labels=c("Z < 2", "Z > 2", "Z > 5", "Z > 10"))
# plot  
p = ggplot(pres, aes(x = module, y = array, fill = Z)) + 
  geom_tile() + 
  scale_x_continuous(breaks = c(1:15), labels = c(1:15), expand = c(0, 0)) + 
  scale_y_discrete(labels = c("Rat", "Human", "Consensus"), expand = c(0, 0)) +
  scale_fill_manual(name = "Zsummary", values = ryb8[5:8]) +
  labs(x = "Module", y = "Array") + 
  sci_theme + 
  theme(axis.line = element_blank(), 
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "top")
p
ggsave("figures/figure-3A.pdf", p, width = 5.9, height = 3.1, 
       useDingbats = F, units = "cm")

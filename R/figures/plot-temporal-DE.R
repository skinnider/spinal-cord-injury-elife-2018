# Plot module regulation at 1 day, 3 days, and 7 days post-SCI.
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
source("R/theme.R")

# read data
experiments = c("E-GEOD-464", "E-GEOD-5296", "E-GEOD-45006")
contrasts = list(c("E1d", "E3d", "E7d"),
                 c("E24", "E72", "E168", "E672"),
                 c("Ed1", "Ed3", "Ewk1", "Ewk8"))
names(contrasts) = experiments
results = data.frame()
for (experiment in experiments) {
  exptFile = paste("data/expression/DE", experiment, "regulation.tsv", 
                   sep = "/")
  reg = read.delim(exptFile)
  exptContrasts = contrasts[[experiment]]
  reg = reg[reg$condition %in% exptContrasts,]
  reg = reg[reg$sig == T,]
  reg = reg[, c("module", "condition", "p", "sig", "direction")]
  reg$experiment = experiment
  results = rbind(results, reg)
}

# standardize condition labels
labels = paste0("D", c(1, 3, 7, 28))
for (i in seq_len(nrow(results)))
  for (j in 1:4) 
    if (results$condition[i] %in% sapply(contrasts, function(x) x[j]))
      results$condition[i] = labels[j]

# add missing data 
for (i in 1:15)
  for (expt in experiments)
    for (condition in labels)
      if (nrow(results[results$experiment == expt & results$module == i & 
                       results$condition == condition, ]) == 0)
        results = rbind(results, list(i, condition, 1, FALSE, "n.s.", expt))

# plot 
results$condition = factor(results$condition, levels = c(
  "D1", "D3", "D7", "D28"))
results$experiment = fct_recode(results$experiment,
                                "GSE45006" = "E-GEOD-45006",
                                "GSE464" = "E-GEOD-464",
                                "GSE5296" = "E-GEOD-5296")
results = results[!(results$experiment == "GSE464" &
                      results$condition == "D28"),]
plot.df = results[results$experiment != "consensus", ]
p = ggplot(plot.df, aes(x = module, y = experiment, fill = direction)) + 
  geom_tile() + 
  facet_grid(condition ~ ., scales = "free_y", space = "free") + 
  scale_x_continuous(breaks = 1:15, labels = 1:15, expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_manual(name = "Direction", values = ryb8[c(1, 5, 8)]) +
  labs(x = "Module", y = "Experiment") + 
  sci_theme + 
  theme(axis.line = element_blank(), 
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_text(size = 6),
        legend.position = "top")
p
ggsave("figures/figure-3G.pdf", p, width = 6.4, height = 5.8,
       useDingbats = F, units = "cm")

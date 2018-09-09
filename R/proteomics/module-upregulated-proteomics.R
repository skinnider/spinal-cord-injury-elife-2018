# Analyze severity-dependent module regulation in proteomics data from the
# Foster lab experiment (2018). 
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(tidyverse)
library(limma)

# read proteomics data
dat = read.delim("data/proteomics/proteomics-processed.txt")
genes = dat$Ortholog
dat = dat[, -1]
rownames(dat) = genes
# convert to matrix
expr = as.matrix(dat)

# load GTEx modules
modules = read.delim("data/modules/GTEx-modules.tsv")
# load targets
targets = read.delim("data/proteomics/targets.txt")
targets$binary = ifelse(targets$label == "sham", "control", "SCI")

# fit linear severity-dependent design
fit = eBayes(lmFit(expr, model.matrix(~ targets$group)))

# do geneSetTest for module upregulation 
diff = data.frame()
for (module in unique(modules$module)) {
  if (module == 0)
    next
  moduleGenes = modules$gene[modules$module == module]
  moduleOnChip = intersect(moduleGenes, rownames(expr))
  for (alt in c("up", "down")) {
    test = geneSetTest(moduleOnChip, fit$t[, 2], type = "t",
                       alternative = alt, ranks.only = T)
    diff = rbind(diff, list(module, test, alt))
  }
} 
colnames(diff) = c("module", "p", "direction")
diff$p = diff$p * 2 * length(unique(modules$module)) # Bonferroni correction
diff$sig = diff$p < 0.05
diff = diff[order(diff$p), ]

# write output
write.table(diff, "data/proteomics/regulation.txt",
            sep = "\t", quote = F, row.names = F)

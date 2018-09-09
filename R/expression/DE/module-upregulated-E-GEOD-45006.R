# Test whether M3 and M7 genes are differentially expressed following spinal 
# cord injury in the BMC Genomics dataset at E-GEOD-45006 using limma.
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(limma)
source("R/expression/DE/DE-functions.R")

# read modules
modules = read.delim("data/modules/GTEx-modules.tsv")

# load preprocessed expression data
expt = "E-GEOD-45006" # BMC Genomics 2013
array = "A-AFFY-43"
load(paste0("data/expression/DE/", expt, "/", expt, "-human.Rdata"))

# read experiment targets
expt.dir = paste("data/expression/arrayexpress/cel", array, expt, sep = "/")
targets = readTargets(paste("data/expression/DE", expt, "targets.tsv", 
                            sep = "/"))

# create design matrices
f1 = factor(targets$condition) # injury vs. control
f2 = factor(paste0(targets$condition, targets$time)) # as above + time
d1 = model.matrix(~0 + f1)
colnames(d1) = levels(f1)
d2 = model.matrix(~0 + f2)
colnames(d2) = levels(f2)

# define a contrast matrix
cm1 = makeContrasts(injury="SCI-sham", levels = d1)
cm2 = makeContrasts(Ed1 = "SCId1-sham", Ed3 = "SCId3-sham", 
  Ewk1 = "SCIwk1-sham", Ewk2 = "SCIwk2-sham",  Ewk8 = "SCIwk8-sham",
  levels = d2)
matrix.list = list(cm1, cm2)
designs = list(d1, d2)

# run geneSetTest
diff = data.frame(module = integer(), condition = character(), p = numeric(),
                   direction = character())
for (i in 1:length(matrix.list)) {
  design = designs[[i]]
  # fit the linear model
  fit = lmFit(exprEnsembl, design)
  # define matrix
  cont.matrix = matrix.list[[i]]
  # extract linear model fit for contrasts
  fit2 = contrasts.fit(fit, cont.matrix)
  fit2 = eBayes(fit2)
  # finally, do gene set enrichment analysis for each module
  for (module in unique(modules$module)) {
    if (module == 0)
      next 
    moduleGenes = modules$gene[modules$module == module]
    moduleOnChip = intersect(moduleGenes, rownames(exprEnsembl))
    conditions = colnames(fit2$t)
    for (i in seq_len(length(conditions))) {
      condition = conditions[i]
      for (alt in c("up", "down")) {
        test = geneSetTest(moduleOnChip, fit2$t[, i], type = "t", 
                          alternative = alt, ranks.only = T)
        diff = rbind(diff, list(module, condition, test, alt))
      }
    }
  }
  colnames(diff) = c("module", "condition", "p", "direction")
}
diff$p = diff$p * 2 * max(modules$module) # Bonferroni correction
diff$sig = diff$p < 0.05
diff = diff[order(diff$module), ]

# write output 
outFile = paste("data/expression/DE", expt, "regulation.tsv", sep = "/")
write.table(diff, outFile, sep = "\t", quote = F, row.names = F)

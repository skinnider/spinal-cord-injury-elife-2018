# Test whether M3 and M7 genes are differentially expressed following spinal 
# cord injury in E-GEOD-69334 using limma.
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(limma)
library(hgu133a.db)
source("R/expression/DE/DE-functions.R")

# load preprocessed expression data
array = "A-AFFY-43"
expt = "E-GEOD-69334" # time and severity series data
load(paste0("data/expression/DE/", expt, "/", expt, 
            "-human.Rdata"), verbose = T)

# read experiment targets
expt.dir = paste("data/expression/arrayexpress/cel", array, expt, sep = "/")
targets = readTargets(paste("data/expression/DE", expt, "targets.tsv", 
                            sep = "/"))
middle.targets = targets[targets$tissue == "Middle", ]

# create injury matrix
f1 = factor(targets$group)
d1 = model.matrix(~ 0 + f1)
colnames(d1) = levels(f1)
cm1 = makeContrasts(injury = "Control-Uninjured", levels = d1)

# create NT-3 matrix
f2 = factor(targets$group[targets$tissue == "Middle"])
d2 = model.matrix(~ 0 + f2)
colnames(d2) = levels(f2)
cm2 = makeContrasts(treatmentTmiddle = "chitosan-Tube", levels = d2)

# load GTEx modules and seed proteins
modules = read.delim("data/modules/GTEx-modules.tsv")
sci = read.delim("data/literature_review/Table_S1.txt")
genes = unique(sci$ensembl)

# finally, analyze regulation of each module
diff = data.frame(module = integer(), condition = character(), p = numeric(),
                  direction = character())
matrix.list = list(cm1, cm2)
designs = list(d1, d2)
for (matrix_idx in seq_along(matrix.list)) {
  cm = matrix.list[[matrix_idx]]
  design = designs[[matrix_idx]]
  # fit the linear model
  if (matrix_idx == 1) {
    fit = lmFit(exprEnsembl, design)
  } else if (matrix_idx == 2) {
    fit = lmFit(exprEnsembl[, colnames(exprEnsembl) %in% 
                              middle.targets$filename], design)
  }
  fit2 = contrasts.fit(fit, cm)
  fit2 = eBayes(fit2)
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
}
colnames(diff) = c("module", "condition", "p", "direction")
diff$p = diff$p * 2 * max(modules$module) # Bonferroni correction
diff$sig = diff$p < 0.05
diff = diff[order(diff$module), ]

# write output 
outFile = paste("data/expression/DE", expt, "regulation.tsv", 
                sep = "/")
write.table(diff, outFile, sep = "\t", quote = F, row.names = F, col.names = T)

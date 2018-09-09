# Analyze module regulation in RNA-seq data from GSE76097.
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(limma)

# load preprocessed expression data
load("data/expression/geo/GSE76097/GSE76097-human.Rdata")
# read experiment targets
targets = read.delim("data/expression/geo/GSE76097/GSE76097-targets.txt")

# remove control rats
exprEnsembl = exprEnsembl[, targets$condition != "Uninjured"]

# create design matrix
f1 = factor(targets$genotype[targets$condition != "Uninjured"])
d1 = model.matrix(~ 0 + f1)
colnames(d1) = levels(f1)
cm1 = makeContrasts(STAT3_KO = "WT-STAT3", levels = d1)

# fit the linear model
fit = lmFit(exprEnsembl, d1)
fit2 = contrasts.fit(fit, cm1)
fit2 = eBayes(fit2)

# load GTEx modules
modules = read.delim("data/modules/GTEx-modules.tsv")

# finally, analyze regulation of each module
diff = data.frame(module = integer(), condition = character(), p = numeric(),
                  direction = character())
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
diff$p = diff$p * 2 * max(modules$module) # Bonferroni correction
diff$sig = diff$p < 0.05
diff = diff[order(diff$module), ]

# write output
write.table(diff, "data/expression/geo/GSE76097/regulation.tsv", sep = "\t", 
            quote = F, row.names = F, col.names = T)

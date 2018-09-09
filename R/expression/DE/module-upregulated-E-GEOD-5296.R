# Test whether M3 and M7 genes are differentially expressed following spinal 
# cord injury in the unpublished dataset at E-GEOD-5296 using limma.
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(limma)
source("R/expression/DE/DE-functions.R")

# read modules
modules = read.delim("data/modules/GTEx-modules.tsv")

# load preprocessed expression data
array = "A-AFFY-45"
expt = "E-GEOD-5296" # unpublished dataset
load(paste0("data/expression/DE/", expt, "/", expt, "-human.Rdata"))

# read experiment targets
expt.dir = paste("data/expression/arrayexpress/cel", array, expt, sep = "/")
targets = readTargets(paste("data/expression/DE", expt, "targets.tsv", 
                             sep = "/"))
targets$location[targets$location == "at site"] = "site"
exprEnsembl = exprEnsembl[, colnames(exprEnsembl) %in% targets$filename]

# create design matrices
f1 = factor(with(targets, paste0(severity, location))) # injury vs. control
f2 = factor(with(targets, paste0(severity, location, time))) # as above + time
d1 = model.matrix(~0 + f1)
colnames(d1) = levels(f1)
d2 = model.matrix(~0 + f2)
colnames(d2) = levels(f2)

# define a contrast matrix
times = unique(targets$time)
levels = unique(targets$location)
cm1 = makeContrasts(injury="moderatesite-shamsite", levels = d1)
args2 = paste0("moderatesite", times, "-shamsite", times)
names(args2) = paste0("E", times)
args2 = as.list(args2)
args2$levels = d2
cm2 = do.call(makeContrasts, args2)

# analyze injury vs. control, plus injury vs. control at timepoints
matrix.list = list(cm1, cm2)
designs = list(d1, d2)
diff = data.frame(module = integer(), condition = character(), p = numeric(), 
                   direction = character())
for (i in 1:length(matrix.list)) {
  # define design
  design = designs[[i]]
  # fit the linear model
  fit = lmFit(exprEnsembl, design)
  # define matrix
  cont.matrix = matrix.list[[i]]
  # extract linear model fit for contrasts
  fit2 = contrasts.fit(fit, cont.matrix)
  fit2 = eBayes(fit2)
  # finally, do geneSetTest for each module
  for (module in unique(modules$module)) {
    if (module == 0)
    next 
    moduleGenes = modules$gene[modules$module == module]
    moduleOnChip = intersect(moduleGenes, rownames(exprEnsembl))
    conditions = colnames(fit2$t)
    for (i in seq_len(length(conditions))) {
      condition = conditions[i]
      for (alt in c("up", "down")) {
        test = geneSetTest(moduleOnChip, fit2$t[, i], type="t", 
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
outFile = paste("data/expression/DE", expt, "regulation.tsv", sep = "/")
write.table(diff, outFile, sep="\t", quote=F, row.names=F)

# Evaluate spinal cord gene coexpression module preservation in our proteomics
# dataset.
setwd("~/git/spinal-cord-injury-elife-2018")
options(stringsAsFactors = F)
library(WGCNA)

# load GTEx expression data
load("data/modules/GTEx-expression.Rdata")
gtex = expr
# load GTEx modules
modules = read.table("data/modules/GTEx-modules.tsv", header = T)
modules = modules[, 2]

# read preprocessed proteomics data
expr = read.delim("data/proteomics/proteomics-processed.txt")
genes = expr$Ortholog
expr = expr[, -1]
rownames(expr) = genes
# transpose for WGCNA and log transform
expr = t(expr)

# set up module preservation
setLabels = c("GTEx", "Proteomics")
multiExpr = list(GTEx = list(data = gtex), Proteomics = list(data = expr))
multiColor = list(GTEx = modules)
# run module preservation
mp = modulePreservation(multiExpr, multiColor, referenceNetworks = 1, 
                        nPermutations = 100, randomSeed = 1, quickCor = 0, 
                        verbose = 3, networkType = "signed", 
                        corOptions = "use = 'p', method = 'spearman'")
# produce Z-summary table
ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1],
                 mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1],
               mp$preservation$Z[[ref]][[test]][, -1])
table = cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
              signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2))
table$module = rownames(table)
table = table[, c("module", "Zsummary.pres", "medianRank.pres",
                  "medianRank.qual", "Zsummary.qual")] # reorder
# Write results data frame
write.table(table, "data/proteomics/preservation.txt", sep = "\t",
            quote = F, row.names = F)
